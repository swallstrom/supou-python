'''
A python version of Brandon Kellys' IDL script supou.pro, as described in Kelly et al. 2011

FIT A SUPERPOSITION OF ORNSTEIN-UHLENBECK (sup-OU) PROCESSES TO A
TIME SERIES SAMPLED AT IRREGULARY INTERVALS WITH MEASUREMENT
ERROR. BAYESIAN INFERENCE IS PERFORMED. THE OU PROCESSES ARE ASSUMED
TO HAVE THE SAME MEAN AND WHITE NOISE AMPLITUDE (SIGMA), BUT
DIFFERENT VALUE OF THE ANGULAR BREAK FREQUENCY, WHICH FALL ON A
REGULAR GRID FROM OMEGA_1 TO OMEGA_2.
'''

import numpy as np
from scipy.optimize import minimize
import emcee
import corner
import matplotlib.pyplot as plt
import time
import glob

""" 
*** Remember to change rootdir and filenamelist at end of script! ***
"""

def logit(x,inverse=False):
    '''
    Compute the logit or inverse logit function of x. X must be between 0 and 1.
    '''
    if inverse == False:
        if np.any(x < 0.) or np.any(x > 1.):
            print('All elements of x must be between 0 and 1')
            return 0
    
        logit = np.log(x / (1 - x))
    
        return logit

    else:
        inv_logit = np.exp(x) / (1 + np.exp(x))
        
        return inv_logit


def kalman_params_supou(y, time, yvar, thetai, nou, counts):
    '''
    Compute the mean and variance of the joint likelihood function for (y - mu) using the Kalman recursions
    
    Inputs:
    y - timeseries fluxes
    time - timeseries days since first observation
    yvar - timeseries flux error
    thetai - starting parameters [omega1, omega2, mu, ysigma, slope]. Omega are the break frequencies of the sup-OU process. Mu is the mean level of the timeseries, ysigma the variance. Slope is the slope between omega1 and omega2 (between -2 and 0)
    nou - number of OU processes in the sum used to model y
    counts - boolean
    '''
    ny = len(y)
    omega1 = np.exp(thetai[0])
    omega2 = np.exp(thetai[1])
    mu = thetai[2]
    ysigma = thetai[3]
    slope = 2. * logit(thetai[4], inverse=True) 
    
    omgrid = np.arange(nou) / (nou - 1.) * (thetai[1] - thetai[0]) + thetai[0]
    omgrid = np.exp(omgrid)
    
    weights = omgrid**(1. - slope / 2.) / np.sqrt(np.sum(omgrid**(2. - slope)))
        
    sigsqr = 2. * ysigma**2 / np.sum(weights**2 / omgrid)
    
    yhat = np.zeros(ny)
    yhvar = np.zeros(ny)
    yhat[0] = mu
    if counts:
        yhvar[0] = ysigma**2 + 1. / np.exp(yhat[0])
    else:
        yhvar[0] = ysigma**2 + yvar[0]
    
    xhat = np.zeros(nou)
    xcovar = np.diag( sigsqr / (2. * omgrid) )
    
    for i in range(1,ny):
    
        dt = time[i] - time[i-1]
        alpha = np.exp(-1. * omgrid * dt)
        
        ouvar = sigsqr / (2. * omgrid) * (1. - alpha**2)
        
        psi = alpha * np.matmul(weights, xcovar)
                       
        xhat = alpha * xhat + psi / yhvar[i-1] * (y[i-1] - mu - np.sum(weights * xhat))
                       
        yhat[i] = mu + np.sum(weights * xhat)
                       
        xcovar = np.outer(alpha,alpha) * xcovar + np.diag(ouvar) - np.outer(psi, psi) / yhvar[i-1]
        if counts:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + 1. / np.exp(yhat[i])
        else:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + yvar[i]
    
    return yhat, yhvar


def lnlike(y, yhat, yhvar, counts=False):
    '''Compute log-likelihood of the data'''
    lnlike = -0.5 * np.log( 2. * np.pi * yhvar ) - 0.5 * (y - yhat)**2 / yhvar
    return np.sum(lnlike)

def lnprior(thetai, omega_max,omega_min, counts=False):
    '''
    Compute logarithm of the prior
    
    omega1 and omega2 must be between omega_min and omega_max
    mu must be between -100 and 100
    ysigma must be > 0
    '''
    if np.log(omega_min) < thetai[0] < np.log(omega_max) and np.log(omega_min) < thetai[1] < np.log(omega_max) and thetai[0] < thetai[1] and -100 < thetai[2] < 100 and thetai[3] > 0:
        return -1. * np.log(np.log(omega_max) - thetai[0]) + np.log( logit(thetai[4], inverse=True) * (1. - logit(thetai[4], inverse=True)) )
    return -np.inf

def lnprob(thetai,y,time,yvar,nou,omega_max, omega_min, counts=False):
    '''Compute logarithm of the posterior'''
    yhat,yhvar = kalman_params_supou(y, time, yvar, thetai, nou, counts) 
    #print(yhat,yhvar)
    lp = lnprior(thetai, omega_max,omega_min)
    if lp == -np.inf:
        return lp
    if np.isnan(lp):
        print(thetai,yhat,yhvar,lp,np.log(np.log(omega_max) - thetai[0]),logit(thetai[4], inverse=True),(1. - logit(thetai[4], inverse=True)))
        exit()
    print(lp)
    ll = lnlike(y, yhat, yhvar)
    if np.isnan(ll):
        print(thetai,yhat,yhvar,ll,(y - yhat)**2)
        exit()
    print(ll)
    print(lp + ll)
    return lp + ll

def minlnprob(thetai,y,time,yvar,nou,omega_max,omega_min, counts=False):
    return -1*lnprob(thetai,y,time,yvar,nou,omega_max,omega_min)

def supou(y, time, yvar, name, rootdir, nou=32, miniter=20000, maxiter=50000,
          burniter=100, nwalkers=50, silent=False, counts=False, mle=False,
          **kwargs
        ):

    ny = len(y)
    nt = len(time)
    nyvar = len(yvar)

    if (ny != nt) or (nt != nyvar) or (ny != nyvar):
        print("time, y and yvar must all be the same size")
        return

    dt = time[1:] - time[0:-1]
    omega_min = 1. / (10. * (time[-1] - time[0]))
    omega_max = 1. / min(dt / 10.)

    ## First get initial guess values
    mu = np.mean(y) + np.std(y) / np.sqrt(ny)
    ysig0 = (np.var(y) - np.median(yvar))
    if ysig0 < 0: ysig0 = 0.
    ysig0 = np.sqrt(ysig0)
    if ysig0 < 0.1 * np.std(y):  ysig0 = 0.1 * np.std(y)
    ysig = ysig0
    print(ysig)

    omega1 = np.log10(omega_max / omega_min)*0.5 + np.log10(omega_min) -1 
    omega1 = 10**omega1
    omega2 = np.log10(omega_max / omega_min)*0.5 + np.log10(omega_min)
    omega2 = 10.**omega2

    slope = 0.9
    
    ## Get maximum-likelihood estimate for initial guess
    theta = np.array([[np.log(omega1)], [np.log(omega2)], [mu], [ysig], [logit(slope / 2.)]]).T

    print(theta)
    
    if not silent:
        print('Getting MLE for initial guesses...')

    res = minimize(minlnprob,
                   theta,
                   args = (y,time,yvar,nou,omega_max,omega_min),
                   bounds=[(np.log(2.*omega_min),np.log(omega_max/2.)),
                           (np.log(2.*omega_min),np.log(omega_max/2.)),
                           (np.min(y), np.max(y)),
                           (1e-4*np.std(y), 1e300),
                           (-2,2)]
                   )

    print(res.message)
    print(res.x)
    print(res.success)

    sampler=emcee.EnsembleSampler(nwalkers,5,lnprob,args=(y,time,yvar,nou,omega_max,omega_min))
    pos = [res.x * 1e-3*np.random.randn(5) for i in range(nwalkers)]

    ## Do a burn-in first, then reset and start the proper run from the burn-in position
    pos_burn, prob_burn, state_burn = sampler.run_mcmc(pos,burniter)
    sampler.reset()
    sampler.run_mcmc(pos_burn, 1000)

    samples=sampler.chain.reshape((-1,5)) 
    omega1_mcmc, omega2_mcmc, mu_mcmc, ysig_mcmc, slope_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))


    ### NOTE: want exp(omega1) and exp(omega2) to later take the log10 when plotting
    with open("{}/{}-python.sup".format(rootdir, name), 'w') as f:
        f.write('#omega1,omega2,mu,ysig,slope\n')
        for i,v in enumerate(samples):
            f.write('{},{},{},{},{}\n'.format(np.exp(samples[i][0]), np.exp(samples[i][1]), samples[i][2], samples[i][3], samples[i][4],))
    f.closed


    print(omega1_mcmc, omega2_mcmc, mu_mcmc, ysig_mcmc, slope_mcmc)

    ## Plot the walker paths
    fig, axes = plt.subplots(5, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)    
    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
    axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
    axes[4].plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)

    axes[0].set_ylabel(r"$\Omega 1$")
    axes[1].set_ylabel(r"$\Omega 2$")
    axes[2].set_ylabel(r"$\mu$")
    axes[3].set_ylabel(r"$\sigma_y$")
    axes[4].set_ylabel(r"$slope$")

    fig.tight_layout(h_pad=0.0)
    fig.savefig("{}/{}_line-time.pdf".format(rootdir, name))
    plt.clf()
    fig2=corner.corner(samples,labels=[r"$\log \Omega_1$", r"$\log \Omega_2$", r"$\mu$", r"$\sigma_y$", r"$\mathrm{slope}$"])
    #plt.show()
    fig2.savefig("{}/{}_triangle.pdf".format(rootdir, name))
    return 



if __name__ == "__main__":
    #Time the script
    startTime = time.time()

    # Directory of input data and output results
    rootdir = '/Users/sofiaw/data/JCMT/SCUBA2/QuasarVariability/results'

    # Input files = light curves of the form: days since first obs, flux, flux uncertainty
    filenamelist = glob.glob('{}/*_allfluxes_supou.csv'.format(rootdir))

    for filename in filenamelist:
        day, y, yvar = np.loadtxt(filename, delimiter=', ', unpack=True)
        supou(y, day, yvar, filename.split('/')[-1].split('_')[0], rootdir)


    print(" -----\n Script took {:.1f} seconds ({:.1f} minutes) \n -----".format(time.time()-startTime, (time.time()-startTime)/60.))

