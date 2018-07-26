
import numpy as np
from scipy.optimize import minimize

""" TO DO:
Doublecheck the matrix multiplication!
"""
def logit(x,inverse=False):
    
    if inverse == False:
        #xvals = x[(x < 0) or (x > 1)]
        
        #if True in xvals:
        #    print('All elements of x must be between 0 and 1')
        #    return 0
    
        logit = np.log(x / (1 - x))
    
        return logit

    else:
        inv_logit = np.exp(x) / (1 + np.exp(x))
        
        return inv_logit


def kalman_params_supou( y, time, yvar, thetai, nou, counts): #, yhat, yhvar

    ny = len(y)
    omega1 = np.exp(thetai[0])
    omega2 = np.exp(thetai[1])
    mu = thetai[2]
    ysigma = thetai[3]
    slope = 2. * logit(thetai[4], inverse=True) 
    #;slope = 1d
    
    omgrid = np.arange(nou) / (nou - 1.) * (thetai[1] - thetai[0]) + thetai[0]
    omgrid = np.exp(omgrid)
    
    weights = omgrid**(1. - slope / 2.) / np.sqrt(np.sum(omgrid**(2. - slope)))
    
    #;omgrid = dindgen(nou) / (nou - 1d) * (omega2 - omega1) + omega1
    
    sigsqr = 2. * ysigma**2 / np.sum(weights**2 / omgrid)
    
    yhat = np.empty(ny)
    yhvar = np.empty(ny) 
    yhat[0] = mu
    if counts:
        yhvar[0] = ysigma**2 + 1. / np.exp(yhat[0])
    else:
        yhvar[0] = ysigma**2 + yvar[0]
    
    xhat = np.empty(nou)
    xcovar = np.diag( sigsqr / (2. * omgrid) )
    
    for i in range(1,ny-1):#= 1, ny - 1:
    
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
    
    return yhat, yhvar#, xcovar #What is this supposed to be returning??? 


def lnlike(y, yhat, yhvar, counts=False):
    
    lnlike = -0.5 * np.log( 2. * np.pi * yhvar ) - 0.5 * (y - yhat)**2 / yhvar
    return np.sum(lnlike)

def lnprior(thetai, omega_max, counts=False):
    return -1. * np.log(np.log(omega_max) - thetai[0]) + np.log( logit(thetai[4], inverse=True) * (1. - logit(thetai[4], inverse=True)) )

def lnprob(thetai,y,time,yvar,nou,omega_max, counts=False):
    yhat,yhvar = kalman_params_supou(y, time, yvar, thetai, nou, counts) #, yhat, yhvar
    return lnprior(thetai, omega_max) * lnlike(y, yhat, yhvar)

def minlnprob(thetai,y,time,yvar,nou,omega_max, counts=False):
    return -1*lnprob(thetai,y,time,yvar,nou,omega_max)

def supou(y, time, yvar, nou=32, miniter=20000, maxiter=50000, 
          burniter=10000, nwalkers=100, silent=False, counts=False, mle=False,
          #rhat=rhat, that=that, yhat=yhat, yhvar=yhvar,
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

    mu = np.mean(y) + np.std(y) / np.sqrt(ny)# * randomn(seed, nchains)
    ysig0 = np.sqrt((np.var(y) - np.median(yvar)) > 0.)
    ysig0 = ysig0 > 0.1 * np.std(y)
    ysig = ysig0# * (ny - 1) / randomchi(seed, ny - 1, nchains) )

    omega1 = np.log10(omega_max / omega_min) + np.log10(omega_min) - 1 #* randomu(seed, nchains) 
    omega1 = 10**omega1
    omega2 = np.log10(omega_max / omega_min) + np.log10(omega_min) # * randomu(seed, nchains)
    omega2 = 10.**omega2

    #if omega1 > omega2:
    #    a = np.copy(omega2)
    #    omega2 = np.copy(omega1)
    #    omega1 = np.copy(a)

    slope = 0.9 + 0.2 #*np.random.random
    

    theta = np.array([[np.log(omega1)], [np.log(omega2)], [mu], [ysig], [logit(slope / 2.)]]).T
    
    if not silent:
        print('Getting MLE for initial guesses...')

    res = minimize(minlnprob,
                   theta,
                   args = (y,time,yvar,nou,omega_max),
                   method="SLSQP",
                   bounds=[(np.log(2.*omega_min),np.log(omega_max/2.)),
                           (np.log(2.*omega_min),np.log(omega_max/2.)),
                           (-1e300, 1e300),
                           (1e-4*np.std(y), 1e300),
                           (-2,2)]
                   )#, constraints=[{},{},{},{}]

    if res.success:
        print(res.x)
        
    return #post


time, y, yvar = np.loadtxt('example_lc.txt', delimiter=', ', unpack=True)
supou(y, time, yvar)
