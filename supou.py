
import numpy as np

""" TO DO:
Doublecheck the matrix multiplication!
"""
def logit(x,inverse=False):
    
    if inverse == False:
        xvals == x[(x < 0) or (x > 1)]
        
        if True in xvals:
            print('All elements of x must be between 0 and 1')
            return, 0
    
        logit = np.log(x / (1 - x))
    
        return, logit

    else:
        inv_logit = exp(x) / (1 + exp(x))
        
        return, inv_logit


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
                       
        xcovar = np.outer(alpha*alpha) * xcovar + np.diag(ouvar) - np.outer(psi, psi) / yhvar[i-1]
        if counts:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + 1. / np.exp(yhat[i])
        else:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + yvar[i]
    
    return yhat, yhvar#, xcovar #What is this supposed to be returning??? 


def lnlike(y, yhat, yhvar):
    
    lnlike = -0.5 * np.log( 2. * np.pi * yhvar ) - 0.5 * (y - yhat)**2 / yhvar
    return np.sum(lnlike)

def lnprior(thetai, omega_max):
    return -1. * np.log(np.log(omega_max) - thetai[0]) + np.log( logit(thetai[4], inverse=True) * (1. - logit(thetai[4], inverse=True)) )

def lnprob(thetai,y,time,yvar,nou,omega_max):
    yhat,yhvar = kalman_params_supou(y, time, yvar, thetai, nou, counts) #, yhat, yhvar
    return lnprior(thetai, omega_max) * lnlike(y, yhat, yhvar)
