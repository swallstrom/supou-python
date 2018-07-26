
import numpy as np

""" TO DO:
Check: diag_matrix - I think that's the same as np.diag
logit/inv_logit

"""
def kalman_params_supou( y, time, yvar, theta, nou, counts): #, yhat, yhvar

    ny = len(y)
    omega1 = np.exp(theta[0])
    omega2 = np.exp(theta[1])
    mu = theta[2]
    ysigma = theta[3]
    slope = 2. * inv_logit(theta[4]) #???
    #;slope = 1d
    
    omgrid = np.linspace(nou) / (nou - 1.) * (theta[1] - theta[0]) + theta[0]
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
        yhvar[0] = ysigma^2 + yvar[0]
    
    xhat = np.empty(nou)
    xcovar = np.diag( sigsqr / (2. * omgrid) )
    
    for i in range(1,ny-1):#= 1, ny - 1:
    
        dt = time[i] - time[i-1]
        alpha = np.exp(-1. * omgrid * dt)
        
        ouvar = sigsqr / (2. * omgrid) * (1. - alpha**2)
        
        psi = alpha * np.matmul(weights, xcovar)
                       
        xhat = alpha * xhat + psi / yhvar[i-1] * (y[i-1] - mu - np.sum(weights * xhat))
                       
        yhat[i] = mu + total(weights * xhat)
                       
        xcovar = np.matmul(alpha*alpha) * xcovar + np.diag(ouvar) - np.matmul(psi, psi) / yhvar[i-1]
        if counts:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + 1. / np.exp(yhat[i])
        else:
            yhvar[i] = np.matmul(weights,np.matmul(xcovar,weights)) + yvar[i]
    
        #endfor

    return yhat, yhvar#, xcovar #What is this supposed to be returning??? 
#end

def lnlike(theta, y, yhat, yhvar):
    
    lnlike = -0.5 * np.log( 2. * np.pi * yhvar ) - 0.5 * (y - yhat)^2 / yhvar
    return np.sum(lnlike)

def lnprior(theta, omega_max):
    return -1. * np.log(np.log(omega_max) - theta[0]) + np.log( inv_logit(theta[4]) * (1. - inv_logit(theta[4])) )

def lnprob(theta,y,time,yvar,nou,omega_max):
    yhat,yhvar = kalman_params_supou(y, time, yvar, theta, nou, counts) #, yhat, yhvar
    return lnprior(theta, omega_max) * lnlike(theta, y, yhat, yhvar)
