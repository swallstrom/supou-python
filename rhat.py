import numpy as np

###
###  Calculate Rhat, the potential scale reduction factor
###
###  Author: Brandon C. Kelly, Steward Observatory, Jan 2008
###  Converted from IDL -> Python by J.P. Marshall, ASIAA, Jul 2018
###
###  Inputs: Theta: [niter,nchains,npar] array containing parameter
###                 values.
###  Outputs: Rhat: npar element array containing Rhat values.


def rhat_function(theta):
	ntheta = theta.shape()
	niter = ntheta[0]
	nchains = ntheta[1]
	npar = ntheta[2]

	Bvar = np.array()
	Wvar = np.array()

	for j in range(0, npar):
		thetabarj = np.sum(theta[:,:,j],1) / niter
		thetabar  = np.mean(thetabarj)

		sjsqr = 0.

		for i in range(0,nchains):
			sjsqr = sjsqr + np.sum(theta[:,i,j] - thetabarj[i]**2) / (niter - 1)

        np.insert(Bvar,j,niter / (nchains - 1.0) * np.sum( (thetabarj - thetabar)^2 ))
        np.insert(Wvar,j,sjsqr / nchains)

	varplus = (1 - 1.0 / niter) * Wvar + Bvar / niter
	rhat = np.sqrt(varplus/Wvar)

	return, rhat
