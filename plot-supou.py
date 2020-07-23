#Plot distribution of lower and upper break frequencies

import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit, leastsq
from scipy import stats


def plot_distr(filename, p0_init, col, fit=False, **kwargs):
    plt.clf()
    print("File in: "+filename)
    data = np.loadtxt(filename, delimiter=',', usecols=col)
    sourcename = filename.split('.')[0]
    savename = sourcename + '_supou.pdf'
	
    lowhis, lowbins_edge = np.histogram(np.log10(data[:,0]), bins=100, density=True)
    highhis, highbins_edge = np.histogram(np.log10(data[:,1]), bins=100, density=True)
    

    lowbins = []
    for i in range(len(lowbins_edge)-1):
        lowbins.append((lowbins_edge[i]+lowbins_edge[i+1])/2)
    highbins = []
    for i in range(len(highbins_edge)-1):
        highbins.append((highbins_edge[i]+highbins_edge[i+1])/2)


    #Find max of each distr by gauss fit to max 20ish% of data
    def gauss(x, *p):
        off, A1, mu1, sigma1= p
        return off+A1*np.exp(-(x-mu1)**2/(2*sigma1**2))
    def parab(x, *p):
        a, b, c = p
        return a*x**2+b*x+c

    indx = [i for i,x in enumerate(lowhis) if x > 0.8*max(lowhis)]
    dat_lowhis = []
    dat_lowbins = []
    for each in indx:
        dat_lowhis.append(lowhis[each])
        dat_lowbins.append(lowbins[each])

    indx2 = [i for i,x in enumerate(highhis) if x > 0.4*max(highhis)]
    dat_highhis = []
    dat_highbins = []
    for each in indx2:
        dat_highhis.append(highhis[each])
        dat_highbins.append(highbins[each])


    #Plot histograms
    plt.hist(np.log10(data[:,0]), color='cyan', histtype='step', bins=100, normed=True)
    plt.hist(np.log10(data[:,1]), color='orange', histtype='step', bins=100, normed=True)

    if fit:
        p0 = [-5e3, -1e4, -1e4] #initial guess for parab
        fres_low, var_low = curve_fit(parab,dat_lowbins,dat_lowhis, p0=p0)
        fitx_low = np.linspace(min(dat_lowbins), max(dat_lowbins), 1000) #grid for fit
        ffit_low = parab(fitx_low, *fres_low)
        #plt.plot(fitx_low, ffit_low, 'b--')
        omega_low = fitx_low[np.argmax(ffit_low)]
        taulow = 1/(10**omega_low)

        p0 = p0_init #initial guess for gauss
        fres_high, var_high = curve_fit(gauss,dat_highbins,dat_highhis, p0=p0)
        fitx_high = np.linspace(min(dat_highbins), max(dat_highbins), 1000) #grid for fit
        ffit_high = gauss(fitx_high, *fres_high)
        plt.plot(fitx_high, ffit_high, 'r--')
        tauhigh = 1/(10**fres_high[2])
        tauhigh_err = tauhigh*0.5*fres_high[3]/0.434
        print('Results of gauss fit: ', fres_high)

        plt.figtext(0.7,0.83,'tau$_L$ = {:.0f} d'.format(taulow), color='blue')
        plt.figtext(0.7,0.8,'tau$_H$ = {:.0f} +/- {:.0f} d'.format(tauhigh, tauhigh_err), color='red')

    plt.xlabel('log($\omega$)')
    plt.title(sourcename)
    plt.savefig(savename)
    plt.show()
    
    return


p0 = [0.0, 1.0, -2.5, 0.5]


filename = '0238+166-python.sup'
cols = [0,1]

plot_distr(filename, p0, cols, fit=True)







