#Import Packages

import numpy as np

import pylab as pl

import matplotlib as mpl

import matplotlib.dates as d

import datetime as dt

import scipy.stats

import tkinter.filedialog as tk

from scipy.stats import chi2

from math import exp

from scipy.stats import norm, rayleigh, expon, moyal, chisquare

from scipy.optimize import curve_fit

#Initialize Vectors
stations = []

#Function used to fit the saturation plot

def linfit (x, c1, c2):
	
    return c1*np.float64(x)

def tanhfit (x, c1, c2):
    
    return -0.5*np.tanh(c1*np.float64(x)+c2)+0.5

#General Functions
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    val = array[idx]
    return idx, val

def chisq (ph, phfit, nprm=1,sd=None):
	"""Gives either the reduced chi squared of phfit with respect to ph if sd is given or normal chi squared if sd is none."""

	if sd==None:
		chisq = np.sum(np.power(ph-phfit,2))  
	else:
		chisq = np.sum(np.power(ph-phfit, 2)/np.power(sd,2))  
	return chisq 

def pvalue (rcsq, deg):
    """Returns p-value when given a reduced chi square statistic and degrees of freedom"""    		

    pv = chi2.sf(rcsq, deg)

    return pv
    
#Declare Classes

class saturation:
    """A class used for finding the satuartion point of a station."""
    #Initialization of Class	
    def __init__(self, station):
            self.phs = [station[5], station[6], station[7], station[8]]
            self.pis = [station[9], station[10], station[11], station[12]]
            self.name = station[0]

	#given a set of bin edges returns the average pulse height and the average pulse integral of pulse heights between given edges. Empty bins are discarded. 
	#Also returns the standard deviations in both the pulse intregrals and pulse heights
    def avgs(self, bins = (np.linspace(0,150000, num=120)).tolist()):
        pihs = []
        for i, ph in enumerate(self.phs):
            pih = list(map(lambda x, y:(x,y), self.pis[i], ph))
            pih.sort(key=lambda tup: tup[0])
            pihs.append(pih.copy())
        
        avgphs = []
        stdphs = []
        avgpis = []
        stdpis = []
        
        for pih in pihs:
            temp = [[],[]]
            avgph = []
            stdph = []
            avgpi = []
            stdpi = []
            counter = 1
            for pi in pih:
                while counter<len(bins):
                    if pi[0] < bins[counter]:					
                        temp[0].append(pi[0])					
                        temp[1].append(pi[1])
                        break
                    elif temp[0] and temp[1]:
                        counter += 1
                        """avgpi.append(np.average(temp[0]))
                        stdpi.append(np.std(temp[0]))
                        avgph.append(np.average(temp[1]))
                        stdph.append(np.sqrt(np.sum(np.power(np.sqrt(temp[1]),2))))
                        temp = [[],[]]
                        continue"""
                        if len(temp[0]) > 5:
                            avgpi.append(np.average(temp[0]))
                            stdpi.append(np.std(temp[0]))
                            avgph.append(np.average(temp[1]))
                            stdph.append(np.std(temp[1]))
                            temp = [[],[]]
                            continue
                        else:
                            avgpi.append(temp[0][0])
                            stdpi.append(1)
                            avgph.append(temp[1][0])
                            stdph.append(np.sqrt(np.average(temp[1])))
                            temp = [[],[]]
                            continue				
                    else:
                        counter += 1
                        continue
            avgphs.append(avgph)
            avgph = []
            stdphs.append(stdph)
            stdph = []
            avgpis.append(avgpi)
            avgpi = []
            stdpis.append(stdpi)
            stdpi = []
        return avgpis, stdpis, avgphs, stdphs

	#Gives parameters for fit using func. Takes the function used for the fit, a list of pulse intregrals, a list of pulse heights and a list of uncertainties in the pulse heights.
	#Returns the parameters of the fit and the chi square with p-value of the fit.
    def giveFitParam(self, func, x, y, yerr):	
        paramc = curve_fit(func, x, y, sigma=yerr, method='lm')
        param = paramc[0]
        csq = chisq(y, func(x, *param), sd=yerr)
        return param, csq
	
    def findChiSqrs(self, bins = (np.linspace(0,150000, num=120)).tolist(), startPoint = 2):
        pis, spis, phs, sphs = self.avgs(bins = bins)
        cpis = []
        cphs = []
        csqs = []
        pvs = []
        prms = []
        for j, ph in enumerate(phs):
            i = 4
            i += 1
            cpi = []
            cph = []
            csq = []
            pv = []
            prm = []
            if(ph and pis[j]):
                while i <= len(ph):
                    parm, c = self.giveFitParam(linfit, pis[j][startPoint:(startPoint + i)], ph[startPoint:(startPoint + i)], sphs[j][startPoint:(startPoint + i)])
                    cpi.append(pis[j][i-1])
                    cph.append(ph[i-1])
                    csq.append(c)
                    pv.append(pvalue(c, len(ph[startPoint:(startPoint + i)])-1))
                    prm.append(parm)
                    i += 1
                    
            cpis.append(cpi)
            cpi = []
            cphs.append(cph)
            cph = []
            csqs.append(csq)
            csq = []
            pvs.append(pv)
            pv = []
            prms.append(prm)
            prm = []
        return cpis, cphs, csqs, pvs, prms
    
    def findSatPoint(self):
        """Finds saturation point for this station station"""
        cpis, cphs, csqs, pvs, prms = self.findChiSqrs()
        
        prms = []
        npvs = []
        perrs = []
        satpts = []
        saterrs = []
        x =  np.linspace(0, 4500, 250)
        for i, pv in enumerate(pvs):
            if pv:
                #Make fit to pvalue data and calculate pvalue of fit
                param, cov = curve_fit(tanhfit, cphs[i], pv, p0=[0.005, -10], sigma=0.01*np.ones(len(pv)))
                perr = np.sqrt(np.diag(cov))
                perrs.append(perr)
                prms.append(param)
                ncsq = chisq(pv, tanhfit(cphs[i], *param))
                npv = pvalue(ncsq, len(pv)-2)
                npvs.append(npv)
                
                #Find saturation point
                idx, val =  find_nearest(tanhfit(x, *param), 0.9)
                satpt = x[idx]
                satpts.append(satpt)
                
                #Find saturation point error
                param_plus = param + perr
                param_min = param - perr
                idx_plus, val_plus =  find_nearest(tanhfit(x, *param_plus), 0.9)
                saterr_plus = np.abs(x[idx] - x[idx_plus])
                idx_min, val_min =  find_nearest(tanhfit(x, *param_min), 0.9)
                saterr_min = np.abs(x[idx] - x[idx_min])
                saterr = [saterr_min, saterr_plus]
                saterrs.append(saterr)
                
                """
                #Make plot of p-values and fit
                pl.subplot(1,2,i+1)
                pl.plot(x, tanhfit(x, *param))
                pl.plot(x, tanhfit(x, *param_min), '--k')
                pl.plot(x, tanhfit(x, *param_plus), '--k')
                pl.plot(cphs[i], pv, '.')
                pl.xlabel('Pulse Height')
                pl.ylabel('P-Value')
                pl.hlines(0.9, 0, 4500, 'r')
                pl.xlim(0,4500)
                pl.ylim(-.1,1.1)
                
        pl.show()
        """
        return satpts, saterrs, npvs
        
        
#Filling Vectors with Data

#Loop Over the Data Files
root = tk.Tk()
root.withdraw()
files = tk.askopenfilenames(parent=root, title='Choose any files')
root.destroy()

for file in files:
    
    with open(file) as fa:

	#Read all lines
        station = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
        for line_aa in fa.readlines():

            if line_aa.startswith("# Station:"):
                a = line_aa.split(':',2)
                station[0].append(a[1])
            elif len(line_aa.split('\t',23))==23:
                a = line_aa.strip()

		#read each single column

                cols = line_aa.split('\t',23)
		
		#Now you have all the colums and you can do all your analysis. Col1 through to Col4 countain data about the time of events. Col5 through to Col8 contain pulseheights of events. Col9 			through to 12 contain pulse integrals. Col 13 through to 16 contain the estimate for the number of particles detected. Col17 through to 20 contain the relative time between detections. Col 			21 contains trigger times. Col 22 and 23 the zenith and azimuth of the shower (currently not useable).
                for i, col in enumerate(cols):
                    if(i>3):
                        station[i+1].append(float(col))
                    else:
                        station[i+1].append(col)
            else:
                continue
        stations.append(station.copy())
        station.clear()

"""

#Create comparison plot of all stations, assuming data given was from different stations.
        
#First create a list of all instances of the saturation class
classes = []
for station in stations:
    classes.append(saturation(station))

#Create list of saturation points and their errors 
satpts = []
saterrs = [[],[]]
for cl in classes:
    satpt, saterr, pv = cl.findSatPoint()
    satpts.append(satpt)
    for sater in saterr:
        saterrs[0].append(sater[0])
        saterrs[1].append(sater[1])
    
satpts_flat = [val for sublist in satpts for val in sublist]

#Create list of station names
names = []
for i, satpt in enumerate(satpts):
    for j, pt in enumerate(satpt):
        names.append(stations[i][0][0] + str( j + 1))


#Create actual plot
x = np.arange(len(satpts_flat))
pl.errorbar(x, satpts_flat, yerr=saterrs, fmt='none', capsize=5)
pl.plot(x, satpts_flat, '.r')
pl.xticks(x, names)
pl.xlabel('Stations')
pl.ylabel('Saturation Point (Pulse Height)')
pl.show()

"""

#Create comparison plot of station, assuming data given was from the same station at different times.

#First create a list of all instances of the saturation class
classes = []
stations.sort(key=lambda x: x[1])
for station in stations:
    classes.append(saturation(station))

#Create list of saturation points and their errors 
satpts = []
saterrs = [[],[]]
for cl in classes:
    satpt, saterr, pv = cl.findSatPoint()
    satpts.append(satpt)
    for sater in saterr:
        saterrs[0].append(sater[0])
        saterrs[1].append(sater[1])
    
satpts_flat = [val for sublist in satpts for val in sublist]

#Create list of station names
names = []
for i, satpt in enumerate(satpts):
    for j, pt in enumerate(satpt):
        names.append(stations[i][1][0])

#Create actual plot
for i, satpt in enumerate(satpts[0]):
    pts = satpts_flat[i::len(satpts[0])]
    errs = [saterrs[0][i::len(satpts[0])], saterrs[1][i::len(satpts[0])]]
    x = np.arange(len(pts))
    pl.subplot(len(satpts[0])/2, 2, i+1)
    pl.errorbar(x, pts, yerr=errs, fmt='none', capsize=5)
    pl.plot(x, pts, '.r')
    pl.xticks(x, names[i::len(satpts[0])], rotation='vertical')
    pl.xlabel('Dates')
    pl.ylabel('Saturation Point (Pulse Height)')
    pl.title('Detector Nr.:' + str(i+1))
pl.show()

"""

#Testing Area
x = saturation(stations[0])
pi, spi, ph, sph = x.avgs()
cpi, cph, csq, pv, prm = x.findChiSqrs()
satpts, saterrs, npvs = x.findSatPoint()

print(satpts)


pl.subplot(221)
pl.errorbar(pi[0],ph[0], yerr=sph[0], fmt='none')
pl.plot(pi[0],ph[0], '.r')
pl.xlabel('Pulse Integral')
pl.ylabel('Pulse Height')
pl.xlim(0,150000)
pl.ylim(0,4500)
pl.plot(np.linspace(0,150000,200),linfit(np.linspace(0,150000,200), *prm[0][find_nearest(pi[0],20000)[0]]))

pl.subplot(222)
pl.plot(stations[0][9],stations[0][5], '.r')
pl.xlabel('Pulse Integral')
pl.ylabel('Pulse Height')
pl.xlim(0,150000)
pl.ylim(0,4500)
pl.plot(np.linspace(0,150000,200),linfit(np.linspace(0,150000,200), *prm[0][20]))

pl.subplot(223)
pl.plot(cpi[0], csq[0], '.')
pl.xlabel('Pulse Integral')
pl.ylabel('Reduced Chi Sqr')
pl.xlim(0,150000)

pl.subplot(224)
pl.plot(cph[0], pv[0], '.')
pl.xlabel('Pulse Integral')
pl.ylabel('P-Value')
pl.xlim(0,4500)

pl.show()
"""








































