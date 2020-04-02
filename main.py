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
	
    return c1*np.float64(x)+c2

#General Functions

def chisq (ph, phfit, nprm=2,sd=None):
	"""Gives either the reduced chi squared of phfit with respect to ph if sd is given or normal chi squared if sd is none."""

	if sd==None:
		chisq = np.sum(np.power(ph-phfit,2))  
	else:
		chisq = np.sum(np.power(ph-phfit, 2)/np.power(sd,2))  
	return chisq 

def pvalue (rcsq, deg):
    """Returns p-value when given a reduced chi square statistic and degrees of freedom"""    		

    if(rcsq>1):
        pv = chi2.sf(deg*rcsq, deg)
    else:
        pv = chi2.cdf(deg*rcsq, deg)
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
        return avgphs, stdphs, avgpis, stdpis

	#Gives parameters for fit using func. Takes the function used for the fit, a list of pulse intregrals, a list of pulse heights and a list of uncertainties in the pulse heights.
	#Returns the parameters of the fit and the chi square with p-value of the fit.
    def giveFitParam(self, func, pi, ph, sph):	
        paramc = curve_fit(func, pi, ph, sigma=sph, method='lm')
        param = paramc[0]
        csq = chisq(ph, func(pi, *param), sd=sph)
        return param, csq
	
    def findChiSqrs(self, bins = (np.linspace(0,150000, num=120)).tolist(), startPoint = 2):
        phs, sphs, pis, spis = self.avgs(bins = bins)
        cpis = []
        cphs = []
        csqs = []
        pvs = []
        for j, ph in enumerate(phs):
            i = 4
            i += 1
            cpi = []
            cph = []
            csq = []
            pv = []
            if(ph and pis[j]):
                while i <= len(ph):
                    prm, c = self.giveFitParam(linfit, pis[j][startPoint:(startPoint + i)], ph[startPoint:(startPoint + i)], sphs[j][startPoint:(startPoint + i)])
                    cpi.append(pis[j][i-1])
                    cph.append(ph[i-1])
                    csq.append(c)
                    pv.append(pvalue(c, len(ph[startPoint:(startPoint + i)])-2))
                    i += 1
            else:
                print('Could not find saturation point: at least one data set is empty!')
            cpis.append(cpi)
            cpi = []
            cphs.append(cph)
            cph = []
            csqs.append(csq)
            csq = []
            pvs.append(pv)
            pv = []
        return cpis, cphs, csqs, pvs
        
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

#Testing Area
x = saturation(stations[0])
pi, spi, ph, sph = x.avgs()
cpi, cph, csq, pv = x.findChiSqrs()









































