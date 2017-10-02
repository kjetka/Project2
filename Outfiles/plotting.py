from __future__ import division
from numpy import * 
from numpy.linalg import *
from matplotlib.pyplot import*
import matplotlib.ticker
from scipy.optimize import curve_fit

import glob, os

def func(x,a0,e0,x0):
   return (e0+(a0*(x-x0)**2))

def filer_funk(stikkord): #Stikkord = hvilke filer, Skanntype = n

	filer = []
	filnavn = []
	for file in glob.glob("*.txt"):
		filer.append(file)
	#antall_filer = len(filer)
	filer_brukes = []
	for i in range(len(filer)):
		if stikkord in filer[i]:
			filer_brukes.append(filer[i])
	
	return filer_brukes

filer = filer_funk("time")
fil = filer[0]
n = []
iterations =[]
with open(fil) as infile:
	for i in range(2):
		firstline = infile.readline()
	for line in infile:
		thisline = line.split('&')
		n.append(float(thisline[0]))
		iterations.append(float(thisline[1]))

p0=[30,-12,2.8]
popt,pcov = curve_fit(func,n,iterations,p0)
a0 = popt[2]
beta = popt[1]
B0 = 2/9*(1/a0)*beta

x_plot = np.linspace(n[0],n[-1],100)
y_fit = popt[1]+popt[0]*(x_plot-popt[2])**2
y_fit = popt[1]+popt[0]*(x_plot)**2


figure()

plot(n, iterations, 'r*-', label ='Calculated values')
plot(x_plot,y_fit, label ='Approximated function')
ticklabel_format(style='sci', axis='y', scilimits=(0,0))

xlabel('Mesh points N')
ylabel('No. Similarity transforms')
legend(loc = 2)
title('Increase in Similarity transforms')
savefig('fit.pdf')



print "a0 = ",popt[0], 'x0 = ', popt[2]