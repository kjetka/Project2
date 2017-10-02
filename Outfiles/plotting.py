from __future__ import division
from numpy import * 
from numpy.linalg import *
from matplotlib.pyplot import*
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


figure()
plot(n, iterations, 'r*-', label ='Calculated values')

plot(x_plot,y_fit, label ='Approximated function')
xlabel('Mesh points N')
ylabel('No. Similarity transforms')
legend()
savefig('fit.pdf')
show()

figure()

plot(log10(n),log10(iterations), label = 'Time arma')
legend()
show()

"""

for k in range(len(filer)):

	nummerisk = []
	analytisk = []
	feil =[]
	nummerisk.append(0)
	analytisk.append(0)

	with open(filer[k]) as infile:
		for i in range(3):
			firstline = infile.readline()
		
		for line in infile:
			thisline = line.split()
			nummerisk.append(thisline[0])
			analytisk.append(thisline[1])
			feil.append(thisline[2])
	nummerisk.append(0)
	analytisk.append(0)
	n = len(nummerisk)-2
	print n

	x = linspace(0,1,n+2)
	figure(k)
	plot(x,nummerisk, label = "Numerical")
	plot(x, analytisk, label = "Analytical")
	legend()
	title("Plot with the special algorithm, h = %.0e" %(1/float(n)))
	xlabel("x"); ylabel("y")
	savefig("n_%.0f.pdf" %n)
"""

