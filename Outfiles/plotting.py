from __future__ import division
from numpy import * 
from numpy.linalg import *
from matplotlib.pyplot import*
import glob, os


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
print filer


n,time_jac, time_arma = loadtxt(filer[0],delimiter = '&', unpack=True, skiprows=2)

time_arma = time_arma*1e-3

semilogy(n,time_arma, label = 'Time arma')
semilogy(n, time_jac, label = 'Time jacobi')
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

