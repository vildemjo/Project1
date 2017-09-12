# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 13:49:32 2017

@author: vilde
"""

from __future__ import division
from numpy import * 
from numpy.linalg import *
from matplotlib.pyplot import*
import glob, os


def filer_funk(): #Stikkord = hvilke filer, Skanntype = n
	skanntype = "error"
	filer = []
	filnavn = []
	for file in glob.glob("error.txt"):
		filer.append(file)
	antall_filer = len(filer)
	
	
	# filer_med_stikkord =[]
	# legendnavn_med_stikkord = []
	# for fil in filer:
	# 	if stikkord in fil:
	# 		filer_med_stikkord.append(fil)
	# 		legendnavn = fil[fil.find(skanntype):fil.find(skanntype)+4]
	# 		bytte = legendnavn.find('_')
	# 		legendnavn =  legendnavn[0:bytte]+ '='+ legendnavn[bytte+1:]
	# 		legendnavn_med_stikkord.append(legendnavn)
	# return filer_med_stikkord, legendnavn_med_stikkord	#, plotnavn, antall_filer 		#filnavn = filer - .txt

	return filer

print filer_funk()
filer = filer_funk()

for k in range(len(filer)):

	h = []
	error = []
	with open(filer[k]) as infile:
		for i in range(1):
			firstline = infile.readline()
		
		for line in infile:
			thisline = line.split()
			h.append(thisline[0])
			error.append(thisline[1])

	n = len(nummerisk)-2
	print n

	figure(k)
	plot(h, error)
	legend()
	title("Plot of the developement of the error")
	xlabel("log(h)"); ylabel("log(Error)")
	savefig("ErrorDevelopement.pdf")
