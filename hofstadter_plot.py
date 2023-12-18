#!/usr/bin/python3
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import functools
from matplotlib.ticker import MultipleLocator
from scipy.linalg import eigh
import matplotlib as mpl
import time 


mpl.use("GTK3Cairo")


plt.rc("text", usetex=True)
plt.rc("font", family="serif")


#entries from V of non-commutative Torus on off-diagonals
def offdiag(i):
	if i==(size):
		return (0)
	elif i==-1:
		return(size-1)
	else:
		return(i)
		
def M(alpha,centre,u,v):
	M=np.zeros((size,size))
	for i in range(0,size):
		M[i,i]=u*2*np.cos((2*np.pi*(i+1)*alpha)-centre)
		M[i,offdiag(i+1)] = v
		M[i,offdiag(i-1)] = v
	return M

if __name__== "__main__":


	u=1/6
	v=1/6
	w=1/6

	size= 180
	alpharr = np.zeros(size+1)

	for i in range(1,size+1):
		alpharr[i]=i/size


	centre_arr= np.linspace(0,2*np.pi,size)
	theta=np.zeros(size)
	x_vals=np.zeros(size)
	
	fig = plt.figure(figsize=(11,11))
	ax = fig.add_subplot(111, aspect=1.0)
	
	
	for alpha in alpharr:
		theta[:]=alpha
		for centre in centre_arr:
			ew = np.array(np.linalg.eigvalsh(H(alpha,centre,u,v)))
			ax.plot(theta*np.ones_like(ew), (ew - 2*w*np.cos(2*np.pi*theta))/2, ",", color="blue", markersize=0.001)

	ax.set_xlim(0-1/40,1+1/40)
	#ax.set_xlim(0,1)
	ax.set_ylim(-0.5-1/40,0.5+1/40)
	ax.set_xticks([0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1])
	ax.set_xticklabels([r"0","",r"$0.25$","",r"$0.5$","",r"$0.75$","",r"1"])

	ax.set_yticks([-0.5,-0.375,-0.25,-0.125,0,0.125,0.25,0.375,0.5])
	ax.set_yticklabels([r"$-1$",r"",r"$-0.5$",r"",r"0",r"",r"0.5",r"",r"1"])


	ax.set_xlabel(r"$\theta=\frac{p}{q}$",fontsize=20)
	ax.set_ylabel(r"$\sigma(\hat{H}_{\theta})-2\pi\alpha_z\cos(2\pi\theta)$",fontsize=20)
	ax.set_title(r"Modified Hofstadter Butterfly",fontsize=22,pad=30,x=0.55)
	ax.tick_params(axis='both', which='major', width=1,length=6,
		                   direction="in", labelsize=16,zorder=5)
	ax.tick_params(axis='both', which='minor', width=1, length=2,
		                   direction="in", labelsize=4, zorder=5)
	ax.tick_params(which="both", top=True,labeltop=True,right=True,
		                   labelright=True, zorder=5)
	#ax.xaxis.set_major_locator(MultipleLocator(1))
	ax.xaxis.set_minor_locator(MultipleLocator(1/40))
	#ax.yaxis.set_major_locator(MultipleLocator(2))
	ax.yaxis.set_minor_locator(MultipleLocator(1/40))
	fig.subplots_adjust(left=0.15,right=0.85,top=0.85, bottom=0.15)
	plt.savefig("Hof16_16_16.pdf")

	#plt.show()
