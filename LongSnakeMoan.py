#!/usr/bin/env python

# Long Snake Moan
#  - a Python Tighbinding simulation of torsionally disordered conjugated polymer in vacuum
# Jarvist Moore Frost 2012-2014

#Originally:
## Version 1.5 - originally pyTB example 'one dimensional chain'
## Copyright under GNU General Public License 2010, 2012
## by Sinisa Coh and David Vanderbilt (see gpl-pytb.txt)

from IPython import embed# we do this so we can drop to interactive python for debugging. 
 # Major python coolio
 # --> embed() <-- TA DA!

from pythtb import * # import TB model class --> now pythtb 1.6.2 version

import matplotlib.pyplot as pl

import random
import datetime
from math import *

#Pretty colours; via http://blog.olgabotvinnik.com/post/58941062205/prettyplotlib-painlessly-create-beautiful-matplotlib
try:
    import brewer2mpl #'pip install brewer2mpl'
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    colours = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
except ImportError: #If no brewer2mpl library
    #Otherwise, boring built in ones...
    colours='brgcmkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk' # Aaah, we fade to grey (fade to grey)
    print "Hey - no brew2mpl. Thus a simple palette."

# Matplotlib setup for publication quality figures...
pl.ion()

fig_size=[240.7/72, 240.7/72] #square aspect ratio, size of Latex 2 column width
params = {'backend': 'ps',
           'axes.labelsize': 9,
           'text.fontsize': 9,
           'legend.fontsize': 9,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'text.usetex': True,
           'figure.figsize': fig_size}
#pl.rcParams.update(params)
# ^- This is all currently borked; so don't use it.

print "'You ought to hear my Long Snake Moan' ~ PJ Harvey"
#print random.gauss(0,0.1)

# Model setup...
# Huckel type molecule...
# SiteE = Ionisation Potential
# J     = electron hopping integral
siteE=5.853060  #5.853060
J0=0.6 #Maximal hopping integral (inter-torsional theta=0)
trapE=0.176 #Derived from the DFT calculation on octamers
siteSigma=0.05
ThetaSigma=0.1 #0.5
ThetaMin=pi/4 #pi/4 #Pi/4
nsites=100

#Change data collection / plotting
DoS_averages=100
histbins=500

DoS=[]

# SPECIFY MODEL
# Define Lattice Vectors
lat=[[1.0]]
# Define coordinates of orbitals
orb=[[0.0]]
#1D model please (k,r)
my_model=tbmodel(1,1,lat,orb)
# set on-site energies
my_model.set_sites([siteE])
# constant J hop between sites (tridiagonal density matrix)
my_model.add_hop(J0, 0, 0, [1])
# Show me please
my_model.display()

# SOLVE MODEL
# generate k-points
path=[[-0.5],[0.5]]
kpts=k_path(path,nsites)
print "kpts, ",kpts

(evals,eigs)=my_model.solve_all(kpts,eig_vectors=True)

print "Infinite model:\n",evals

finite_model=my_model.cut_piece(nsites,0,glue_edgs=False) #glue_edgs makes system periodic

print finite_model._hoppings

for i in range(nsites-1):
    finite_model._site_energies[i]=finite_model._site_energies[i]+random.gauss(0,siteSigma)
    finite_model._hoppings[i][0]=finite_model._hoppings[i][0]*(cos(random.gauss(ThetaMin,ThetaSigma))**2) #random around 45 degrees

finite_model.display()

#Solve the finite model!
(evals,eigs)=finite_model.solve_all(eig_vectors=True)

# start outputting useful info

print("Eigenvalues, Finite model"),evals

#Initialise our figures...
fig=pl.figure()
pl.axes().set_aspect('equal')

ax1=pl.subplot(511) #5 subplots stacked on top of one another
pl.title("LongSnakeMoan - Disordered\nSigma=%s ThetaSigma=%s"%(siteSigma,ThetaSigma))

energies=[]

#F=0.01 #V/site, roughly equiv. to V/nm. 
#Thomas K says 1V/100nm is about right for an Organic Solar Cell at Short Circuit
#Apply electric field gradient across model
#for i in range(1,100):
#    finite_model._site_energies[i]=finite_model._site_energies[i]-(i*F)


def edit_hamiltonian(J0, ThetaMin, ThetaSigma, siteE, siteSigma, BetaSegment):
    for k in range(nsites-1):
        finite_model._site_energies[k]=siteE+random.gauss(0,siteSigma)
        finite_model._hoppings[k][0]=J0*(cos(random.gauss(ThetaMin,ThetaSigma))**2) #random around Theta Min
    if (BetaSegment>0):
        for j in range(nsites/2-BetaSegment/2,nsites/2+BetaSegment/2): #filthy
#        finite_model._site_energies[j]=finite_model._site_energies[j]-trapE #From NWCHEM / BNL calc on PFO
            finite_model._hoppings[j][0]=J0 #Perfect transfer integrals


betas=[0,8] # ,2,4,6,8] #Width of beta phase segments

# TODO: (re)write this as a function, so we can script for some data goodness.
# ( I'm not going to do this now as I'll make a mess )

for betasegments, colour in zip(betas,colours): #I'm ashamed of this nasty hack. JMF
    # this means you step over the iterator value, and have a different colour plot for each loop
    print "Iterate value betasegments=",betasegments," colour value",colour

    DoS.extend([4.50]) # nasty hack to get Beta + alpha phase DoS to plot with same x-axes
    DoS.extend([7.5])   #  ^- which doesn't work on old MatPlotLib...
#    print DoS  # Nb: very verbose - all eigenvalues

    for a in range(DoS_averages):
        print "DoS average ",a
        edit_hamiltonian(J0,ThetaMin,ThetaSigma,siteE,siteSigma,betasegments)
        evals=finite_model.solve_all(eig_vectors=False) # only DoS here please
        #NB: internally pythtb just uses 'eigvalsh' which is the standard Hermition solver
        # there is also scipy.sparse.linalg.eigsh which may be sig. faster
        DoS.extend(list(evals))

    #Just one sample of eigenvectors + eigenvalues for plotting wave functions
    (evals,eigs)=finite_model.solve_all(eig_vectors=True)
   
    energies.append(evals[0]) #log of first eigenvalues; for trap-depth calculation

    #Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
    pl.subplot(511)
    for j in range(0,5): #Number of eigenvalues plotted (electron wavefns)
        pl.fill_between(range(nsites),evals[j],evals[j]+eigs[j]*eigs[j], facecolor=colour)
    pl.ylabel("Occ %")
    #pl.ylim((3.8,5.0))
    pl.yticks(fontsize=9)
    pl.xticks(visible=False)

    #Plot Hamiltonian
    pl.subplot(512, sharex=ax1)
    #pl.ylim((5.5,6.0))
    pl.ylabel("S (eV)")
    pl.yticks(fontsize=9)
    pl.xticks(visible=False)

    pl.plot(finite_model._site_energies,color=colour)

    #Plot Eigenvalues - not sure how useful this display is
    pl.subplot(513,sharex=ax1)
   # pl.ylim((0.0,1.0))
    pl.yticks(fontsize=9)
    pl.xticks(visible=False)

    pl.ylabel("J (eV)")
    pl.plot(zip(*finite_model._hoppings)[0])

#    pl.ylabel("Eigenvalues\n(eV)")
#    pl.plot(evals,color=colour)

    #Plot DoS
    if (betasegments>0): # this is a beta phase region .'. 5th plot
        pl.subplot(515,sharex=ax2) #Nb: order matters to share axis...
        pl.xlabel("Eigenvalue (eV)")
    else:
        ax2=pl.subplot(514) # non-beta-phase on the 4th plot
        pl.xticks(visible=False)

    pl.ylabel("DoS")
    pl.yticks(fontsize=9)

    #embed() #iPython interactive session - squash some bugs!

    print "Colour: ", colour
    # Histogram of TB Eigenvalues (i.e. DoS)
    #pl.hist(evals,50,color=colour)  # Colour broken with OS X matplotlib - complains about number of colours vs. number of data in set
    pl.hist(DoS,histbins,histtype='stepfilled',color=colour)
    # Histogram of Hopping Integrals (funky zip command to rearrange J magnitude into 1D vector)
#    pl.hist(zip(*finite_model._hoppings)[0],50,color=colour)
    
    DoS=[]

# OK; wrap up time! Save those plots

pl.tight_layout(pad=0.3) #, w_pad=0.5, h_pad=1.0) # Magic incantation for non-terrible plots
pl.show()
#For some reason this finally holds the window open
## - I don't understand this at all, might be something to do with a dated version of matplotlib on Ubuntu 10.04

print "Lowest Eigenvalues:\n", energies
print "Trap Depth (Min - max):",energies[0]-energies[-1]

print "Saving figures...(one moment please)"
now=datetime.datetime.now().strftime("%Y-%m-%d-%H:%M") #String of standardised year-leading time
pl.annotate("JMF %s"%now,xy=(0.75,0.02),xycoords='figure fraction') #Date Stamp in corner
pl.show() #Show the plot

fig.savefig("%s-LongSnakeMoan.pdf"%now) #Save figures as both PDF and easy viewing PNG (perfect for talks)
fig.savefig("%s-LongSnakeMoan.png"%now)
