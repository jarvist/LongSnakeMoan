#!/usr/bin/env python

# Version 1.5
# one dimensional chain

# Copyright under GNU General Public License 2010, 2012
# by Sinisa Coh and David Vanderbilt (see gpl-pytb.txt)

from pytb import * # import TB model class
import pylab as pl
import random
import datetime
from math import *

pl.ion()

#print random.gauss(0,0.1)

# Huckel type molecule...
# SiteE = Ionisation Potential
# J     = electron hopping integral
siteE=5.853060  #5.853060
J=0.6
trapE=0.176
siteSigma=0.05
JSigma=0.3 #0.5

# specify model
# Define Lattice Vectors
lat=[[1.0]]
# Define coordinates of orbitals
orb=[[0.0]]
#1D model please (k,r)
my_model=tbmodel(1,1,lat,orb)

# set on-site energies
my_model.set_sites([siteE])

my_model.add_hop(-J, 0, 0, [1])

my_model.display()

# solve model
# generate k-points
path=[[-0.5],[0.5]]
kpts=k_path(path,100)
#print "kpts, ",kpts

(evals,eigs)=my_model.solve_all(kpts,eig_vectors=True)

finite_model=my_model.cut_piece(100,0,glue_edgs=False) #glue_edgs makes system periodic

print finite_model._hoppings

for i in range(99):
    finite_model._site_energies[i]=finite_model._site_energies[i]+random.gauss(0,siteSigma)
    finite_model._hoppings[i][0]=finite_model._hoppings[i][0]*(cos(random.gauss(0,JSigma))**2)

finite_model.display()

#Solve the finite model!
(evals,eigs)=finite_model.solve_all(eig_vectors=True)

# start outputting useful info


print("Eigenvalues, Finite model"),evals
#Initialise our figures...
fig=pl.figure()

pl.subplot(411)
pl.title("Beta-PFO 1D pyTB - Disordered\nSigma=%s JSigma=%s"%(siteSigma,JSigma))

energies=[]

for i in range(75,85,1):
    finite_model._site_energies[i]=finite_model._site_energies[i]-trapE

F=0.01 #V/site, roughly equiv. to V/nm. 
#Thomas K says 1V/100nm is about right for an Organic Solar Cell at Short Circuit

#Apply electric field gradient across model
#for i in range(1,100):
#    finite_model._site_energies[i]=finite_model._site_energies[i]-(i*F)

colours='bgrcmykkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk'
  # Aaah, we fade to grey (fade to grey)
b=10
for i, colour in zip(range(25,45,b),colours):
    print "Iterate value i=",i," colour value",colour
    (evals,eigs)=finite_model.solve_all(eig_vectors=True)
    print evals[0]
    energies.append(evals[0]) #log of first eigenvalues

    #Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
    pl.subplot(411)
    for j in range(0,5):
        pl.fill_between(range(100),evals[j],evals[j]+eigs[j]*eigs[j], facecolor=colour)
    pl.ylabel("Occ %")
    pl.yticks(fontsize=9)

    #Plot Hamiltonian
    pl.subplot(412)
    pl.ylim(-1.1,-1.0)
    pl.ylabel("Hamiltonian\n(eV)")
    pl.yticks(fontsize=9)

    pl.plot(finite_model._site_energies,color=colour)

    #Plot Eigenvalues - not sure how useful this display is
    pl.subplot(413)
    pl.yticks(fontsize=9)

    pl.ylabel("Js (eV)")
    pl.plot(zip(*finite_model._hoppings)[0])

#    pl.ylabel("Eigenvalues\n(eV)")
#    pl.plot(evals,color=colour)

    #Plot DoS
    pl.subplot(414)
    pl.ylabel("DoS (#, eV)")
    pl.yticks(fontsize=9)
    pl.xlabel("Tight Binding Site (#) / DoS Histogram (eV)") #xLabel for shared axes

    pl.hist(evals,50,color=colour)

    pl.draw()

    for j in range(i,i+b): #filthy
        finite_model._site_energies[j]=finite_model._site_energies[j]-trapE #From NWCHEM / BNL calc on PFO
        finite_model._hoppings[j][0]=-J-0.3
    # See: /work/jmf02/NWCHEM/BETA-PHASE/tune_beta_alpha
    #finite_model._site_energies[10]=finite_model._site_energies[10]-0.01
 

#pl.subplot(413)
#pl.plot(energies)

#Top Subplot
#pl.plot(evals)

#For some reason this finally holds the window open
pl.show()

print "Saving figures...(one moment please)"
now=datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
fig.savefig("%s-BetaPhasePyTB.pdf"%now)
fig.savefig("%s-BetaPhasePyTB.png"%now)
