#!/usr/bin/env python

# Version 1.5
# one dimensional chain

# Copyright under GNU General Public License 2010, 2012
# by Sinisa Coh and David Vanderbilt (see gpl-pytb.txt)

from pytb import * # import TB model class
import pylab as pl

# Huckel type molecule...
# SiteE = Ionisation Potential
# J     = electron hopping integral
siteE=5.853060  #5.853060
J=0.6
trapE=0.176

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

(evals,eigs)=finite_model.solve_all(eig_vectors=True)
print ("Eigenvalues, non beta phase"),evals

#OK, a bith filthy here, there's no access method for the site_energies themselves...
#self._site_energies=nu.array(en_list)
#finite_model._site_energies[10]=-1.1

finite_model.display()

#Solve the finite model!
(evals,eigs)=finite_model.solve_all(eig_vectors=True)

# start outputting useful info

print("Eigenvalues, Finite model"),evals

fig=pl.figure()
pl.subplot(311)

(evals,eigs)=finite_model.solve_all(eig_vectors=True)
#pl.plot(eigs[0]*eigs[0])

pl.title("Beta Phase - 1D Tight Binding Model")

energies=[]

#OK; now we have a 'Alpha phase' bit of Huckel tight binding
(evals,eigs)=finite_model.solve_all(eig_vectors=True)
print evals[0]
energies.append(evals[0])
pl.subplot(311)
pl.plot(eigs[0]*eigs[0])
pl.subplot(312)
pl.plot(finite_model._site_energies)

for i in range(75,85,1):
    finite_model._site_energies[i]=finite_model._site_energies[i]-trapE

for i in range(25,32,1):
    finite_model._site_energies[i]=finite_model._site_energies[i]-trapE #From NWCHEM / BNL calc on PFO
    # See: /work/jmf02/NWCHEM/BETA-PHASE/tune_beta_alpha
    #finite_model._site_energies[10]=finite_model._site_energies[10]-0.01
    (evals,eigs)=finite_model.solve_all(eig_vectors=True)
    print evals[0]
    energies.append(evals[0]) #log of first eigenvalues
    pl.subplot(311)
    pl.plot(eigs[0]*eigs[0])
    pl.plot(eigs[1]*eigs[1])
    pl.subplot(312)
    pl.ylim(-1.1,-1.0)
    pl.plot(finite_model._site_energies)
# plot band structure

pl.subplot(313)
pl.plot(energies)

#Top Subplot
#pl.hist(evals,50)
#pl.plot(evals)


#Bottom Subplot
#pl.subplot(212)
#pl.ylim(-1.,1.)
#pl.plot(eigs[0])
#pl.plot(eigs[0]*eigs[0])
#pl.plot(finite_model._site_energies)
#print("Eigenvectors, first band, orbital"),eigs[0,:]


pl.show()
#pl.savefig("beta.pdf")


#fig=pl.figure()
#ax=fig.add_subplot(111)
#x=evals
#line,=ax.plot(x,eigs[0]*eigs[0])
#def animate(i):
#    finite_model._site_energies[10]=finite_model._site_energies[10]-0.01
#    (evals,eigs)=finite_model.solve_all(eig_vectors=True)
#    print finite_model._site_energies[10],evals[0]
#    line.set_ydata(eigs[0]*eigs[0])
#    return line,
#def init():
#    line.set_ydata(np.ma.array(x,mask=True))
#    return line,
#ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init, interval=25, blit=True)
#pl.show()

#print("Eigenvectors, Finite model"),eigs


