from scipy.linalg import eigvalsh 
from scipy.linalg import eigh
from numpy import *
import matplotlib.numerix as nx
import pylab as p

Har = 27.21 #1 Hartree in eV

LAMBDA = 0.5 #lattice reorganisation energy in units eV
K_BOLTZMANN = 8.6173857e-5 #in units of eV
T = 300.0 #in Kelvin

THERMAL=(1.0/(4*LAMBDA*K_BOLTZMANN*T)) #frac{1}{4 \lambda k_{B} T}

# the function takes as argument:
# n-> number of sites
# d-> the distance of the PCBM molecule (in Bohr Radii)
# a-> the size of each monomer (in Bohr Radii)
# eps-> the relative dielectric constant
# t -> the transfer integral along the chain (in eV)

def Hamiltonian(n,d,a,eps,t):
	H = zeros( (n,n) )
	for i in range(0,n):
		# the distance along the chain 
		xd = a* abs( (n-1)/2. - i) 
		# print xd
		# the distance of the monomer to PCBM
		r = sqrt(xd*xd+d*d)
		#the coulomb potential at that point (in eV)
		E = - Har /(eps *r)
		#set diagonal element:	
		H[i,i] = E

		#set off diagonal element
		if ( i>1 and i<n-2):
			H[i,i+1] = -t
			H[i+1,i] = -t
	return H

def BetaHamiltonian(n,b,E,j,J):
    H = zeros( (n,n) )
    for i in range(0,n):
        H[i,i] = E
        #set off diagonal element
        if ( i>0 and i<n-2):
            if (i<(n/2.)-(b/2.) or i >= (n/2.)+(b/2.)):
                H[i,i+1] = j
                H[i+1,i] = j
            else:
                H[i,i+1] = j #J
                H[i+1,i] = j #J
                #messing up the site energy / trace for BETA. Should really be both this site + next one
                H[i][i] = H[i][i] - J 
    H[0][0]=0 #Zero potential at end of chain...?
    H[n-1][n-1]=0
    print 
    print H
    return H


def MasterEquation(n,H):
    M=zeros( (n,n) )
    for i in range(0,n):
        for j in range(0,n):
            dG=H[i,i]-H[j,j] #Site energy difference, from Hamiltonian
            dE=-LAMBDA+dG-(j-i)*0.01 #0.01 is ElectricField
            M[i,j]=H[i,j]*H[i,j]*exp(dE*dE*THERMAL)   #MARCUS HOPPING RATE
    for i in range (0,n):
        M[i,i]=0.0
        M[i,i]=-sum(M[0:n,i]) #slice, to get [i,i] trace = negative sum of rates
    return M

def GetLowestEV(H):
    vals = eigvalsh(H)
    print  
    print vals 
    return vals[0]

def GetLowestEVrinf(n, t):
	H = zeros( (n,n) )
	for i in range(0,n):
		#set off diagonal element
		if ( i!=n-1):
			H[i,i+1] = -t
			H[i+1,i] = -t
	vals = eigvalsh(H)
	return vals[0]
	

# following lines are just quick tests
#H=Hamiltonian(4, 20,10,1,0.1)
#val = eigvalsh(H)
#print H
#print val

#val=GetLowestEV(Hamiltonian(10,20,10,1,0.1))
#print val

distance = nx.arange(6.,20.,1.) # distance from 6 to 20 angstroms 
dAU = distance / 0.527       # convert to bohrs (I think?)
aAU =      7.   / 0.527       # lenght of a monomer in bohr radii
eps =      4.
n   =      21  
t   =      1.


n=4 #needs to be even?
E=0.1
j=-0.06
J=-0.6

#energies  = zeros(n)
#for i in range(2,n):
#    energies[i] =  GetLowestEV(BetaHamiltonian(i,2,E,j,J))

print "Betaphase"
print "Eigenvectors"
BetaH=BetaHamiltonian(n,4,E,j,J)
#w,v=eigh(BetaH)
#print "Eigenvalues", w
#print "Eigenvectors", v
#print "first Eigenvector..." 
#print v[0]
#print "Eigenvector squared..."
#print v[0]*v[0]

M=MasterEquation(n,BetaH)

print "Hamiltonian:\n",BetaH
print "Master Equation Rates:\n",M

p.plot (
        range(0,n),v[0]*v[0])
#        range(0,n),v[1]*v[1],
#        range(0,n),v[1]*v[1]+v[0]*v[0])
        #        range(0,n),v[1]*v[1],
#        range(0,n),v[1]*v[1]+v[0]*v[0])
#,
#        range(0,n),v[2]*v[2], 
#        range(0,n),energies)
#p.legend(('Hamiltonian','Beta - Eigenvector_0^2','Beta - Eigenvector_1^2','Beta - Eigenvector_0^2 + Eigenvector_1^2'))
p.show()



