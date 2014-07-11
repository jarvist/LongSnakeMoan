# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

#using Roots, Polynomial

using Calculus

println("Sturm und Drang: DoS by Sturm sequences")

N=1000000

# Calculates number of eigenvalues less than 'sigma' in tridiagonal matrix 
# described by: diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1)
# Nb: Off digonal elements E are supplied squared for computational speed
function sturm(D,E,sigma)
    t=0.0
    countnegatives=0
    
    t=D[1]-sigma #first term of sequence, to avoid needing E[0], t[0]
    if t<0.0
        countnegatives=countnegatives+1
    end

    for i=2:length(D)
        t=D[i]-sigma-E[i-1]/t   # Sturm sequence, overwriting temporary files...
        if t<0.0                     # if t<0, we've found another eigenvalue
            countnegatives=countnegatives+1 
        end
    end

    return countnegatives
end

#units eV
kB=8.6173324E-5

# SQUEEZED PFO FUNCTIONS
# Define potential function
P=0.05

#U(theta)=( 1.0+cos(4*theta) + P*abs(theta) ) #defined as a function for further maths
#U(theta)=( 0.05* (cos(4*theta) + cos(2*theta)) + P*sin(theta)^2 ) # More physical form - increasing cost for theta=90 degrees as fn of pressure
#U(theta)=( 0.05 * sin(theta)^2 ) #P3HT like

#Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) #Attempting to calculate partition function directly; this is correct

#println("Partition function Z=",Z)

#end

# TODO: Some clever physics here.

function randexp(N) # random exponential with the ln(X) 0<X<1 method
    return (log(rand(N)))
end

#p=Poly([])
#fzero(p)

function randH(disorder,B,Z,U)
# Random Trace / diagonal elements
    D=5.0 + disorder*randn(N)
# Random Off-diag elements
#E=0.1+0.05*randn(N-1)

#Generate thetas...
#    thetas=randexp(N-1)
    thetas=Float64[]
    for i=1:N-1 #number of thetas we need for off-diagonal element
        theta=0.0
        while true # this is a do-while loop, Julia styleee
            theta=2*pi*rand()     #random theta
            p=exp(-U(theta)*B)/Z  #probability by stat mech
            p>rand() && break     #rejection sampling of distribution
        end
        push!(thetas,theta)       #tack theta onto end of list of thetas
    end

#Transfer integral from
    J0=0.800 #Max Transfer integral
    E=(J0*cos(thetas)).^2 
# Squared (element wise) for the Sturm algorithm...
    E= E.^2
    return (D,E)
end

outfile=open("DoS.dat","w+")
onsetfile=open("onset.dat","w+")
potfile=open("potential.dat","w+")

P=0.0
disorder=0.0

# OK; Raos FF paper (Moreno et al. J.Phys.Chem.B 2010)
# 'full' potential energy Fig 4.a. puts a barrier at 90 degress of ~3.0 kCal / mol = 126 meV
E0=0.126

@sync @parallel for T=300.0:0.1:300.0 #T=100.0:100:400 #:0.1:1
    B=1/(T*kB) #300K * k_B in eV
    U(theta)=( E0 * sin(theta)^2 ) #P3HT like
    Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) # recalculate Z now that P is changing
    println("Partition function for Z(E0=",E0,")=",Z)

# Following checks the partition function code, outputting p(robability) as a fn(theta) for varying P
    for t = -pi:(pi/180.0):pi
#    println("Partition function Z=",Z)
        p=exp(-U(t)*B)/Z
#        println(t," ",p)
        @printf(potfile,"%f %f %f\n",t,U(t),p)
    end
e
    
    D,E=randH(disorder,B,Z,U)
    @printf("Calculated with P= %f Z= %f\n",P,Z)
#println("STURM sequence method...")
    sigma=4.0
    pDoS=0
    pold=0
    onset=false
    for sigma=6:0.001:6.4 #sigma=3:0.010:7 #sigma=6:0.001:6.4 #6 to 6.4 covers right lobe
        pDoS=sturm(D,E,sigma)
        @printf("%f %f %f %f\n", T, sigma , pDoS, pDoS-pold)
        @printf(outfile,"%f %f %f %f\n",T,sigma,pDoS, pDoS-pold)
        if (pDoS-pold>10.0 && onset==false)
            @printf(onsetfile,"%f %f %f\n",P,sigma,pDoS-pold)
            onset=true
        end
        pold=pDoS
    end
end

close(outfile)
close(onsetfile)

end

#println("Elements...(offdiag^2, diag))");
#println([E;D])

#println("Full square matrix H");
H=diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1) #build full matrix from diagonal elements; for comparison
println(H)

println("Eigenvalues")
println(eigvals(H))
println("Min Eigenvalue")
println(eigmin(H))

println("I time a single histogram (Es<sigma) point:")
@time sturm(D,E,5.0)
println("I time eigvals(H)")
@time eigvals(H)
