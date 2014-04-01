# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

#using Roots, Polynomial

using Calculus

println("Sturm und Drang: DoS by Sturm sequences")

N=100

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
B=1/(300*kB) #300K * k_B in eV

U(theta)=( 1.0+cos(4*theta) + P*abs(theta) ) #defined as a function for further maths

Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) #Attempting to calculate partition function directly; this is correct

println("Partition function Z=",Z)

# Following checks the partition function code, outputting p(robability) as a fn(theta) for varying P
outfile=open("data.dat","w+")
for P=0:1.0:5
Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) # recalculate Z now that P is changing
    for t = -pi:(pi/1800.0):pi
#    println("Partition function Z=",Z)
        p=exp(-U(t)*B)/Z
#        println(t," ",p)
        @printf(outfile,"%f %f %f\n",t,U(t),p)
    end
end

close(outfile)

#end

# TODO: Some clever physics here.

function randexp(N) # random exponential with the ln(X) 0<X<1 method
    return (log(rand(N)))
end

#p=Poly([])
#fzero(p)

function randH()
# Random Trace / diagonal elements
    D=5.0 + 0.0*randn(N)
# Random Off-diag elements
#E=0.1+0.05*randn(N-1)

#Generate thetas...
#    thetas=randexp(N-1)
    thetas=Float64[]
    for i=1:N-1
        theta=0.0
        while true 
            theta=2*pi*rand()
            p=exp(-U(theta)*B)/Z
            p>rand() && break
        end
        push!(thetas,theta)
    end

#Transfer integral from
    J0=0.500 #Max Transfer integral
    E=J0*cos(thetas).^2
# Squared (element wise) for the Sturm algorithm...
    E= E.^2
    return (D,E)
end

D,E=randH()

#println(D)
#println(E)

for P=0:1.0:5
    Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) # recalculate Z now that P is changing
    D,E=randH()
    @printf("Calculated with P= %f Z= %f\n",P,Z)
#println("STURM sequence method...")
    sigma=4.0
    for sigma=4.0:0.2:5.5
        @printf("%f %f\n", sigma , sturm(D,E,sigma))
    end
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
