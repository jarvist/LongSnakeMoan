# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

#using Roots, Polynomial

println("Sturm und Drang: DoS by Sturm sequences")

N=100

# Calculates number of eigenvalues less than 'sigma' in tridiagonal matrix 
# described by: diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1)
# Nb: Off digonal elements E are supplied squared for computational speed
function sturm(D,E,sigma)
    n=length(D)
    t=zeros(n)
    countnegatives=0
    
    t[1]=D[1]-sigma #first term of sequence, to avoid needing E[0], t[0]
    if t[1]<0.0
        countnegatives=countnegatives+1
    end

    for i=2:n
        t[i]=D[i]-sigma-E[i-1]/t[i-1]   # Sturm sequence
        if t[i]<0.0                     # if t<0, we've found another eigenvalue
            countnegatives=countnegatives+1 
        end
    end

    return countnegatives
end

# SQUEEZED PFO FUNCTIONS
# Define potential function
P=0.05
B=0.025 #300K * k_B in eV

U(theta)=cos(4*theta)+P*theta #defined as a function for further maths
#Z=sum(exp(-U/B),[0:2*pi]) #Attempting to calculate partition function directly; this is not correct

# Random Trace / diagonal elements
D=5.0+0.1*randn(N)
# Random Off-diag elements
#E=0.1+0.05*randn(N-1)

function randexp(N) # random exponential with the ln(X) 0<X<1 method
    return (log(rand(N)))
end

#p=Poly([])
#fzero(p)

#Generate thetas...
thetas=randexp(N-1)
#Transfer integral from
J0=0.500 #Max Transfer integral
E=J0*cos(thetas).^2
# Squared (element wise) for the Sturm algorithm...
E= E.^2

println("STURM sequence method...")
sigma=4.0
while sigma<=5.5
    @printf("%f %f\n", sigma , sturm(D,E,sigma))
    sigma=sigma+0.02
end

#println("Elements...(offdiag^2, diag))");
#println([E;D])

#println("Full square matrix H");
H=diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1) #build full matrix from diagonal elements; for comparison
#println(H)

println("Eigenvalues")
println(eigvals(H))
println("Min Eigenvalue")
println(eigmin(H))

println("I time a single histogram (Es<sigma) point:")
@time sturm(D,E,5.0)
println("I time eigvals(H)")
@time eigvals(H)
