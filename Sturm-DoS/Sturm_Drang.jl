# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

println("Sturm und Drang: DoS by Sturm sequences")

N=10000

function sturm(D,E,sigma)
    n=length(D)
    t=zeros(n)
    countnegatives=0
    
    t[1]=D[1]-sigma #first term of sequence, to avoid needing E[0], t[0]
    if t[1]<0.0
        countnegatives=countnegatives+1
    end

    for i=2:n
        t[i]=D[i]-sigma-E[i-1]/t[i-1]
        if t[i]<0.0
            countnegatives=countnegatives+1
        end
    end

    return countnegatives
end

# Random Trace / diagonal elements
D=5.0+0.1*randn(N)
# Random Off-diag elements
E=0.1+0.05*randn(N-1)
# Squared (element wise) for the Sturm algorithm...
E= E.^2

println("STURM sequence method...")
sigma=4.5
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
