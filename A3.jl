# setting working directory
# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")
# mkpath("Assignment3")
# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment3")
# mkpath("graphs")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl

# define Functions

u1(x) = log(x)
u2(x) = (x)^(1/2)
u3(x,s) = x^(1-s)/(1-s)
u4(x) = u3(x,2)
u5(x) = u3(x,5)
u6(x) = u3(x,10)

# define derivatives
u1p(x) = 1/x
u2p(x) = (1/2)*x^(-1/2)
u3p(x,s) = x^(-s)

xax = range(0.05,2;length=1000)

u1.(xax)

# Polynomial Interpolation: Newton Basis
# Code from Giray Otken, First semester in numerical analysis with Julia

    function diff(x::Array,y::Array)
        m = length(x) #here m is the number of data points. #the degree of the polynomial is m-1 a=Array{Float64}(undef,m)
        a = [y[i] for i in 1:m]
        # for i in 1:m
        #     a[i]=y[i]
        # end
        #a = [(y[i]-y[i-1])/(x[i]-x[i-j+1]) for j in 2:m, i in reverse(collect(j:m))]
        for j in 2:m
            for i in reverse(collect(j:m))
                a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
            end
        end
        return a
    end

    function newton(x::Array,y::Array,z)
        m=length(x) #here m is the number of data points, not the degree # of the polynomial
        a=diff(x,y)
        sum=a[1]
        pr=1.0
        for j in 1:(m-1)
            pr=pr*(z-x[j])
            sum=sum+a[j+1]*pr
        end
        return sum
    end

    # Runge's Example with Newton Polynomial for different n


nom = ["Log","Square","CES σ=2", "CES σ=5", "CES σ=10"]
funs = [u1,u2,u4,u5,u6]
for i in 1:5
    name = nom[i]
    for n = (4,6,11,21)
        fn = funs[i].(xax)
        # Grid of nodes for interpolation
        xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
        yi = funs[i].(xi) # the corresponding y-coordinates
        # Interpolation
        interp=map(z->newton(xi,yi,z),xax) # Interpolating poly for the data
        # Plot
        gr()
        plot(title="Interpolation $name n=$n - Newton Polynomial")
        plot!(xax,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(xax,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("graphs/Newton $name _n_$n")

    end
end
