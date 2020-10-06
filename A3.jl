# setting working directory
# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")
# mkpath("Assignment3")
#cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment3")
# mkpath("graphs")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
# Pkg.add("ForwardDiff")
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
Pkg.add("LaTeXStrings")
using LaTeXStrings
using Latexify

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

function poly_inter_err(a,b,n,f::Function)
    a , b  = min(a,b), max(a,b)
    h = (b-a)/n
    x = [a+(i-1)*h for i in 1:(n+1)]
    Π = h^(n+1)*factorial(n)/4
    ξ = abs(maximum(f.(x)))
    err = Π*ξ/factorial(n+1)
    return err
end


err=[poly_inter_err.(0.05,2,i,funs) for i in [4,6,11]]
Err=[err[1] err[2] err[3]]

copy_to_clipboard(true)
mdtable(Err, side=nom, head=["n=4","n=6","n=11"],fmt = FancyNumberFormatter(4))

## Cubic Splines - Natural

#The function CubicNatural takes the x and y-coordinates of the data as input, and
# computes the natural cubic spline interpolating the data, by solving the resulting matrix
# equation. The code is based on Algorithm 3.4 of Burden, Faires, Burden [4]. The output
# is the coefficients of the m 􀀀 1 cubic polynomials, ai; bi; ci; di; i = 1; :::;m 􀀀 1 where m is
# the number of data points. These coefficients are stored in the arrays a; b; c; d, which are
# declared global, so that we can access these arrays later to evaluate the spline for a given
# value w.

function CubicNatural(x::Array,y::Array)
    m=length(x) # m is the number of data points
    n=m-1
    a=Array{Float64}(undef,m)
    b=Array{Float64}(undef,n)
    c=Array{Float64}(undef,m)
    d=Array{Float64}(undef,n)
    a = y
    # for i in 1:m
    #     a[i]=y[i]
    # end
    h = [x[i+1]-x[i] for i in 1:n]
    # h=Array{Float64}(undef,n)
    # for i in 1:n
    #     h[i]=x[i+1]-x[i]
    # end
    u=Array{Float64}(undef,n)
    u[1]=0
    for i in 2:n
        u[i]=3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1]
    end
    s=Array{Float64}(undef,m)
    z=Array{Float64}(undef,m)
    t=Array{Float64}(undef,n)
    s[1]=1
    z[1]=0
    t[1]=0
    for i in 2:n
        s[i]=2*(x[i+1]-x[i-1])-h[i-1]*t[i-1]
        t[i]=h[i]/s[i]
        z[i]=(u[i]-h[i-1]*z[i-1])/s[i]
    end
    s[m]=1
    z[m]=0
    c[m]=0
    for i in reverse(1:n)
        c[i]=z[i]-t[i]*c[i+1]
        b[i]=(a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2*c[i])/3
        d[i]=(c[i+1]-c[i])/(3*h[i])
    end
    return a, b, c, d
end

# CubicNaturalEval. The inputs are the value at which
# the spline is evaluated, w, and the x-coordinates of the data. The function first finds the
# interval [xi; xi+1]; i = 1; :::;m-1; w belongs to, and then evaluates the spline at w using the
# corresponding cubic polynomial.

function CubicNaturalEval(x::Array,y::Array,w)
    m=length(x)
    if w<x[1]||w>x[m]
        return print("error: spline evaluated outside its domain")
    end
    n=m-1
    p=1
    for i in 1:n
        if w<=x[i+1]
            break
        else
            p=p+1
        end
    end
    # p is the number of the subinterval w falls into, i.e., p=i means
    # w falls into the ith subinterval $(x_i,x_{i+1}), and therefore
    # the value of the spline at w is
    # a_i+b_i*(w-x_i)+c_i*(w-x_i)^2+d_i*(w-x_i)^3.
    a, b, c, d = CubicNatural(x,y)
    return a[p]+b[p]*(w-x[p])+c[p]*(w-x[p])^2+d[p]*(w-x[p])^3
end

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
        interp=map(z->CubicNaturalEval(xi,yi,z),xax) # Interpolating poly for the data
        # Plot
        gr()
        plot(title="Interpolation $name n=$n - Cubic Spline Natural")
        plot!(xax,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(xax,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("graphs/CSN_$name _n_$n")
    end
end
