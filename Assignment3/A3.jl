
# setting working directory
# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")
# mkpath("Assignment3")

# mkpath("graphs")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
Pkg.add("Dierckx")
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
Pkg.add("Interpolations")
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
Pkg.add("ForwardDiff")
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
Pkg.add("LaTeXStrings")
using LaTeXStrings
Pkg.add("Latexify")
using Latexify
Pkg.add("PyPlot")


cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment3")

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
u4p(x) = x^(-2)
u5p(x) = x^(-5)
u6p(x) = x^(-10)

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
mdtable(Err, side=nom, head=["n=4","n=6","n=11"],fmt = FancyNumberFormatter(4))|> print
mdtable(Err, side=nom, head=["n=4","n=6","n=11"],fmt = FancyNumberFormatter(4))|> display

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


## Schumaker shape preserving Spline

##find_xi funtion that takes in x and y coordinates
##and slopes si at each xi and returns the ξ parameter
function find_xi(t::Array,z::Array,s::Array)
    #Based on Algorithm Tudd pag 233
    n = length(t)
    δ = [(z[i+1]-z[i])/(t[i+1]-t[i]) for i in 1:(n-1)]
    ξ = zeros(n-1)
    #Check if lemma 6.11.1 applies
    for i in 1:(n-1)
        if (s[i] + s[i+1])/2 == δ[i]
            ξ[i] = t[i+1]
        elseif (s[i]-δ[i])*(s[i+1]-δ[i]) ≥ 0
            ξ[i] = (t[i]+t[i+1])/2
        elseif abs(s[i+1]-δ[i]) < abs(s[i]-δ[i])
            ξub = t[i]+ (2*(t[i+1]-t[i])*(s[i+1]-δ[i]))/(s[i+1]-s[i])
            ξ[i] = (t[i]+ξub)/2
        else # abs(s[i+1]-δ[i]) ≥ abs(s[i]-δ[i])
            ξlb = t[i+1]+ (2*(t[i+1]-t[i])*(s[i]-δ[i]))/(s[i+1]-s[i])
            ξ[i] = (t[i+1]+ξlb)/2
        end
    end
    return ξ, δ
end

#Schumaker takes in x, y coordinates, slopes

function Schumaker(t::Array,z::Array,s::Array,w)
    m = length(t)
    if w < t[1] || w > t[m]
        return print("error: spline evaluated outside its domain")
    end
    n = m-1
    ξ, δ = find_xi(t,z,s)
    i = 1
    for j in 1:n
        if w ≤ t[j+1]
            break
        else
            i=i+1
        end
    end
    α = ξ[i]-t[i]
    β = t[i+1]-ξ[i]
    sb = (2*(z[i+1]-z[i])-(α*s[i]+β*s[i+1]))/(t[i+1]-t[i])
    A1, B1, C1 = z[i], s[i], (sb-s[i])/(2*α)
    A2, B2, C2 = A1 + α*B1 + C1*(α^2), sb, (s[i+1]-sb)/(2*β)
    if w ≤ ξ[i]
        w_t = A1 + B1*(w - t[i]) + C1*(w-t[i])^2
    else
        w_t = A2 + B2*(w - ξ[i]) + C2*(w-ξ[i])^2
    end
    return w_t, ξ, δ
end


# logfun = u1.(xax)
# # Grid of nodes for interpolation
# xi = collect(range(0.05,2;length=6)) ; # Collect makes it an array instead of a collection
# yi = u1.(xi) # the corresponding y-coordinates
# si = u1p.(xi) #the corresponding slopes
# xip = find_xi(xi,yi,si)[1]
#
# plot(xax,logfun)
# vline!(xi, color=:blue)
# vline!(xip, color=:red)
#
# # Interpolation
# interp=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
# # Plot
# gr()
# plot(xax,logfun,linewidth=3,legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
# plot!(xax,interp,linewidth=3,label="Interpolation")
# plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
#
# xi = [0.05,0.10]
# yi = u1.(xi)
# si = u1p.(xi)
#
# find_xi(xi,yi,si)
# Schumaker(xi,yi,si,0.075)

nom = ["Log","Square Root","CES σ=2", "CES σ=5", "CES σ=10"]
funs = [u1,u2,u4,u5,u6]
slopes = [u1p,u2p,u4p,u5p,u6p]
for i in 1:5
    name = nom[i]
    for n = (4,6,11,21)
        fn = funs[i].(xax)
        # Grid of nodes for interpolation
        global xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
        global yi = funs[i].(xi) # the corresponding y-coordinates
        global si = slopes[i].(xi) #the corresponding slopes
        # Interpolation
        interp=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
        # Plot
        gr()
        plot(title="Interpolation $name n=$n - Schumaker Shape Preserving Spline")
        plot!(xax,fn,linewidth=3,label = "Function: $name",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(xax,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("graphs/Schu_$name _n_$n")
    end
end

## graphs
##
fnom = ["Log","Square Root","CES "*latexinline("σ=2"), "CES "*latexinline("σ=5"), "CES "*latexinline("σ=10")]
funs = [u1,u2,u4,u5,u6]
slopes = [u1p,u2p,u4p,u5p,u6p]
method = [newton, CubicNaturalEval, Schumaker]
gridpts = [4,6,11,21]
nom = ["Log","SQR","CES_2", "CES_5", "CES_10"]
titlenom = ["Log","Square Root","CES σ=2", "CES σ=5", "CES σ=10"]
#titlenom = ["Log_n=6","Square_Root n=6","CES σ=2 n=6", "CES σ=5 n=6", "CES σ=10 n=6"]
for i in 1:length(funs)
    plotname = nom[i]
    title = titlenom[i]
    name = fnom[i]
    fn = funs[i].(xax)
    IE = Array{Float64}(undef, 3, 4)
    local interp = Array{Float64}(undef,1000,3)
    for n = gridpts
        for j in 1:length(method)
            met = method[j]
            # Grid of nodes for interpolation
            global xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
            global yi = funs[i].(xi) # the corresponding y-coordinates
            # Interpolation
            if j ≤ 2
                global interp[:,j] = map(z->method[j](xi,yi,z)[1],xax)
            else
                global si = slopes[i].(xi) #the corresponding slopes
                global interp[:,j]=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
            end
            #IE[j,n] = interp_err_SN(fn,interp)
            #IE[j,n] = sqrt(sum((fn .- interp).^2))
        end
        global fmet = ["$name","Newton Basis Polynomials","Natural Cubic Splines","Shape-preserving Schumaker Splines"]
        # Plot
        gr()
        plot(title="Interpolation $title Grid size n=$n")
        plot!(xax,fn,linewidth=3,label = "Function: $title",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(xax,interp[:,1],linewidth=3,label=fmet[2], linestyle=:dash)
        plot!(xax,interp[:,2],linewidth=3,label=fmet[3],linestyle=:dot)
        plot!(xax,interp[:,3],linewidth=3,label=fmet[4], linestyle=:dashdot)
        plot!(xi,yi,linetype=:scatter,marker=(:circle,6),markercolor=:blue,label = "Data")
        savefig("graphs/$plotname $n.png")
    end
end


## Errors
##

fnom = ["Log","Square Root","CES "*latexinline("σ=2"), "CES "*latexinline("σ=5"), "CES "*latexinline("σ=10")]
funs = [u1,u2,u4,u5,u6]
slopes = [u1p,u2p,u4p,u5p,u6p]
method = [newton, CubicNaturalEval, Schumaker]
gridpts = [4,6,11,21]
nom = ["Log","SQR","CES_2", "CES_5", "CES_10"]
titlenom = ["Log","Square Root","CES σ=2", "CES σ=5", "CES σ=10"]
for i in 1:length(funs)
    plotname = nom[i]
    title = titlenom[i]
    name = fnom[i]
    fn = funs[i].(xax)
    IE = Array{Float64}(undef, 3, 4)
    for j in 1:length(method)
        for n in 1:length(gridpts)
            # Grid of nodes for interpolation
            global xi = collect(range(0.05,2;length=gridpts[n])) ; # Collect makes it an array instead of a collection
            global yi = funs[i].(xi) # the corresponding y-coordinates
            # Interpolation
            if j ≤ 2
                global interp = map(z->method[j](xi,yi,z)[1],xax)
            else
                global si = slopes[i].(xi) #the corresponding slopes
                global interp=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
            end
            IE[j,n] = sqrt(sum((fn .- interp).^2))
        end
    end
    global fmet = ["$name","Newton Basis Polynomials","Natural Cubic Splines","Shape-preserving Schumaker Splines"]
    mdtable(IE, side=fmet, head=["n=4","n=6","n=11","n=21"],fmt = FancyNumberFormatter(4))|> print
    # Plot
    gr()
    plot(title=title)
    plot!(gridpts,IE[1,:],linewidth=3,label = fmet[2], yscale=:log10, legend=:topright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(gridpts,IE[2,:],linewidth=3,label = fmet[3], yscale=:log10)
    plot!(gridpts,IE[3,:],linewidth=3,label = fmet[4], yscale=:log10)
    savefig("graphs/IntError_$plotname")
end


## Grid spacing
##

# # No need for ScaledInterpolations if not using PolyRange
#grid = range(0,1,length=n)
# grid = a + (b-a)*grid^θ PolYnomial scaling
# a = 0.01, b = 0.5, θ = 2, interval = [0,1], number of points = 6
#    ci = 0.01 .+ (0.50-0.01).*range(0,1;length=6).^2

function polygrid(n::Integer,a::Float64,b::Float64,θ::Float64)
    grid = collect(range(0,1; length=n))
    xi = a .+ (b-a).*grid.^θ
    return xi
end

polygrid(6,0.05,2.0,2.0)

function interp_err_SN(F::Array,I::Array)
    error = maximum(abs.(I.-F))/maximum(abs.(F))
end

maximum(abs.(interp.-u1.(xax)))/maximum(abs.(u1.(xax)))

interp_err_SN(u1.(xax),interp)

fnom = ["Log","Square Root","CES "*latexinline("σ=2"), "CES "*latexinline("σ=5"), "CES "*latexinline("σ=10")]
funs = [u1,u2,u4,u5,u6]
slopes = [u1p,u2p,u4p,u5p,u6p]
method = [newton, CubicNaturalEval, Schumaker]
curvature = [1.0,2.0,3.0,4.0]
nom = ["Log","SQR","CES_2", "CES_5", "CES_10"]
titlenom = ["Log","Square Root","CES σ=2", "CES σ=5", "CES σ=10"]
#titlenom = ["Log_n=6","Square_Root n=6","CES σ=2 n=6", "CES σ=5 n=6", "CES σ=10 n=6"]
for i in 1:length(funs)
    plotname = nom[i]
    title = titlenom[i]
    name = fnom[i]
    fn = funs[i].(xax)
    IE = Array{Float64}(undef, 3, 4)
    for j in 1:length(method)
        met = method[j]
        for n in 1:length(curvature)
            curv = curvature[n]
            # Grid of nodes for interpolation
            global xi = polygrid(6,0.05,2.0,curvature[n]) ; # Collect makes it an array instead of a collection
            global yi = funs[i].(xi) # the corresponding y-coordinates
            # Interpolation
            if j ≤ 2
                global interp = map(z->method[j](xi,yi,z)[1],xax)
            else
                global si = slopes[i].(xi) #the corresponding slopes
                global interp=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
            end
            #IE[j,n] = interp_err_SN(fn,interp)
            IE[j,n] = sqrt(sum((fn .- interp).^2))
        end
    end
    global fmet = ["$name","Newton Basis Polynomials","Natural Cubic Splines","Shape-preserving Schumaker Splines"]
    mdtable(IE, side=fmet, head=latexinline(["θ=1","θ=2","θ=3","θ=4"]),fmt = FancyNumberFormatter(4))|> print
    # Plot
    gr()
    #pyplot()
    plot(title=title, xticks=curvature)
    plot!(curvature,IE[1,:],linewidth=3,label = fmet[2], yscale=:log10,  legend=:topleft,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(curvature,IE[2,:],linewidth=3,label = fmet[3],yscale=:log10,)
    plot!(curvature,IE[3,:],linewidth=3,label = fmet[4],yscale=:log10,)
    xlabel!("Curvature")
    savefig("graphs/IntError_Curv_$plotname")
end


fnom = ["Log","Square Root","CES "*latexinline("σ=2"), "CES "*latexinline("σ=5"), "CES "*latexinline("σ=10")]
funs = [u1,u2,u4,u5,u6]
slopes = [u1p,u2p,u4p,u5p,u6p]
method = [newton, CubicNaturalEval, Schumaker]
curvature = [1.0,1.5,2.0,3.0]
nom = ["Log","SQR","CES_2", "CES_5", "CES_10"]
titlenom = ["Log","Square Root","CES σ=2", "CES σ=5", "CES σ=10"]
#titlenom = ["Log_n=6","Square_Root n=6","CES σ=2 n=6", "CES σ=5 n=6", "CES σ=10 n=6"]
for i in 1:length(funs)
    plotname = nom[i]
    title = titlenom[i]
    name = fnom[i]
    fn = funs[i].(xax)
    IE = Array{Float64}(undef, 3, 4)
    local interp = Array{Float64}(undef,1000,3)
    for n in 1:length(curvature)
        curv = curvature[n]
        for j in 1:length(method)
            met = method[j]
            # Grid of nodes for interpolation
            global xi = polygrid(6,0.05,2.0,curvature[n]) ; # Collect makes it an array instead of a collection
            global yi = funs[i].(xi) # the corresponding y-coordinates
            # Interpolation
            if j ≤ 2
                global interp[:,j] = map(z->method[j](xi,yi,z)[1],xax)
            else
                global si = slopes[i].(xi) #the corresponding slopes
                global interp[:,j]=map(z->Schumaker(xi,yi,si,z)[1],xax) # Interpolating poly for the data
            end
            #IE[j,n] = interp_err_SN(fn,interp)
            #IE[j,n] = sqrt(sum((fn .- interp).^2))
        end
        global fmet = ["$name","Newton Basis Polynomials","Natural Cubic Splines","Shape-preserving Schumaker Splines"]
        # Plot
        gr()
        plot(title="Interpolation n=6 - $title, θ=$curv")
        plot!(log.(xax),fn,linewidth=3,label = "Function: $title",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(log.(xax),interp[:,1],linewidth=3,label=fmet[2], linestyle=:dash)
        plot!(log.(xax),interp[:,2],linewidth=3,label=fmet[3],linestyle=:dot)
        plot!(log.(xax),interp[:,3],linewidth=3,label=fmet[4], linestyle=:dashdot)
        plot!(log.(xi),yi,linetype=:scatter,marker=(:circle,6),markercolor=:blue,label = "Data")
        savefig("graphs/$plotname $curv.png")
    end

    # mdtable(IE, side=fmet, head=latexinline(["θ=1","θ=2","θ=3","θ=4"]),fmt = FancyNumberFormatter(4))|> print
    # # Plot
    # gr()
    # #pyplot()
    # plot(title=title, xticks=curvature)
    # plot!(curvature,IE[1,:],linewidth=3,label = fmet[2], yscale=:log10,  legend=:topleft,foreground_color_legend = nothing,background_color_legend = nothing)
    # plot!(curvature,IE[2,:],linewidth=3,label = fmet[3],yscale=:log10,)
    # plot!(curvature,IE[3,:],linewidth=3,label = fmet[4],yscale=:log10,)
    # xlabel!("Curvature")
    # savefig("graphs/IntError_Curv_$plotname")
end




# function curv(n::Float64,a::Float64,b::Float64,F::function,T::function)
#     error = 1.0
#     θ = 1
#     xi = polygrid(n,a,b,θ)
#     yi = F.(xi)
#     fn = F.(xax)
#     interp = map(z->T(xi,yi,z),xax)
#     while error > 10E-5
