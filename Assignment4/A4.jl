# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")
# mkpath("Assignment4")
cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment4")
# mkpath("graphs")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
#Pkg.add("Optim")
using Optim #  # https://julianlsolvers.github.io/Optim.jl/stable/
    using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
# Call Scaled Interpolation Functions
    include("./Scaled_Interpolation_Functions.jl")
## Minimum Bracketing Function
    # Minimum Bracketing Function
        # This function comes from Numerical Recipes, Section 10.1
        # Input: numbers a,b that initiliaze bracket and objective function F
            # Optional: Bounds for the brackeet so that a,b,c ∈ [x_min,x_max]
        # Output: numbers (a,b,c) and images (fa,fb,fc) such that a<b<c and fb<fa,fc
        function mnbrak(a,b,F::Function,x_min=-Inf,x_max=Inf)
            # Auxiliariy variables
            Tiny    = 1E-20
            G_limit = 100
            # Define golden ratio
            γ = 1.618034
            # Evaluate function at end points
            fa, fb = F(a), F(b)
            # Swap a and b so that we can go downhill in the direction from a to b.
                # This way we know for certain that fb<fa, we only need to find c such that fb<fc
            if fb>fa
            a , b = b , a
            fa,fb = fb, fa
            end
            # Guess for value c expanding bracket away from b -> (a,b,c) or (c,b,a)
                # Check bounds
            c  = max( min( b + γ*(b-a) , x_max ) , x_min)
            fc = F(c)
            # While fb>fc we need to keep on bracketing
            while fb>fc
                # Compute u (new candidate) by parabolic extrapolation from a, b, c.
                    # TINY is used to prevent any possible division by zero.
                r = (b-a)*(fb-fc)
                q = (b-c)*(fb-fa)
                u = min(max(b-((b-c)*q-(b-a)*r)/(2*sign(q-r)*max(abs(q-r),Tiny)),x_min),x_max)
                u_lim = min(max(b+G_limit*(c-b),x_min),x_max)
                # Test various cases for new candidate
                if ((b-u)*(u-c) > 0) # Parabolic u is between b and c
                    fu=F(u)
                    if (fu < fc) # Got a minimum between b and c so bounds are (b,u,c)
                        a, fa, b, fb = b, fb, u, fu
                        break
                    elseif (fu > fb) # Got a minimum between a and u so bounds are (a,b,u)
                        c, fc = u, fu
                        break
                    end
                    # If you got here is because candidate u failed
                    # Get new candidate by expanding interval with golden ratio
                    # Check bounds
                    u  = max(min( c+γ*(c-b) , x_max ),x_min)
                    fu = F(u)
                elseif ((c-u)*(u-u_lim) > 0.0) # Parabolic u is between c and its limit (ulim)
                    fu=F(u)
                    if (fu < fc) # Drop c and replace it with u, get new u with golden expansion
                        b, c, fb, fc = c, u, fc, fu
                        u  = max(min( c+γ*(c-b) , x_max ),x_min)
                        fu = F(u)
                    end
                elseif ((u-u_lim)*(u_lim-c) >= 0.0) # U is above its limit, reign it in!
                    u  = u_lim
                    fu = F(u)
                else # Nothing worked, use golden expansion
                    u  = max(min( c+γ*(c-b) , x_max ),x_min)
                    fu = F(u)
                end
                # If no break then forget the oldest point and go onto next iteration
                    # This means assigning b->a, c->b, u-> and forgetting about a
                a, b, c, fa, fb, fc = b, c, u, fb, fc, fu
            end
            # Return solution once out of the loop
                # Swap a and c so that a<b<c
                if a>c
                a , c  = c, a
                fa, fc = fc, fa
                end
            # Minimum bracketed in [a,c] with intermediate point b so that fb<fa,fc
            # println("Minimum bracketed in [$a,$c] with intermediate point b=$b \n Function values: F(a)=$fa, F(b)=$fb, F(c)=$fc")
            return a,c,b,fa,fc,fb
        end



## VFI- Continuous

# defining parameters

# Paramters
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        z::Float64 = 1    ; # Productivity
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.98 ; # Discount factor
        σ::Float64 = 2 ;
        η::Float64 = 1 ;
        δ::Float64 = 0.05 ;
        l::Float64 = 0.4 ;
        M = ((z*α*β)/(1-β+β*δ))^(1/(1-α))
        N = z*M^α-δ*M
        χ = (z*(1-α)*M^α)/(l^(σ+η)*N^σ)
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9   ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Minimum consumption for numerical optimization
        c_min::Float64    = 1E-16;
        l_min::Float64    = 1E-16;
        l_max::Float64    = 1;
    end

    p = Par()

    #Steady state constants
    # @unpack z, α, β, σ, η, l, δ = p
    # global M = ((z*α*β)/(1-β+β*δ))^(1/(1-α))
    # global N = z*M^α-δ*M
    # # Getting the right χ for l_ss= 0.4
    # global χ = (z*(1-α)*M^α)/(l^(σ+η)*N^σ)
## Grid
    function Make_K_Grid(n_k,θ_k,p::Par,scale_type="Poly")
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        # Get k_grid
        if θ_k≠1
            if scale_type=="Poly"
            k_grid = PolyRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
            elseif scale_type=="Exp"
            # k_grid = ExpRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
            else
            error("scale_type must be either Poly or Exp")
            end
        else
        k_grid = range(1E-5,2*k_ss,length=n_k)
        end
        # Return
        return k_grid
    end

## Model Parameters

# Generate structure of model objects
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model paramters in their own structure
        # Grids
        θ_k::Float64    = 1     # Curvature of k_grid
        n_k::Int64      = 20    # Size of k_grid
        n_k_fine::Int64 = 1000  # Size of fine grid for interpolation
        k_grid          = Make_K_Grid(n_k,θ_k,p)    # k_grid for model solution
        k_grid_fine     = Make_K_Grid(n_k_fine,1,p) # Fine grid for interpolation
        # Value and policy functions
        V         = Array{Float64}(undef,n_k)       # Value Function
        G_kp      = Array{Float64}(undef,n_k)       # Policy Function
        G_c       = Array{Float64}(undef,n_k)       # Policy Function
        G_l       = Array{Float64}(undef,n_k)       # Policy Function
        V_fine    = Array{Float64}(undef,n_k_fine)  # Value Function on fine grid
        G_kp_fine = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
        G_l_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
        # Anaytical Solutions
        V_a       = Array{Float64}(undef,n_k_fine)  # Analytical Value Function on fine grid
        G_kp_a    = Array{Float64}(undef,n_k_fine)  # Analytical Policy Function on fine grid
        G_c_a     = Array{Float64}(undef,n_k_fine)  # Analytical Policy Function on fine grid
        Euler     = Array{Float64}(undef,n_k_fine)  # Errors in Euler equation
    end

    M = Model()


## Utility function

    function utility(k,kp,l,p::Par)
        @unpack z, α, η, σ, χ, c_min, l_min, δ= p
        y = (a,b) -> z*(a^α)*b^(1-α)+(1-δ)*a
        c = (a,b,d) -> y(a,b)-d
        lb = (a) -> (χ*a^(η+1)/(η+1))
        if l ≤ l_min
            lab = lb(l_min)
            if c(k,kp,l_min) ≤ c_min
                cons = c_min^(1-σ)/(1-σ)
            else
                cons = c(k,kp,l_min)^(1-σ)/(1-σ)
            end
        else
            lab  = lb(l)
            if c(k,kp,l) ≤ c_min
                cons = c_min^(1-σ)/(1-σ)
            else
                cons = c(k,kp,l)^(1-σ)/(1-σ)
            end
        end
        return cons-lab
    end

    # Derivative of utility function wrt labor
    function d_utility_l(k,kp,l,p::Par)
        @unpack α,δ,σ,z,η,c_min,χ = p
        c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp
        d_u = 0
        if c>c_min
            d_u = c^(-σ)
        else
            d_u = c_min^(-σ)
        end
        return d_u*z*(k^α)*(1-α)*(l^(-α))-χ*(l^η)
    end

## Steady State values fn
    # Steady state values (funciton)
    function SS_values(p::Par)
        # This function takes in parameters and provides steady state values
        # Parameters: productivity (z), returns to scale (a) and discount factor (b)
        # Output: values for capital, production, consumption, rental rate, wage
        @unpack z, α, β, l, M, N = p
        l_ss = l
        k_ss = M*l_ss
        c_ss = N*l_ss
        y_ss = z*l_ss*M^α
        r_ss = α*z*M^(α-1)
        w_ss = (1-α)*z*M^α
        return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
    end

    # Test steady state function
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    println(" ")
    println("------------------------")
    println(" Steady State variables")
    println("   Quantities: l = $l_ss; k = $k_ss; y = $y_ss; c = $c_ss;")
    println("   Prices:     r = $r_ss; w = $w_ss;")
    println("------------------------")
    println(" ")

    # Path to steady state
    function path(k,kp,l,M::Model)
        # This function takes in k, k', l and parameters and provides
        # Output: values for capital, production, consumption, rental rate, wage
        @unpack p, n_k = M
        @unpack z, α, β, δ = p
        w = (a,b) -> (1-α)*(a^α)*(b^(-α))
        r = (a,b) -> α*(a^(a-1))*(b^(1-α))
        y = (a,b) -> (a^α)*b^(1-α)+(1-δ)*a
        c = (a,b,d) -> y(a,b)-d
        wage = w.(k,l)
        int = r.(k,l)
        cons = c.(k,kp,l)
        output = y.(k,l)
        labor = get_l.(k,kp,fill(M,length(kp)))
        return wage, int, cons, labor, output
    end


## Optimal labot- Continuous

    # Function to find optimal labor
    function get_l(k,kp,M::Model)::Float64
        @unpack p = M
        L(x) = -utility(k,kp,x,p)
        result = maximize(L,0,1)
        converged(result) || error("Failed to converge in $(iterations(result)) iterations")
        l = maximizer(result)
        return l
    end

## Euler Errors

    # # Euler Equation
    # function Euler_Error(k,kp,kpp,l,lp,p::Par)
    #     # Return percentage error in Euler equation
    #     @unpack z, α, β, σ = p
    #     LHS = 1 ./(z.*k.^α.*l.^(1 .-α).-kp).^σ
    #     RHS = β.*α.*z.*kp.^(α.-1).*lp.^(1 .-α)./(z.*kp.^α.*lp.^(1 .-α).-kpp).^σ
    #     return (RHS/LHS-1)*100
    # end

    # Function that returns the percentage error in the euler equation
    function Euler_Error(k,kp,kpp,l,lp,p::Par)
        # Return percentage error in Euler equation
        @unpack z, α, β, σ, δ = p
        k_in = k.>0
        kp_in = kp.>0
        kpp_in = kpp.>0
        l_in = l.>0
        lp_in = lp.>0
        k = k[k_in.*kp_in.*kpp_in.*l_in.*lp_in]
        kp = kp[k_in.*kp_in.*kpp_in.*l_in.*lp_in]
        kpp = kpp[k_in.*kp_in.*kpp_in.*l_in.*lp_in]
        l = l[k_in.*kp_in.*kpp_in.*l_in.*lp_in]
        lp = lp[k_in.*kp_in.*kpp_in.*l_in.*lp_in]
        LHS = (z.*(k.^α).*(l.^(1-α)).+(1-δ).*k.-kp).^(-σ)
        RHS = β.*(α.*z.*((lp./kp).^(1-α)).+(1-δ)).*((z.*(kp.^α).*(lp.^(1-α)).+(1-δ).*kp.-kpp).^(-σ))
        return (RHS./LHS.-1).*100
    end

    # Euler error by number of points in grid
    function Euler_grid(M::Model)
        @unpack n_k_fine, k_grid_fine, p = M
        # Interpolation of G_kp
        G_kp_ip = ScaledInterpolations(M.k_grid,M.G_kp, BSpline(Cubic(Line(OnGrid()))))
        G_l_ip = ScaledInterpolations(M.k_grid,M.G_l, BSpline(Cubic(Line(OnGrid()))))
        # Euler error of numerical policy function on grid
        Euler = zeros(n_k_fine)
        for i=1:n_k_fine
            k   = k_grid_fine[i]
            kp  = G_kp_ip(k)
            kpp = G_kp_ip(kp)
            l   = G_l_ip(k)
            lp  = G_l_ip(kp)
            Euler[i] = Euler_Error(k,kp,kpp,l,lp,p)
        end
        # Return
        return Euler
    end

## VFI fixed point
    function VFI_Fixed_Point(T::Function,M::Model)
        # Unpack model structure
        @unpack p, n_k, θ_k, k_grid = M
        # VFI paramters
        @unpack max_iter, dist_tol = p
        # Initialize variables for loop
        V_old  = zeros(n_k)     ; # Initialize value function, here I just do 0, not the best option
        # V_old  = utility.(collect(k_grid),zeros(n_k),p) ; # Start at utility with zero savings
        V_dist = 1              ; # Initialize distance
        println(" ")
        println("------------------------")
        println("VFI - n_k=$n_k - θ_k=$θ_k")
        for iter=1:max_iter
            # Update value function
            V_new, G_kp, G_c, G_l = T(Model(M,V=copy(V_old)))
                # println("T(V) = $V_new")
                # println("  V  = $V_old")
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,100)==0
                println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
            end
            # Check convergence and return results
            if V_dist<=dist_tol
                println("VFI - n_k=$n_k - θ_k=$θ_k")
                println("Iterations = $iter and Distance = ",100*V_dist,"%")
                println("------------------------")
                println(" ")
                # Interpolate to fine grid
                V_ip = ScaledInterpolations(M.k_grid,V_new, BSpline(Cubic(Line(OnGrid()))))
                    V_fine = V_ip.(collect(M.k_grid_fine))
                G_kp_ip = ScaledInterpolations(M.k_grid,G_kp, BSpline(Cubic(Line(OnGrid()))))
                    G_kp_fine = G_kp_ip.(collect(M.k_grid_fine))
                G_c_ip = ScaledInterpolations(M.k_grid,G_c, BSpline(Cubic(Line(OnGrid()))))
                    G_c_fine = G_c_ip.(collect(M.k_grid_fine))
                G_l_ip = ScaledInterpolations(M.k_grid,G_l, BSpline(Cubic(Line(OnGrid()))))
                    G_l_fine = G_l_ip.(collect(M.k_grid_fine))
                # Update model
                M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,G_l=G_l,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine,G_l_fine=G_l_fine)
                # Euler error
                Euler = Euler_Error(M.k_grid_fine,G_kp_fine,G_kp_ip.(collect(G_kp_fine)),G_l_fine,G_l_ip.(collect(G_kp_fine)),p)
                # Update model
                M = Model(M; Euler=Euler)
                # Return results
                return M::Model
            end
        end
        # If loop ends there was no convergence -> Error!
        error("Error in VFI - Solution not found")
    end

## Bellman Max
# Bellman operator - Continuous Choice- multivariable optimization
function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack z, α, β, c_min, δ, l_min = p
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # println("Initial V = $V")
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
        # println("Interpolation test: Vp(k_min)=",Vp(k_grid[1])," Vp(k/2)=",Vp((k_grid[1]+k_grid[end])/2))
    # Define the derivative (might not be used, depends on algorithm)
        # For many methods we don't need to use ForwardDiff, we have direct measures,
        # for example for cubic splines. I am using ForwardDiff as an example
    # dVp(x) = ForwardDiff.derivative(Vp,x)

        # println("Interpolation test: Obj_Fun(k_min)=",Obj_Fun(k_grid[1])," Obj_Fun(k/2)=",Obj_Fun((k_grid[1]+k_grid[end])/2))
    get_kp_max(k) = min(z*(k^α)*(1^(1-α)) + (1-δ)*k - c_min , k_grid[end])
    kp_min = k_grid[1]
    l_max = 1.0
    Obj_Fun(k,x1,x2) = -utility(k,x1,x2,p) - β*Vp(x1)
    for i = 1:n_k
        # Min and max kp given current k
        kp_max = get_kp_max(k_grid[i])
        k_grid[i] < k_ss || (kp_min = k_ss)
        k_grid[i] > k_ss || (kp_max = k_ss)
            # Note: we can be smarter here.
            # We know that if k_grid[i]<k_ss then the kp_max is actually k_ss
            # Similarly, if k_grid[i]>k_ss then kp_min = k_ss
        # Check if lower bound binds
        k_0, l_0, k_corner, l_corner = corner(k_grid[i],Vp,Model(n_k=100))
        #x0 = [k_gr, l_0]
            # No need in this problem because without capital there is no production
            # also no consumption. So value would be -Inf.
            # I am leaving two ways of getting the derivative
        # dObj_min = ForwardDiff.derivative(Obj_Fun,kp_min)
        # dObj_min = -d_utility(k_grid[i],kp_min,p) - β*dVp(kp_min)
        # dObj_min = -1
        if k_corner
            # Define objective function: Right hand side of Bellman operator
                # Optimizer will minimize so we need the negative of the function
             V[i]    = Obj_Fun(k_grid[i],k_0,l_0)
             G_kp[i] = k_0
             G_l[i] = l_0
             println("\n corner solution for capital")
             println("k=$k_0 ; l=$l_0; VF=$(V[i])")
        elseif l_corner
            #println("\n corner solution for labor")
            min_result = optimize(x->Obj_Fun(k_grid[i],x,l_0),kp_min,kp_max)
            # Record results
                V[i]     = -min_result.minimum
                G_kp[i]  = min_result.minimizer
                G_l[i]   = l_0
                println("\n corner solution for labor")
                println("k=$(min_result.minimizer); l=$l_0; VF=$(V[i])")
        # # Check if upper bound binds
        #     # This can bind but shouldn't I will leave code to check but comment it
        #     # I am leaving two ways of getting the derivative
        # # dObj_max = ForwardDiff.derivative(Obj_Fun,kp_max)
        # # dObj_max = -d_utility(k_grid[i],kp_max,p) - β*dVp(kp_max)
        # dObj_max = 1
        # if dObj_max<0
        #     V[i]    = -Obj_Fun(kp_max)
        #     G_kp[i] = kp_max
        # else
        # Bracket the solution for k=k_grid[i]
            # In this case this is not necessary, we know the min is bracketed by [kp_min and kp_max]
            # We can still try to use mnbrak to get a better bracket, but in this case the algorithm
            # stops in the initial evaluation, it just verifies that we do have a bracket
            #kp_min, kp_max = mnbrak(kp_min,(kp_max+kp_min)/2,Obj_Fun,kp_min,kp_max)
        else
            # Maximize
            upper = [kp_max,1]
            lower = [kp_min,0]
            init_val = [(kp_min+kp_max)/2,l_ss]
            #inner_optimizer = NelderMead()
            min_result = Optim.optimize(x->Obj_Fun(k_grid[i],x[1],x[2]),lower,upper,init_val,NelderMead(), Optim.Options(iterations=10^6))
                #min_result = optimize(Obj_Fun,kp_min,kp_max)
                #println("maximum = $(-min_result.minimum) with argmin = $(min_result.minimizer) in "*"$(min_result.iterations) iterations")
            # Check result
                converged(min_result) || error("Failed to solve Bellman max in $(iterations(min_result)) iterations")
            # Record results
                V[i]     = -min_result.minimum
                G_kp[i]  = min_result.minimizer[1]
                G_l[i]   = min_result.minimizer[2] #get_l(k_grid[i],G_kp[i],M)
                # println("   Maximization result k_grid[$i] - kp=",min_result.minimizer," V(k)=",min_result.minimum," in $(iterations(min_result)) iterations")
        end # Upper bound if  
        # Lower bound if
    end # loop of k_grid[i]
    # Fill in policies
    y = (a,b) -> (a^α)*b^(1-α)+(1-δ)*a
    c = (a,b,d) -> y(a,b)-d
    G_c = c.(k_grid,G_kp,G_l)
    # Return Results
        # println("T(V) = $V")
    return V, G_kp, G_c, G_l
end

## Graphs
# Graphs
function VFI_Graphs(M::Model,VFI_Type)
    gr()
    # Value Function Analytical vs 200
        plot(M.k_grid_fine,M.k_grid_fine,linewidth=4,label = "Analytical Solution",title="Value Function",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot(M.k_grid,M.V,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        plot!(M.k_grid_fine,M.V_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        xlabel!("Capital")
        ylabel!("Value")
        savefig("./graphs/VFI_"*VFI_Type*"_V_$(M.n_k)_$(M.θ_k).png")
    # Capital Policy Function Analytical vs 200
        plot(M.k_grid_fine,M.k_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title="Policy Function - K",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        #plot!(M.k_grid_fine,M.G_kp_a,linewidth=3,label = "Analytical Solution",title="Policy Function - K",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid,M.G_kp,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        plot!(M.k_grid_fine,M.G_kp_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        xlabel!("Capital")
        ylabel!("Capital")
    savefig("./graphs/VFI_"*VFI_Type*"_G_kp_$(M.n_k)_$(M.θ_k).png")
        # Euler Error 200
        plot(M.k_grid_fine,zeros(M.n_k_fine),lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid_fine,M.Euler,linetype=:scatter,label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        xlabel!("Capital")
        ylabel!("Percentage Points")
    savefig("./graphs/VFI_"*VFI_Type*"_Euler_$(M.n_k)_$(M.θ_k).png")

    println("\n     Graphs Completed for VFI_$VFI_Type - n_k=$(M.n_k) - θ_k=$(M.θ_k) \n")
end

## Initial Guess

function corner(k,V::AbstractInterpolation,M::Model)
    @unpack n_k, k_grid = M
    l = fill(0.0,n_k)
    V_aux = fill(0.0,n_k)
    for i in 1:n_k
        l[i] = get_l(k,k_grid[i],M)
        V_aux[i] = utility(k,k_grid[i],l[i],p)+V(k_grid[i])
    end
    index = findmax(V_aux)[2]
    k_0 = k_grid[index]
    l_0 = l[index]
    if k_0 == k_grid[end] || k_0 == k_grid[1]
        #print("Corner solution for capital \n")
        k_corner = true
    else
        k_corner = false
    end
    if l_0 == 1 || l_0 == 0
        #print("Corner solution for labor \n")
        l_corner = true
    else
        l_corner = false
    end
    return k_0, l_0, k_corner, l_corner
end

## Bellman operator univariate

# Bellman operator for the continuous choice capital only where labor is solved for with the intratemporal FOC with bisection
function T_univariate_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack β,α,z,δ,η,σ,c_min = p
    # Interpolate current iteration of value function v(k') so I can evaluate it at any k'
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid())))) # Impose monotonicity because I know it is increasing in capital
    # Define boundaries on capital tomorrow k' and on labor l
    l_min = 1E-16
    l_max = 1.0
    get_kp_max(k,l) = z*(k^α)*(l^(1-α)) + (1-δ)*k - c_min # Function because the max k' depends on the value of capital today k and labor l
    # Define function that finds optimal labor l given (k,k') and returns objective function
    function Obj_fn(k,kp,p::Par)
        min_result = optimize(x->d_utility_l(k,kp,x,p).^2,l_min,l_max,Brent())
        l = min_result.minimizer
        return -utility(k,kp,l,p) - β*Vp.(kp),l
    end
    # optimization of the objective function for each capital level in the grid
    # Initialization
    #V_new = zeros(n_k)
    kp_min = 0.0001
    # Maximize
    for i in 1:n_k
        kp_max = min(get_kp_max(k_grid[i],1.0),0.9999*k_grid[end])
        min_result = optimize(x->Obj_fn(k_grid[i],x,p)[1],kp_min,kp_max,Brent())
        # Check result
        #converged(min_result) || error("Failed to solve Bellman max for capital =" k_grid[i]" in $(iterations(min_result)) iterations")
        # Record results
        V[i] = -min_result.minimum
        G_kp[i] = min_result.minimizer
        G_l[i] = Obj_fn(k_grid[i],G_kp[i],p::Par)[2]
    end
    # Fill in policy for consumption
    G_c = z.*(collect(k_grid).^α).*(G_l.^(1-α)) .- G_kp
    # Return results
    return V, G_kp, G_c, G_l
end

## Solve

    # Execute VFI and plot graphs

    # Execute Numerical VFI - equally spaced grid
        @time M_20  = VFI_Fixed_Point(T_cts_max,Model(n_k=20))
        @time M_50  = VFI_Fixed_Point(T_cts_max,Model(n_k=50))
        @time M_100 = VFI_Fixed_Point(T_cts_max,Model(n_k=100))
        # Graphs
            VFI_Graphs(M_20,"cts_max")
            VFI_Graphs(M_50,"cts_max")
            VFI_Graphs(M_100,"cts_max")
    # Execute Numerical VFI - curved spaced grid
        @time M_20c  = VFI_Fixed_Point(T_cts_max,Model(n_k=20,θ_k=2.5))
        @time M_50c  = VFI_Fixed_Point(T_cts_max,Model(n_k=50,θ_k=2.5))
        @time M_100c = VFI_Fixed_Point(T_cts_max,Model(n_k=100,θ_k=3.5))
        # Graphs
            VFI_Graphs(M_20c,"cts_max")
            VFI_Graphs(M_50c,"cts_max")
            VFI_Graphs(M_100c,"cts_max")

    # Execute Numerical VFI - equally spaced grid
        @time M_20u  = VFI_Fixed_Point(T_univariate_max,Model(n_k=20))
        @time M_50u  = VFI_Fixed_Point(T_univariate_max,Model(n_k=50))
        @time M_100u = VFI_Fixed_Point(T_univariate_max,Model(n_k=100))
        # Graphs
            VFI_Graphs(M_20u,"cts_max")
            VFI_Graphs(M_50u,"cts_max")
            VFI_Graphs(M_100u,"cts_max")
    # Execute Numerical VFI - curved spaced grid
        @time M_20cu  = VFI_Fixed_Point(T_univariate_max,Model(n_k=20,θ_k=2.5))
        @time M_50cu  = VFI_Fixed_Point(T_univariate_max,Model(n_k=50,θ_k=2.5))
        @time M_100cu = VFI_Fixed_Point(TT_univariate_max,Model(n_k=100,θ_k=2.5))
        # Graphs
            VFI_Graphs(M_20cu,"cts_max")
            VFI_Graphs(M_50cu,"cts_max")
            VFI_Graphs(M_100cu,"cts_max")

plot(M_20.k_grid_fine,M_20.G_kp_fine)
plot(M_20.k_grid,M_20.G_l)
