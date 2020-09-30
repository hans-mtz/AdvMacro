# setting working directory
# cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")
# mkpath("Assignment2")
cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment2")
# mkpath("Figures")


# Pkg.add("Parameters")
using Plots
using Parameters
# Pkg.add("Roots")

Pkg.add("IntervalArithmetic")
Pkg.add("IntervalRootFinding")
using Roots

# defining parameters

# Paramters
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        z::Float64 = 1    ; # Productivity
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.98 ; # Discount factor
        σ::Float64 = 2
        η::Float64 = 1
        δ::Float64 = 1
        χ::Float64
        l::Float64
        M = ((z*α*β)/(1-β+β*δ))^(1/(1-α))
        N = z*M^(α)-δ*M
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
    end

    # Allocate paramters to object p for future calling

    p = Par(l= 0.4, χ =1.0)

# Getting the right χ
function get_me_chi(l,p::Par)
    @unpack z, α, β, M, N, σ, η = p
    f(x) = x*(N^σ)*(l^(σ+η))-(z*(1-α)*M^α)
    chi = find_zero(f, 1E-9)
    return chi
end

chi = get_me_chi(0.4,p)

p = Par(l= 0.4, χ=40.2)

# Getting the right χ

using IntervalArithmetic, IntervalRootFinding

function get_me_l(k,kp,p::Par)
    @unpack z, α, β, M, N, σ, η, χ = p
    f(l) = z*(1-α)*(k^α)/((z*(k^α)*(l^(1-α))-kp)^σ)-χ*(l^(η+α))
    l = find_zero(f, (0.9), Order2(), atol = 0.1)
    #l = roots(f, 0.1..1)
    #l = fzero(f, 1, order=16)
    if l <= 1
    return l
    else
    return 1
    end
end

get_me_l(k_grid_20[1],k_grid_20[G_kp_20][1],p)


# Steady state values (funciton)
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, M, N, l = p
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


# Steady state values (funciton)
function path(k,kp,l,p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, M, N, l = p
    w = (a,b) -> (1-α)*(a^α)*(b^(-α))
    r = (a,b) -> α*(a^(a-1))*(b^(1-α))
    y = (a,b) -> (a^α)*b^(1-α)
    c = (a,b,d) -> y(a,b)-d
    l = get_me_l.(k,k_p,fill(p,length(k)))
    wage = w.(k,l)
    int = r.(k,l)
    cons = c.(k,kp,l)
    output = y.(k,l)
    labor = l
    return wage, int, cons, labor, output
end

# Utility function
function utility(k,kp,l,p::Par)
    @unpack z, α, χ, η, σ = p
    y = (a,b) -> (a^α)*b^(1-α)
    c = (a,b,d) -> y(a,b)-d
    utils = (c(k,l,kp).^(1-σ))/(1-σ)-(χ*l.^(η+1)/(η+1))
    return utils
end



# VFI with Grid
    # Grid
    function Make_K_Grid(n_k,p::Par)
        # Get SS
        l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        # Get k_grid
        k_grid = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
        # Return
        return k_grid
    end

    # Solve VFI with grid search and loops
    function VFI_grid(T::Function,k_grid,p::Par)
        # VFI paramters
        @unpack max_iter, dist_tol = p
        # Initialize variables for loop
        n_k    = length(k_grid) ; # Number of grid nodes
        V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
        V_dist = 1              ; # Initialize distance
        println(" ")
        println("------------------------")
        println("VFI - Grid Search - n_k=$n_k")
        for iter=1:max_iter
            # Update value function
            V_new, G_kp, G_l = T(V_old)
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            iter  += 1
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,100)==0
                println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
            end
            # Check convergence and return results
            if V_dist<=dist_tol
                println("VFI - Grid Search - n_k=$n_k")
                println("Iterations = $iter and Distance = ",100*V_dist,"%")
                println("------------------------")
                println(" ")
                return V_new, G_kp, G_l
            end
        end
        # If loop ends there was no convergence -> Error!
        error("Error in VFI - Grid Search - Solution not found")
    end

#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI - Grid Search - Loops


# Define function for Value update and policy functions
    function T_grid_loop(V_old,k_grid,p::Par)
        @unpack z, α, β = p
        n_k  = length(k_grid)
        V    = zeros(n_k)
        G_kp = fill(0,n_k)
        G_l = fill(0,n_k)
        for i = 1:n_k
            V_aux = zeros(n_k) ; # Empty vector for auxiliary value of V(i,j)
            G_l_aux = zeros(n_k)
            for j = 1:n_k
                # Evaluate potential value function for combinations of
                # current capital k_i and future capital k_j
                l = get_me_l(k_grid[i],k_grid[j],p)
                G_l_aux[i] = l
                V_aux[j] = utility(k_grid[i],k_grid[j],l,p) + β*V_old[j]
                #println(V_aux[j]," ",k_grid[i]," ",k_grid[j]," ",utility(k_grid[i],k_grid[j],z,a,b) + b*V_old[j])
            end
            # Choose maximum value given current capital k_i
            V[i], G_kp[i] = findmax(V_aux)
            G_l[i] = G_l_aux[G_kp[i]]
        end
        return V, G_kp, G_l
    end




# Solve VFI with grid search and loops
    function Solve_VFI_loop(n_k,p::Par)
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Solve VFI
        V, G_kp, G_l = VFI_grid(x->T_grid_loop(x,k_grid,p),k_grid,p)
        # Return Solution
        return V,G_kp, G_l, k_grid
    end

    # Execute Numerical VFI
    @time V_20, G_kp_20, G_l_20, k_grid_20 = Solve_VFI_loop(20,p)
    @time V_50, G_kp_50, G_l_50, k_grid_50 = Solve_VFI_loop(50,p)
    @time V_200, G_kp_200, G_l_200,k_grid_200 = Solve_VFI_loop(200,p)
    @time V_1000, G_kp_1000, G_l_1000, k_grid_1000 = Solve_VFI_loop(1000,p)

    gr()
    # Value Function Analytical vs 200
        plot(k_grid_20,V_20,linetype=:scatter,marker=(:diamond,3), msw=0 , label="VFI - n_k=50")
        plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=20")
        plot!(k_grid_200,V_200,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=200")
        plot!(k_grid_1000,V_1000,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=1000")
        xlabel!("Capital")
        ylabel!("Value")
        title!("Value Fn VFI")


    # Capital Policy Function Analytical vs 200
        plot(k_grid_200,k_grid_200,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
        plot!(k_grid_20,k_grid_20[G_kp_20],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=20")
        plot!(k_grid_200,k_grid_200[G_kp_200],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=200")
        plot!(k_grid_1000,k_grid_1000[G_kp_1000],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=1000")
        xlabel!("Capital")
        ylabel!("Capital")
        title!("Capital VFI")


    # Euler Equation
    function Euler_Error(k,kp,kpp,p::Par)
        # Return percentage error in Euler equation
        @unpack z, α, β, σ, l = p
        lp = l
        LHS = ((z*kp^α*lp^(1-α)-kpp)/(z*k^α*l^(1-α)-kp))^σ
        RHS = β*α*z*kp^(α-1)*lp^(1-α)
        return (RHS/LHS-1)*100
    end

    # Euler error by number of points in grid
    function Euler_grid(k_grid::AbstractArray,G_kp::AbstractArray)
        # Euler error of numerical policy function on grid
        n_k = length(k_grid)
        Euler = zeros(n_k)
        for i=1:n_k
            k   = k_grid[i]
            kp  = k_grid[G_kp[i]]
            kpp = k_grid[G_kp[G_kp[i]]]
            Euler[i] = Euler_Error(k,kp,kpp,p)
        end
        # Return
        return Euler
    end

    Euler_20 = Euler_grid(k_grid_20,G_kp_20)
    Euler_50 = Euler_grid(k_grid_50,G_kp_50)
    Euler_200 = Euler_grid(k_grid_200,G_kp_200)
    Euler_1000 = Euler_grid(k_grid_1000,G_kp_1000)

    # Graphing Euler Error 200
    plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
    plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
    plot!(k_grid_20,Euler_20,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=20")
    plot!(k_grid_20,Euler_200,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=200")
    plot!(k_grid_20,Euler_1000,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=1000")
    xlabel!("Capital")
    ylabel!("Percentage Points")
    title!("Euler Error VFI")

#-----------------------------------------------------------
    # VFI - Grid Search - Matrices - Howard's Policy Iteration

    # Define function for Value update and policy functions
        function HT_grid_mat(V_old,U_mat,k_grid,p::Par,n_H)
            @unpack z, α, β, H_tol = p
            # Get Policy Function
            n_k    = length(V_old)
            V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
            V_old  = V
            # "Optimal" U for Howard's iteration
                U_vec = U_mat[G_kp]
            # Howard's policy iteration
            # G_kp is a Cartesian Index
            for i=1:n_H
                V = U_vec .+ β*repeat(V_old',n_k,1)[G_kp]
                if maximum(abs.(V./V_old.-1))<=H_tol
                    break
                end
                V_old = V
            end
            # Recover Policy Functions
            G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
            # Return output
            return V, G_kp
        end


    # Solve VFI with Howard's policy iteration
        function Solve_VFI_HPI(n_H,n_k,p::Par)
            # Get Grid
            k_grid = Make_K_Grid(n_k,p)
            # Utility matrix
            U_mat = [utility(k_grid[i],k_grid[j],p) for i in 1:n_k, j in 1:n_k]
            # Solve VFI
            V, G_kp = VFI_grid(x->HT_grid_mat(x,U_mat,k_grid,p,n_H),k_grid,p)
            # Return Solution
            return V,G_kp, k_grid
        end

        # Execute Numerical VFI
        @time V_20, G_kp_20, k_grid_20 = Solve_VFI_HPI(20,20,p)
        @time V_50, G_kp_50, k_grid_50 = Solve_VFI_HPI(20,50,p)
        @time V_200, G_kp_200, k_grid_200 = Solve_VFI_HPI(20,200,p)
        @time V_1000, G_kp_1000, k_grid_1000 = Solve_VFI_HPI(20,1000,p)

        # Value Function Plot
            plot(k_grid_20,V_20,linetype=:scatter,marker=(:diamond,3), msw=0 , label="VFI - n_k=50")
            plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=20")
            plot!(k_grid_200,V_200,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=200")
            plot!(k_grid_1000,V_1000,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=1000")
            xlabel!("Capital")
            ylabel!("Value")
            title!("Value Function VFI-HPI")


        # Capital Policy Function
            plot(k_grid_200,k_grid_200,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
            plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
            plot!(k_grid_20,k_grid_20[G_kp_20],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=20")
            plot!(k_grid_200,k_grid_200[G_kp_200],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=200")
            plot!(k_grid_1000,k_grid_1000[G_kp_1000],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=1000")
            xlabel!("Capital")
            ylabel!("Capital")
            title!("Capital VFI-HPI")

            Euler_20 = Euler_grid(k_grid_20,G_kp_20)
            Euler_50 = Euler_grid(k_grid_50,G_kp_50)
            Euler_200 = Euler_grid(k_grid_200,G_kp_200)
            Euler_1000 = Euler_grid(k_grid_1000,G_kp_1000)

            # Graphing Euler Error 200
            plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
            plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
            plot!(k_grid_20,Euler_20,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=20")
            plot!(k_grid_20,Euler_200,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=200")
            plot!(k_grid_20,Euler_1000,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=1000")
            xlabel!("Capital")
            ylabel!("Percentage Points")
            title!("Euler Error VFI-HPI")
