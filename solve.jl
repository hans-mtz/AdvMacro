#Pkg.add("TimerOutputs")
using TimerOutputs
#Pkg.add("Latexify")
using Latexify

 #cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro/Assignment2")

reset_timer!()

    # Execute Numerical VFI
    @timeit "Plain VFI n_k=20" V_20, G_kp_20, G_l_20, k_grid_20 = Solve_VFI_loop(20,p)
    @timeit "Plain VFI n_k=50"  V_50, G_kp_50, G_l_50, k_grid_50 = Solve_VFI_loop(50,p)
    @timeit "Plain VFI n_k=100"  V_100, G_kp_100, G_l_100, k_grid_100 = Solve_VFI_loop(100,p)
    #@timeit "Plain VFI n_k=200"  V_200, G_kp_200, G_l_200,k_grid_200 = Solve_VFI_loop(200,p)
    #@timeit "Plain VFI n_k=500" V_500, G_kp_500, G_l_500, k_grid_500 = Solve_VFI_loop(500,p)

    gr()
    # Value Function Analytical vs 200
        plot(k_grid_20,V_20,linetype=:scatter,marker=(:diamond,3), msw=0 , label="VFI - n_k=20")
        plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=50")
        plot!(k_grid_100,V_100,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=100")
        #plot!(k_grid_500,V_500,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=500")
        xlabel!("Capital")
        ylabel!("Value")
        title!("Value Fn VFI")
        savefig("graphs/VFI_V")


    # Capital Policy Function Analytical vs 200
        plot(k_grid_100,k_grid_100,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
        plot!(k_grid_20,k_grid_20[G_kp_20],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=20")
        plot!(k_grid_100,k_grid_100[G_kp_100],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=100")
        #plot!(k_grid_500,k_grid_500[G_kp_500],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=500")
        xlabel!("Capital")
        ylabel!("Capital")
        title!("Capital VFI")
        savefig("graphs/VFI_cap")

##-----------------
##-----------------

Euler_20 = Euler_grid(k_grid_20,G_kp_20)
Euler_50 = Euler_grid(k_grid_50,G_kp_50)
Euler_100 = Euler_grid(k_grid_100,G_kp_100)
#Euler_500 = Euler_grid(k_grid_500,G_kp_500)

# Graphing Euler Error 200
plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
plot!(k_grid_20,Euler_20,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=20")
plot!(k_grid_100,Euler_100,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=100")
#plot!(k_grid_500,Euler_500,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=500")
xlabel!("Capital")
ylabel!("Percentage Points")
title!("Euler Error VFI")
savefig("graphs/VFI_Euler")


###-----
# Execute Numerical VFI
@timeit "VFI-HPI n_k=20"  V_20, G_kp_20, G_l_20, k_grid_20 = Solve_VFI_HPI(20,20,p)
@timeit "VFI-HPI n_k=50"  V_50, G_kp_50, G_l_50, k_grid_50 = Solve_VFI_HPI(50,50,p)
@timeit "VFI-HPI n_k=200" V_200, G_kp_200, G_l_200, k_grid_200 = Solve_VFI_HPI(200,200,p)
@timeit "VFI-HPI n_k=500" V_500, G_kp_500, G_l_500, k_grid_500 = Solve_VFI_HPI(500,500,p)

# Value Function Plot
    plot(k_grid_20,V_20,linetype=:scatter,marker=(:diamond,3), msw=0 , label="VFI - n_k=20")
    plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=50")
    plot!(k_grid_200,V_200,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=200")
    plot!(k_grid_500,V_500,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Value")
    title!("Value Function VFI-HPI")
    savefig("graphs/HPI_V")


# Capital Policy Function
    plot(k_grid_200,k_grid_200,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
    plot!(k_grid_20,k_grid_20[G_kp_20],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=20")
    plot!(k_grid_200,k_grid_200[G_kp_200],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=200")
    plot!(k_grid_500,k_grid_500[G_kp_500],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Capital")
    title!("Capital VFI-HPI")
    savefig("graphs/HPI_Cap")

    Euler_20 = Euler_grid(k_grid_20,G_kp_20)
    Euler_50 = Euler_grid(k_grid_50,G_kp_50)
    Euler_200 = Euler_grid(k_grid_200,G_kp_200)
    Euler_500 = Euler_grid(k_grid_500,G_kp_500)

    # Graphing Euler Error 200
    plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
    plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
    plot!(k_grid_20,Euler_20,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=20")
    plot!(k_grid_200,Euler_200,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=200")
    plot!(k_grid_500,Euler_500,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Percentage Points")
    title!("Euler Error VFI-HPI")
    savefig("graphs/HPI_Euler")

##----
##----

# Execute Numerical VFI
@timeit "VFI-MPB n_k=20" V_20, G_kp_20, G_l_20, k_grid_20 = Solve_VFI_MPB(20,p)
@timeit "VFI-MPB n_k=50" V_50, G_kp_50, G_l_50, k_grid_50 = Solve_VFI_MPB(50,p)
@timeit "VFI-MPB n_k=200" V_200, G_kp_200, G_l_200, k_grid_200 = Solve_VFI_MPB(200,p)
@timeit "VFI-MPB n_k=500" V_500, G_kp_500, G_l_500, k_grid_500 = Solve_VFI_MPB(500,p)

print_timer()

Euler_20 = Euler_grid(k_grid_20,G_kp_20)
Euler_50 = Euler_grid(k_grid_50,G_kp_50)
Euler_200 = Euler_grid(k_grid_200,G_kp_200)
Euler_500 = Euler_grid(k_grid_500,G_kp_500)

# Value Function Plot
    plot(k_grid_20,V_20,linetype=:scatter,marker=(:diamond,3), msw=0 , label="VFI - n_k=20")
    plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=50")
    plot!(k_grid_200,V_200,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=200")
    plot!(k_grid_500,V_500,linetype=:scatter,marker=(:diamond,3),msw=0 , label = "VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Value")
    title!("Value Function VFI-MPB")
    savefig("graphs/MPB_V")


# Capital Policy Function
    plot(k_grid_200,k_grid_200,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
    plot!(k_grid_20,k_grid_20[G_kp_20],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=20")
    plot!(k_grid_200,k_grid_200[G_kp_200],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=200")
    plot!(k_grid_500,k_grid_500[G_kp_500],linetype=:scatter,marker=(:diamond,3),msw=0 ,label = "VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Capital")
    title!("Capital VFI-MPB")
    savefig("graphs/MPB_Cap")


    # Graphing Euler Error 200
    plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
    plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=50")
    plot!(k_grid_20,Euler_20,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=20")
    plot!(k_grid_200,Euler_200,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=200")
    plot!(k_grid_500,Euler_500,linetype=:scatter,marker=(:diamond,3),msw=0 ,label="VFI - n_k=500")
    xlabel!("Capital")
    ylabel!("Percentage Points")
    title!("Euler Error VFI-MPB")
    savefig("graphs/MPB_Euler")
