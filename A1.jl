
# setting working directory
cd("/Volumes/SSD Hans/Github/MacroFall2020/AdvMacro")

#Defining vars
a = 1/3
b = 0.6
k_star = (a*b)^(1/(1-a))
r = 1/b

# Initializing array for loop
m, n = 100,2
K = fill(0.0, (m, n))

#defining policy fn
g = x -> a*b*(x^a)

# loop for capital
k = 0.8*k_star
for i in 1:m
    K[i, 1] = k
    K[i, 2] = g(K[i,1])
    global k = K[i, 2]
end

# defining fn's for the rest of the variables
w = x -> (1-a)x^a
y = x -> x^a
c = (p,q) -> y(p)-q
r = x -> a*(x^(a-1))
int = r.(K[1:100,2])
wage = w.(K[1:100,1])
output = y.(K[1:100,1])
cons = c.(K[1:100,1],K[1:100,2])

# importing packages for plotting
import Pkg; Pkg.add("Plots")
using Plots
gr()

# plotting capital
plot(1:10,fill(k_star,10), label="capital SS")
plot!(1:10,K[1:10,1], label="capital path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Capital")
savefig("graphs/capital")

# plotting output

plot(1:10,fill(y(k_star),10), label="output SS")
plot!(1:10,output[1:10], label="output path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Output")
savefig("graphs/output")

# plotting consumption

plot(1:10,fill(c(k_star,k_star),10), label="consumption SS")
plot!(1:10,cons[1:10], label="consumption path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Consumption")
savefig("graphs/consumption")

# plotting wage

plot(1:10,fill(w(k_star),10), label="wage SS")
plot!(1:10,wage[1:10], label="wage path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Wage")
savefig("graphs/wage")

# plotting capital rental rate
plot(1:10,fill(r(k_star),10), label="capital rental rate SS")
plot!(1:10,int[1:10], label="rate path", linestyle = :dot, legend = :topright)
xlabel!("t")
title!("Capital rental rate")
savefig("graphs/capitalr")

# Increase in productivity permanently by 5%

k_starp = (a*b*1.5)^(1/(1-a))
#defining policy fn with productivity increases
g1 = x -> a*b*1.05*(x^a)

# Initializing array for loop
K1 = fill(0.0, (m, n))

# loop for capital
k1 = k_starp
for i in 1:m
    K1[i, 1] = k1
    K1[i, 2] = g(K[i,1])
    global k1 = K[i, 2]
end

# defining fn's for the rest of the variables
w1 = x -> (1-a)*(x^a) #does not change
y1 = x -> 1.5*(x^a)
c1 = (p,q) -> y(p)-q #does not change
r1 = x -> a*(1.5*x^(a-1))
int_p = r1.(K1[1:100,2])
wage_p = w1.(K1[1:100,1])
output_p = y1.(K1[1:100,1])
cons_p = c1.(K1[1:100,1],K1[1:100,2])


# plotting capital
plot(1:10,fill(k_star,10), label="capital SS")
plot!(1:10,K1[1:10,1], label="capital path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Capital (productivity increase)")
savefig("graphs/capitalp")

# plotting output

plot(1:10,fill(y1(k_star),10), label="output SS")
plot!(1:10,output_p[1:10], label="output path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Output (productivity increase)")
savefig("graphs/outputp")

# plotting consumption

plot(1:10,fill(c1(k_star,k_star),10), label="consumption SS")
plot!(1:10,cons_p[1:10], label="consumption path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Consumption (productivity increase)")
savefig("graphs/consumptionp")

# plotting wage

plot(1:10,fill(w1(k_star),10), label="wage SS")
plot!(1:10,wage_p[1:10], label="wage path", linestyle = :dot, legend = :bottomright)
xlabel!("t")
title!("Wage (productivity increase)")
savefig("graphs/wagep")

# plotting capital rental rate
plot(1:10,fill(r1(k_star),10), label="capital rental rate SS")
plot!(1:10,int_p[1:10], label="rate path", linestyle = :dot, legend = :topright)
xlabel!("t")
title!("Capital rental rate (productivity increase)")
savefig("graphs/capitalrp")
