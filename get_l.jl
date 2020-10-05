function get_l(k,kp,p::Par)
    A = range(0,1;length=100)
    B = [utility(k,kp,A[i],p) for i in 1:100]
    index = findmax(B)[2]
    l = A[index]
    return l
end

get_l.(k_grid_20,g20,fill(p,20))

g20 = rand(k_grid_20,20)

gl = get_l.(k_grid_20,g20,fill(p,20))

gml = get_me_l.(k_grid_20,g20,fill(p,20))

plot(gl,gml)

plot(k_grid_20,gl)
plot!(k_grid_20,gml)
