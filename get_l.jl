function get_l(k,kp,p::Par)
    @unpack z, α, β, M, N, σ, η, χ = p
    function u(k,kp,l,p::Par)
        @unpack z, α, χ, η, σ = p
        y = (a,b) -> z*(a^α)*b^(1-α)
        c = (a,b,d) -> y(a,b)-d
        utils = (c(k,l,kp).^(1-σ)/(1-σ))-(χ*l.^(η+1)/(η+1))
        return utils
    end
    A = range(0,1;length=100)
    B = u.(k,kp,A,fill(p,length(k)))
    index = findmax(B)[2]
    l = A[index]
    return l
end

g20 = rand(k_grid_20,20)

gl = get_l.(k_grid_20,g20,fill(p,20))

gml = get_me_l.(k_grid_20,g20,fill(p,20))

plot(gl,gml)

plot(k_grid_20,gl)
plot!(k_grid_20,gml)
