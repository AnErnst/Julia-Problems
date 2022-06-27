
import Pkg; Pkg.add("StaticArrays")
import Pkg; Pkg.add("Plots")


using StaticArrays
const SV = SVector{2,Float64}
const baseMatrix = SMatrix{2,2,Float64}
const GlueID = SVector{1,Int64}


abstract type AbstractParticle end

mutable struct Particle <: AbstractParticle
    pos::SV
    vel::SV
end
Particle(x0, y0, φ0) = Particle(SV(x0, y0), SV(cos(φ0), sin(φ0)))

abstract type Obstacle end

struct Wall <: Obstacle
    sp::SV
    ep::SV
    normal::SV
    A::baseMatrix
    
end

const Billiard = NTuple{N, Obstacle} where N

using LinearAlgebra: dot, normalize

"""
    collision(p::AbstractParticle, o::Obstacle) → t, cp
Find the collision (if any) between given particle and obstacle.
Return the time until collision and the estimated collision point `cp`.
"""
@inline function collision(p::Particle, w::Wall)
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision()
    else
        t = dot(w.sp - p.pos, n)/denom
        return t, p.pos + t * p.vel
    end
end

@inline nocollision() = (Inf, SV(0.0, 0.0))

normalvec(w::Wall, pos) = w.normal

function next_collision(p::AbstractParticle, bd)
    j, ct, cp = 0, Inf, SV(0.0, 0.0)
    for i in eachindex(bd)
        t, c = collision(p, bd[i])
        if t < ct
            j = i
            ct = t
            cp = c
        end
    end
    return j, ct, cp
end



propagate!(p::Particle, pos, t) = (p.pos = pos)

function resolvecollision!(p::AbstractParticle, o::Obstacle)
    n = normalvec(o, p.pos)
    p.vel = (p.vel - 2*dot(n, p.vel)*n)
    

end

"""
    bounce!(p, bd)
Evolve the particle for one collision (in-place).
"""
@inline function bounce!(p::AbstractParticle, bd)
    i::Int, tmin::Float64, cp::SV = next_collision(p, bd)
    if tmin != Inf
        propagate!(p, cp, tmin)
        resolvecollision!(p, bd[i])
        
    end
    return i, tmin, p.pos, p.vel
end

"""
    timeseries!(p::AbstractParticle, bd, n) -> xt, yt, t
Evolve the particle in the billiard `bd` for `n` collisions
and return the position timeseries `xt, yt` along with time vector `t`.
"""
function timeseries!(p::AbstractParticle, bd, n::Int)

    t = [0.0]; xt = [p.pos[1]]; yt = [p.pos[2]]; c = 0

    while c < n
        prevpos = p.pos; prevvel = p.vel
        i, ct = bounce!(p, bd)
        p.pos = (.5,.5)
        println(p.pos,p.vel)
        xs, ys = extrapolate(p, prevpos, prevvel, ct)  
        push!(t, ct) 
        append!(xt, xs)
        append!(yt, ys)
        c += 1
        
    end

    return xt, yt, t
end

extrapolate(p::Particle, prevpos, prevvel, ct) = p.pos

import Pkg; Pkg.add("PyPlot")
x, y, r = 1.0, 1.0, 0.3
A = [-1 0;0 -1]
sp = [0.0,1]; ep = [0.0, 0.0]; n = [1,0]
leftw = Wall(sp, ep, n, A)
sp = [1,0.0]; ep = [1, 1]; n = [-1,0.0]
rightw = Wall(sp, ep, n, A)
sp = [1,1]; ep = [0.0, 1]; n = [0.0,-1]
topw = Wall(sp, ep, n, A)
sp = [1,0.25]; ep = [.75,0]; n = [-1/sqrt(2),1/sqrt(2)] 
diagBotRight = Wall(sp, ep, n, A)
sp = [0.0,0.25]; ep = [.25, 0.0]; n = [1/sqrt(2),1/sqrt(2)] 
diagBotLeft = Wall(sp, ep, n, A)
sp = [1,.75]; ep = [.75,1]; n = [-1/sqrt(2),-1/sqrt(2)] 
diagTopRight = Wall(sp, ep, n, A)
sp = [0.0,0.75]; ep = [.25, 1]; n = [1/sqrt(2),-1/sqrt(2)] 
diagtopLeft = Wall(sp, ep, n, A)
sp = [0.0,0.0]; ep = [1, 0.0]; n = [0.0,1]
botw = Wall(sp, ep, n, A)

bd = (botw, rightw, topw, leftw, diagBotRight, diagBotLeft, diagtopLeft, diagTopRight)
bd isa Billiard

p = Particle(0.5, 0.5, 2π*rand())

xt, yt, t = timeseries!(p, bd, 10)


using PyPlot
import PyPlot: plot


const EDGECOLOR = (0,0.6,0)
function plot(w::Wall)
    PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
        color=EDGECOLOR, lw = 2.0)
end
function plot(bd::Billiard)
    for o ∈ bd; plot(o); end
    gca()[:set_aspect]("equal")
end

"
figure(); plot(bd)

figure(); plot(bd)
xt, yt, t = timeseries!(p, bd, 10)
plot(xt, yt)

"

figure(); plot(bd)
p = Particle(0.4, 0.4, π/4)
for j in 1:10
    xt, yt = timeseries!(p, bd, 20)
    plot(xt, yt, alpha = 0.5)
end




