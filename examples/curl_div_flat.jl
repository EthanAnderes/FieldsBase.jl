#=
In this example (Q,U) represents a regular vector field (not headless)
and E = div(Q,U), B = curl(Q,U)
=#


if VERSION >= v"0.7.0-DEV.1"
    using FFTW
end
FFTW.set_num_threads(6)
BLAS.set_num_threads(6)

############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# the Vmap field is designated as the QU field
struct Vmap{P<:Flat,T<:Real} <: Field{P,S2}
    v1x::Matrix{T}
    v2x::Matrix{T}
    Vmap{P,T}(v1x::Matrix, v2x::Matrix) where {P<:Flat,T<:Real} = new{P,T}(v1x, v2x)
end
has_qu(::Type{Vmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{Vmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Vmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


struct Vfourier{P<:Pix,T<:Real} <: Field{P,S2}
    v1k::Matrix{Complex{T}}
    v2k::Matrix{Complex{T}}
    Vfourier{P,T}(v1k::Matrix, v2k::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(v1k), complex.(v2k))
end
has_qu(::Type{Vfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{Vfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

# H for Helmholtz
struct Hmap{P<:Pix, T<:Real} <: Field{P,S2}
    dx::Matrix{T} #<-- div
    cx::Matrix{T} #<-- curl
    Hmap{P,T}(dx::Matrix, cx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(dx, cx)
end
has_qu(::Type{Hmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Hmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


struct Hfourier{P<:Pix, T<:Real} <: Field{P,S2}
    dk::Matrix{Complex{T}}
    ck::Matrix{Complex{T}}
    Hfourier{P,T}(dk::Matrix, ck::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(dk), complex.(ck))
end
has_qu(::Type{Hfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{Hfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const MyField{P,T} = Union{Vfourier{P,T}, Vmap{P,T}, Hfourier{P,T}, Hmap{P,T}}

############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:MyField{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end

import FieldsBase: ebk_to_quk, quk_to_ebk
function quk_to_ebk(v1k, v2k, g::r𝔽{P,T}) where {P<:Pix, T<:Real}
    dk = v1k .* (im.*g.k[1]) .+ v2k .* (im.*g.k[2])
    ck = v1k .* (im.*g.k[2]) .- v2k .* (im.*g.k[1])
    return dk, ck
end
function ebk_to_quk(dk, bk, g::r𝔽{P,T}) where {P<:Pix, T<:Real}
    invΔk = squash.( .- (g.k[1].^2 .+ g.k[2].^2) )
    v1k = invΔk .* ( (im.*g.k[1]) .* dk .+ (im.*g.k[2]) .* ck )
    v2k = invΔk .* ( (im.*g.k[2]) .* dk .- (im.*g.k[1]) .* ck )
    return v1k, v2k
end
# NOTE: can we get rid of the "Method definition overwritten"?

# (im*g.k[1]) * dk + (im*g.k[2]) * ck
# = - v1k * (g.k[1]^2)  - v2k * (g.k[2].*g.k[1])
#   - v1k * (g.k[2]^2)  + v2k * (g.k[2]*g.k[1])
# = - 2 v1k .* (g.k[1]^2 + g.k[2]^2)
#
# (im*g.k[2]) * dk - (im*g.k[1]) * ck
# =   - v1k * (g.k[2].*g.k[1])  - v2k * (g.k[2]^2)
#   -[- v1k * (g.k[2].*g.k[1])  + v2k * (g.k[1]^2) ]
# = - 2 v2k .* (g.k[1]^2 + g.k[2]^2)







############################################################
#  The fields are ready to go ...
############################################################

nside  = 512
Θpix   = 2.0
P     = Flat{Θpix,nside}
T     = Float32
g     = r𝔽(P,T);

v1x, v2x = rand(T, nside, nside), rand(T, nside, nside)
dx, cx = rand(T, nside, nside), rand(T, nside, nside)
v1k, v2k = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)
dk, ck = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)

p1 = Vmap{P,T}(v1x, v2x)
p2 = Hmap{P,T}(dx, cx)
p3 = Vfourier{P,T}(v1k, v2k)
p4 = Hfourier{P,T}(dk, ck)

p = convert(Vfourier{P,T}, p1)
p.v1k - g * p1.v1x

@time convert(Vfourier{P,T}, p4)
@time convert(Vfourier{P,T}, p1)



p = convert(Hfourier{P,T}, p1)
p = convert(Hfourier{P,T}, p3)
p = convert(Vfourier{P,T}, p3)



using Base.Test

@inferred Vfourier{P,T}(p1)
@inferred Vfourier{P,T}(p2)
@inferred Vfourier{P,T}(p3)
@inferred Vfourier{P,T}(p4)

@inferred 2 * p1 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2
@inferred 2 * p1 - 5.0 * p3
@inferred 2 * p1 - 5.0 * p4
@inferred 2 * p2 - 5.0 * p1
@inferred 2 * p2 - 5.0 * p2
@inferred 2 * p2 - 5.0 * p3
@inferred 2 * p2 - 5.0 * p4
@inferred 2 * p3 - 5.0 * p1
@inferred 2 * p3 - 5.0 * p2
@inferred 2 * p3 - 5.0 * p3
@inferred 2 * p3 - 5.0 * p4
@inferred 2 * p4 - 5.0 * p1
@inferred 2 * p4 - 5.0 * p2
@inferred 2 * p4 - 5.0 * p3
@inferred 2 * p4 - 5.0 * p4



##### Testing dot
r1, r2, r3, r4 = randn(T,nside, nside), randn(T,nside, nside), randn(T,nside, nside), randn(T,nside, nside)
p1, p2 = Vmap{P,T}(r1, r2), Vmap{P,T}(r3, r4)

@test dot(p1, p2) == dot(p2, p1)
@test dot(p1, p2) == (dot(r1,r3) + dot(r2,r4))*r𝔽(P,T).Ωpix
@test dot(p1, p1) > 0

@inferred dot(p1, p2)
@inferred dot(Hfourier{P,T}(p1), Hfourier{P,T}(p2))
@inferred dot(Vfourier{P,T}(p1), p2)
@inferred dot(Vfourier{P,T}(p1), Vfourier{P,T}(p2))
@inferred dot(Hfourier{P,T}(p1), p2)
@inferred dot(p1, Hfourier{P,T}(p2))
@inferred dot(Vmap{P,T}(p1), p2)
@inferred dot(p1, Vmap{P,T}(p2))
@inferred dot(Vmap{P,T}(p1), Vmap{P,T}(p2))
@inferred dot(p1, Vmap{P,T}(Vfourier{P,T}(p2)))

@inferred dot(p1, Hmap{P,T}(p2))
@inferred dot(Hmap{P,T}(p1), Hmap{P,T}(p2))
@inferred dot(Hmap{P,T}(p1), p2)


wn1 = white_noise(g)
wn2 = white_noise(g)
p = Vmap{P,T}(wn1, wn2)
dot(p, p)/2/nside^2 # this should be near 2nside^2


##### Testing DiagOp
p1 = Vmap{P,T}(v1x, v2x)
p2 = Hmap{P,T}(dx, cx)
p3 = Vfourier{P,T}(v1k, v2k)
p4 = Hfourier{P,T}(dk, ck)

@inferred 𝕃(p1)*p1
@inferred 𝕃(p1)*p2
@inferred 𝕃(p1)*p3
@inferred 𝕃(p1)*p4

@inferred 𝕃(p2)*p1
@inferred 𝕃(p2)*p2
@inferred 𝕃(p2)*p3
@inferred 𝕃(p2)*p4

@inferred 𝕃(p3)*p1
@inferred 𝕃(p3)*p2
@inferred 𝕃(p3)*p3
@inferred 𝕃(p3)*p4

@inferred 𝕃(p4)*p1
@inferred 𝕃(p4)*p2
@inferred 𝕃(p4)*p3
@inferred 𝕃(p4)*p4

L1 = 𝕃(p1)^(-5.3)
L2 = 𝕃(p2)^(2)
L3 = 𝕃(p3)^(-.1)
L4 = 𝕃(p4)^(4)
L5 = 𝕃(p1)^(1)
L6 = 𝕃(p2)^(0)
L7 = 𝕃(p1)^(-1)
L8 = 𝕃(p4)^(1.5)
L9 = inv(𝕃(p1))
L10 = inv(𝕃(p2))
L11 = 𝕃(p2)^(-1)

L1*p4
L2*p3
L3*p2
L4*p1
L5*p2
L6*p2
L7*p2
L8*p2
L9*p2 - L7*p2
L9*p3 - L7*p3
L10*p4 - L11*p4
