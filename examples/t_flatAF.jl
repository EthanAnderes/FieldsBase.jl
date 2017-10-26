# ArrayFire tempurature fields


using ArrayFire

# Basic useage for ArrayFire
# a = AFArray(rand(Complex{Float32}, 10,10))
# a = rand(AFArray{Float32}, 100, 100)

############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# AFTmap
struct AFTmap{P<:Flat,T<:Real} <: Field{P,S0}
    tx::AFArray{T,2}
    AFTmap{P,T}(tx::Array) where {P<:Flat,T<:Real} = new{P,T}(AFArray(T.(tx)))
    AFTmap{P,T}(tx::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tx)
end
has_qu(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct AFTfourier{P<:Pix,T<:Real} <: Field{P,S0}
    tk::AFArray{Complex{T},2}
    AFTfourier{P,T}(tk::Array) where {P<:Flat,T<:Real}  = new{P,T}(AFArray(complex.(T.(tk))))
    AFTfourier{P,T}(tk::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tk)
end
has_qu(::Type{AFTfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{AFTfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const AFS0Field{P,T} = Union{AFTfourier{P,T}, AFTmap{P,T}}


############################################################
#  Specify the harmonic transform
############################################################

import Base: *, \, +, -, ^

#  ArrayFire FFT
struct AFr𝔽{P<:Flat,T<:Real} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωpix::T
    period::T
    nyq::T
    k::Vector{AFArray{T,2}}
    x::Vector{AFArray{T,2}}
    sin2ϕk::AFArray{T,2}
    cos2ϕk::AFArray{T,2}
end

# real FFT generated function constructor
@generated function AFr𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωpix   = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [AFArray(T.(reshape(k_side, 1, nside))), AFArray(T.(reshape(k_side[1:nside÷2+1], nside÷2+1, 1)))]
    x      = [AFArray(T.(reshape(x_side, 1, nside))), AFArray(T.(reshape(x_side, nside, 1)))]
    ϕk     = AFArray(T.( atan2.(Array(k[2]), Array(k[1]) )))
    AFr𝔽{P,T}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk))
end

(*)(g::AFr𝔽{P,T}, x) where {P<:Pix,T} = T(g.Ωpix / (2π)) * rfft(x)
(\)(g::AFr𝔽{P,T}, x) where {P<:Pix,T} = T((2π) / g.Ωpix) * irfft(x)

function harmonic_transform(::Type{F}) where F<:AFS0Field{P,T} where {P<:Flat, T<:Real}
    return AFr𝔽(P,T)
end

#= -------------------------------------
The following additional definitions are not usually needed, but in this case
ArrayFire is too eager to promote AFArray{Float32} to AFArray{Float64}:
  e.g. 1.0 .* AFArray{Float32} = AFArray{Float64}
  e.g. AFArray{Float32}.^(5.0) = AFArray{Float64}
=# #-------------------------------

(+)(f::F, n::Number) where F<:AFS0Field{P,T} where {P,T<:Real} = F((data(f) .+ T(n))...)
(+)(n::Number, f::F) where F<:AFS0Field{P,T} where {P,T<:Real} = F((data(f) .+ T(n))...)
(-)(a::F)            where F<:AFS0Field{P,T} where {P,T<:Real} = F((.- data(a))...)
(-)(f::F, n::Number) where F<:AFS0Field{P,T} where {P,T<:Real} = F((data(f) .- T(n))...)
(-)(n::Number, f::F) where F<:AFS0Field{P,T} where {P,T<:Real} = F((T(n) .- data(f))...)
(*)(f::F, n::Number) where F<:AFS0Field{P,T} where {P,T<:Real} = F((T(n) .* data(f))...)
(*)(n::Number, f::F) where F<:AFS0Field{P,T} where {P,T<:Real} = F((T(n) .* data(f))...)
(^)(op::𝕃{F}, a::Number)  where F<:AFS0Field{P,T} where {P,T<:Real} = 𝕃(F((i.^T(a) for i in data(op.f))...))
(^)(op::𝕃{F}, a::Integer) where F<:AFS0Field{P,T} where {P,T<:Real} = 𝕃(F((i.^T(a) for i in data(op.f))...))


# TODO: figure out how to overload _dot so it works for ArrayFire



############################################################
#  The fields are ready to go ...
############################################################


using Base.Test
# using Test # for v0.7

nside  = 512
Θpix   = 2.0
Px     = Flat{Θpix,nside}
Tx     = Float32
g      =  AFr𝔽(Px,Tx)

#### test AFr𝔽
c = AFArray(rand(Complex{Tx}, nside÷2+1,nside))
r = rand(AFArray{Tx}, nside, nside)
g * r
g \ c


tx = AFArray(rand(Tx, nside, nside))
tk = AFArray(rand(Complex{Tx}, nside÷2+1, nside))

t1 = AFTmap{Px,Tx}(tx)
t2 = AFTfourier{Px,Tx}(tk)


2 * t1 - 5 * t1
2 * t1 - Tx(5.0) * t1
2 * t1 - 5.0 * t1


@inferred 2 * t1 - 5.0 * t1
@inferred 2 * t1 - 5.0 * t2
@inferred 2 * t2 - 5.0 * t1
@inferred 2 * t2 - 5.0 * t2



t1x = rand(Tx, nside, nside)
t2x = rand(Tx, nside, nside)
t1 = AFTmap{Px,Tx}(t1x)
t2 = AFTmap{Px,Tx}(t2x)


#NOTE dot has problems ... since vecdot isn't implimented in ArrayFire yet ....

t1 = AFTmap{Px,Tx}(tx)
t2 = AFTfourier{Px,Tx}(tk)

@inferred 𝕃(t1)*t1
@inferred 𝕃(t1)*t2
@inferred 𝕃(t2)*t1
@inferred 𝕃(t2)*t2


L1 = 𝕃(t1)^(-5.3)
L2 = 𝕃(t2)^(2)
L3 = 𝕃(t1)^(-.1)
L4 = 𝕃(t2)^(4)
L5 = 𝕃(t1)^(1)
L6 = 𝕃(t2)^(0)
L7 = 𝕃(t1)^(-1)
L8 = 𝕃(t2)^(1.5)
L9 = inv(𝕃(t1)) #<---- #NOTE error here because of squash
L10 = inv(𝕃(t2)) #<----#NOTE error here because of squash
L11 = 𝕃(t2)^(-1)
#NOTE squash needs to be defined for array fire..

@inferred L1*t1
@inferred L2*t1
@inferred L3*t2
@inferred L4*t2
@inferred L5*t1
@inferred L6*t1
@inferred L7*t2
@inferred L8*t2
@inferred L9*t2 - L7*t2
@inferred L9*t1 - L7*t1
@inferred L10*t1 - L11*t1
