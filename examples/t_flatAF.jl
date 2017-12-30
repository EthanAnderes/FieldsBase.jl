# ArrayFire tempurature fields


using ArrayFire



############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# AFTmap
struct AFTmap{P<:Flat,T<:Real} <: Field{P,T,S0}
    tx::AFArray{T,2}
    AFTmap{P,T}(tx::Array) where {P<:Flat,T<:Real} = new{P,T}(AFArray(T.(tx)))
    AFTmap{P,T}(tx::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tx)
end
has_qu(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{AFTmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct AFTfourier{P<:Pix,T<:Real} <: Field{P,T,S0}
    tk::AFArray{Complex{T},2}
    AFTfourier{P,T}(tk::Array) where {P<:Flat,T<:Real}  = new{P,T}(AFArray(complex.(T.(tk))))
    AFTfourier{P,T}(tk::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tk)
end
has_qu(::Type{AFTfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{AFTfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const AFS0Field{P,T} = Union{AFTfourier{P,T}, AFTmap{P,T}}


#  Specify the harmonic transform
function harmonic_transform(::Type{F}) where F<:AFS0Field{P,T} where {P<:Flat, T<:Real}
    return AFr𝔽(P,T)
end



###########################
# need to define specialized ArrayFire FFT
###########################
#  Specify the harmonic transform
#  ArrayFire FFT
struct AFr𝔽{P<:Flat,T<:Real} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωx::T
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
    Ωx     = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [AFArray(T.(reshape(k_side, 1, nside))), AFArray(T.(reshape(k_side[1:nside÷2+1], nside÷2+1, 1)))]
    x      = [AFArray(T.(reshape(x_side, 1, nside))), AFArray(T.(reshape(x_side, nside, 1)))]
    ϕk     = AFArray(T.( atan2.(Array(k[2]), Array(k[1]) )))
    AFr𝔽{P,T}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk))
end

import Base: *, \
(*)(g::AFr𝔽{P,T}, x) where {P<:Flat,T} = rfft2(x, T(g.Ωx/2/π))
(\)(g::AFr𝔽{P,T}, x) where {P<:Flat{θ,n},T} where {θ,n} = irfft2(x, T(2*π/g.Ωx/n^2))



############################################################
#  The fields are ready to go ...
############################################################


using Base.Test
# using Test # for v0.7

nside  = 1024
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
@time sync((2 * t1 - 5.0 * t2).tk)
function foo(t1, t2)
    a = 2 * t1 - 5 * t2
    b = 2 * t1 - 5 * a
    c = a-b
    sync(c.tk)
    return c
end
@time foo(t1, t2)

#=
# ArrayFire fields (Float32)
0.009777 seconds (114 allocations: 2.813 KiB)
0.009869 seconds (114 allocations: 2.813 KiB)

# Regular arrays (Float32)
0.091618 seconds (69 allocations: 76.135 MiB, 77.88% gc time)
0.025750 seconds (69 allocations: 76.135 MiB, 26.73% gc time)
=#


@inferred 2 * t1 - 5.0 * t1
@inferred 2 * t1 - 5.0 * t2
@inferred 2 * t2 - 5.0 * t1
@inferred 2 * t2 - 5.0 * t2



t1x = rand(Tx, nside, nside)
t2x = rand(Tx, nside, nside)
t1 = AFTmap{Px,Tx}(t1x)
t2 = AFTmap{Px,Tx}(t2x)



##### Testing dot
t1x = rand(Tx, nside, nside)
t2x = rand(Tx, nside, nside)
t1 = AFTmap{Px,Tx}(t1x)
t2 = AFTmap{Px,Tx}(t2x)

@test dot(t1, t2) == dot(t2, t1)
@test dot(t1, t2) == sum(t1x.*t2x)*r𝔽(Px,Tx).Ωx
@test dot(t1, t1) > 0

@inferred dot(t1, t2)
@inferred dot(AFTfourier{Px,Tx}(t1), AFTfourier{Px,Tx}(t2))
@inferred dot(AFTfourier{Px,Tx}(t1), t2)
@inferred dot(t1, AFTfourier{Px,Tx}(t2))
@inferred dot(AFTmap{Px,Tx}(t1), t2)
@inferred dot(t1, AFTmap{Px,Tx}(t2))
@inferred dot(AFTmap{Px,Tx}(t1), AFTmap{Px,Tx}(t2))


wn1 = white_noise(g)
t = AFTmap{Px,Tx}(wn1)
dot(t, t)/nside^2 # this should be near 1.0
dot(AFTfourier{Px,Tx}(t), AFTfourier{Px,Tx}(t))/nside^2 # this should be near nside^2



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
L9 = inv(𝕃(t1)) 
L10 = inv(𝕃(t2))
L11 = 𝕃(t2)^(-1)


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





#=  Basic useage for ArrayFire ... how to get squash to work for ArrayFire?
## I think we need bool ? rtn1 : rtn2 to work for ArrayFire ...
a = rand(Float32, 100, 100)
b = rand(Complex{Float32}, 100, 100)
a[1] = Inf
a[2] = -Inf
a[3] = NaN

b[1] = Inf + im * 0
b[2] = 0 + im * Inf
b[3] = NaN + im * Inf


squash(x::T) where T = ifelse(isfinite(x), x, 0)
squash1(x) = ifelse(isfinite(x), x, eltype(x)(0))
squash2(x::T) where T = ifelse(isfinite(x), x, T(0))
squash3(x::T) where T = ifelse(isfinite(x), x, zero(T))

squash1.(a)
squash1.(b)
ifelse.(isfinite.(b), b, 0)

foo(x) = squash.(x)
foo1(x) = squash1.(x)
foo2(x) = squash2.(x)
foo3(x) = squash3.(x)

@code_warntype foo(a)
@code_warntype foo1(a)
@code_warntype foo2(a)
@code_warntype foo3(a)

using BenchmarkTools
@benchmark foo(b)
@benchmark foo1(b)
@benchmark foo2(b)
@benchmark foo3(b)


a   = rand(Float32, 1_000, 1_000)
b   = rand(Float32, 1_000, 1_000)
c   = rand(Float32, 1_000, 1_000)
aAF = AFArray(a)
bAF = AFArray(b)
cAF = AFArray(c)


ak   = rfft(a)  
akAF = rfft2(aAF)
ak - Array(akAF)

iak   = irfft(ak, 1_000)
iak - a

iakAF = irfft2(akAF, 1/length(a))
iakAF - aAF

plan_rfft(a)
A_mult_B!()


using ArrayFire 
a = AFArray(rand(Float32, 1_000, 1_000))

test1 = rfft2(a, Float32(1/2/π))
a1    = irfft2(test1, Float32(2*π)/length(a)) 
a1 - a



test1 = rfft(a)
a1 = irfft(test1) 
a1 - a


test1 = fft(a)
a1 = ifft(test1)
a1-a



=#