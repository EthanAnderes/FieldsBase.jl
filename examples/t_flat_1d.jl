# test to see if we can use FieldsBase to define 1-d stochastic process.

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

# Tmap
struct Tmap{P<:Flat,T<:Real} <: Field{P,S0}
    tx::Vector{T}
    Tmap{P,T}(tx::Vector) where {P<:Flat,T<:Real} = new{P,T}(tx)
end
has_qu(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct Tfourier{P<:Pix,T<:Real} <: Field{P,S0}
    tk::Vector{Complex{T}}
    Tfourier{P,T}(tk::Vector) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk))
end
has_qu(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const MyField{P,T} = Union{Tfourier{P,T}, Tmap{P,T}}


############################################################
#  Specify the harmonic transform
############################################################

import Base: *, \
import FieldsBase: _dot, white_noise

#  1-d real FFT
struct r𝔽1d{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωpix::T
    period::T
    nyq::T
    k::Vector{T}
    x::Vector{T}
    FFT::F
end

@generated function r𝔽1d(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk
    Ωpix   = Δx
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = k_side[1:nside÷2+1]
    x      = x_side
    dm     = 1 #<-- dimension
    FFT    =  ( (Δx / √(2π))^dm ) * plan_rfft(rand(T,nside))
    r𝔽1d{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, FFT)
end

(*)(g::r𝔽1d{P,T}, x) where {P<:Pix,T} = g.FFT * x
(\)(g::r𝔽1d{P,T}, x) where {P<:Pix,T} = g.FFT \ x

function harmonic_transform(::Type{F}) where F<:MyField{P,T} where {P<:Flat, T<:Real}
    return r𝔽1d(P,T)
end

function white_noise(g::r𝔽1d{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n) ./ sqrt(g.Ωpix)
end

# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:MyField{P,T} where {T,P<:Flat{θ,m}} where {θ,m}
    FT = harmonic_transform(X)
    _complexdot = function (a,b)
        ra, ia = real(a), imag(a)
        rb, ib = real(b), imag(b)
        rtn  = 2*sum(ra.*rb)
        rtn += 2*sum(ia.*ib)
        rtn -= sum(ra[1].*rb[1])
        rtn -= sum(ia[1].*ib[1])
        if iseven(m)
               rtn -= sum(ra[end].*rb[end])
               rtn -= sum(ia[end].*ib[end])
        end
        return rtn
    end
    sum(map(_complexdot, data(a),data(b))) * FT.Ωk
end



############################################################
#  The fields are ready to go ...
############################################################

using Base.Test
nside  = 2^12
Θpix   = 2.0
P     = Flat{Θpix,nside}
T     = Float32
g      =  r𝔽1d(P,T);
tx = rand(T, nside)
tk = rand(Complex{T}, nside÷2+1)



#### test r𝔽1d
c = rand(T, nside÷2+1)
r = rand(T, nside)
g \ (g * r) .- r
g * (g \ c) .- c


t1 = Tmap{P,T}(tx)
t2 = Tfourier{P,T}(tk)

t = convert(Tfourier{P,T}, t1)
@test all(t.tk - g * t1.tx .== 0)

@test typeof(convert(Tfourier{P,T}, t1)) == Tfourier{P,T}
@test typeof(convert(Tfourier{P,T}, t2)) == Tfourier{P,T}
@test typeof(convert(Tmap{P,T}, t1)) == Tmap{P,T}
@test typeof(convert(Tmap{P,T}, t2)) == Tmap{P,T}

@inferred convert(Tfourier{P,T}, t1)
@inferred convert(Tfourier{P,T}, t2)
@inferred convert(Tmap{P,T}, t1)
@inferred convert(Tmap{P,T}, t2)

@inferred 2 * t1 - 5.0 * t1
@inferred 2 * t1 - 5.0 * t2
@inferred 2 * t2 - 5.0 * t1
@inferred 2 * t2 - 5.0 * t2



##### Testing dot
t1x = rand(T, nside)
t2x = rand(T, nside)
t1 = Tmap{P,T}(t1x)
t2 = Tmap{P,T}(t2x)

@test dot(t1, t2) == dot(t2, t1)
@test dot(t1, t2) == dot(t1x,t2x)*r𝔽1d(P,T).Ωpix
@test dot(t1, t1) > 0

@inferred dot(t1, t2)
@inferred dot(Tfourier{P,T}(t1), Tfourier{P,T}(t2))
@inferred dot(Tfourier{P,T}(t1), t2)
@inferred dot(t1, Tfourier{P,T}(t2))
@inferred dot(Tmap{P,T}(t1), t2)
@inferred dot(t1, Tmap{P,T}(t2))
@inferred dot(Tmap{P,T}(t1), Tmap{P,T}(t2))


wn1 = white_noise(g)
t = Tmap{P,T}(wn1)
dot(t, t)/nside # this should be near 1.0
dot(Tfourier{P,T}(t), Tfourier{P,T}(t))/nside # this should be near nside^2


##### Testing DiagOp
t1 = Tmap{P,T}(tx)
t2 = Tfourier{P,T}(tk)

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
