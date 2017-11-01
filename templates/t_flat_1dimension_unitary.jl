#=
    To load use

    ```
    include(joinpath(Pkg.dir("FieldsBase"), "templates", "t_flat_1dimension_unitary.jl"))
    ```
=#



############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# Fmap
struct Fmap{P<:Flat,T<:Real} <: Field{P,T,S0}
    fx::Vector{T}
    Fmap{P,T}(fx::Vector) where {P<:Flat,T<:Real} = new{P,T}(fx)
end
has_qu(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Ffourier
struct Ffourier{P<:Pix,T<:Real} <: Field{P,T,S0}
    fk::Vector{Complex{T}}
    Ffourier{P,T}(fk::Vector) where {P<:Flat,T<:Real} = new{P,T}(complex.(fk))
end
has_qu(::Type{Ffourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Ffourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const FField{P,T} = Union{Ffourier{P,T}, Fmap{P,T}}


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
    FFT    =  (nside^(-dm/2)) * plan_rfft(rand(T,nside)) # unitary normization
    r𝔽1d{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, FFT)
end

(*)(g::r𝔽1d{P,T}, x) where {P<:Pix,T} = g.FFT * x
(\)(g::r𝔽1d{P,T}, x) where {P<:Pix,T} = g.FFT \ x

function harmonic_transform(::Type{F}) where F<:FField{P,T} where {P<:Flat, T<:Real}
    return r𝔽1d(P,T)
end

function white_noise(g::r𝔽1d{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n)  #<----------- needs checking
end

# dot(map, map)
function _dot(a::X, b::X, ::Type{IsMap{true}}) where X<:Field{P} where P<:Flat
    _realdot(a,b) = sum(a.*b)
    sum(map(_realdot, data(a),data(b)))
end


# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:FField{P,T} where {T,P<:Flat{θ,m}} where {θ,m}
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
    sum(map(_complexdot, data(a),data(b))) 
end


############################################################
#  getindex (optional)
############################################################

import Base: getindex

# NOTE: these are not type stable
function getindex(f::FField{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :fk)  ? Ffourier{P,T}(f).fk :
    (sym == :fx)  ? Fmap{P,T}(f).fx :
    error("index is not defined")
end


