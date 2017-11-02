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
using FieldsBase: r𝕌𝔽1
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

import FieldsBase: _dot, white_noise

function harmonic_transform(::Type{F}) where F<:FField{P,T} where {P<:Flat, T<:Real}
    return r𝕌𝔽1(P,T)
end

function white_noise(g::r𝕌𝔽1{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n)  #<----------- needs checking
end

# dot(map, map)
function _dot(a::X, b::X, ::Type{IsMap{true}}) where X<:FField{P} where P<:Flat
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


