# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/f_flat_2dimension_unitary.jl"))


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
using FieldsBase: r𝕌𝔽2
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform


# Fmap
struct Fmap{P<:Flat,T<:Real} <: Field{P,T,S0}
    fx::Matrix{T}
    Fmap{P,T}(fx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(fx)
end
has_qu(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Fmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Ffourier
struct Ffourier{P<:Pix,T<:Real} <: Field{P,T,S0}
    fk::Matrix{Complex{T}}
    Ffourier{P,T}(fk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(fk))
end
has_qu(::Type{Ffourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Ffourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


const FField{P,T} = Union{Ffourier{P,T}, Fmap{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
Fmap{P,T}(f::FField{P,T})     where {P<:Flat,T<:Real} = convert(Fmap{P,T}, f) 
Ffourier{P,T}(f::FField{P,T}) where {P<:Flat,T<:Real} = convert(Ffourier{P,T}, f) 



############################################################
#  Specify the harmonic transform
############################################################

import FieldsBase: _dot, white_noise

function harmonic_transform(::Type{F}) where F<:FField{P,T} where {P<:Flat, T<:Real}
    return r𝕌𝔽2(P,T)
end

function white_noise(g::r𝕌𝔽2{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n,n)  #<----------- needs checking
end

# dot(map, map)
function _dot(a::X, b::X, ::Type{IsMap{true}}) where X<:FField{P} where P<:Flat
    _realdot(a,b) = sum(a.*b)
    return sum(map(_realdot, data(a),data(b))) 
end

# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:FField{P} where P<:Flat
    sum(map(_complexdot, data(a),data(b))) 
end




############################################################
#  getindex (optional)
############################################################

import Base: getindex

# NOTE: these are not type stable
function getindex(f::FField{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :k)  ? Ffourier{P,T}(f).fk :
    (sym == :x)  ? Fmap{P,T}(f).fx :
    error("index is not defined")
end


