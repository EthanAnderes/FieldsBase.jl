# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/s_flat_1dimension_unitary.jl"))



############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
using FieldsBase: r𝕌𝔽1
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform


# Smap
struct Smap{P<:Flat,T<:Real} <: Field{P,T,S0}
    sx::Vector{T}
    Smap{P,T}(sx::Vector) where {P<:Flat,T<:Real} = new{P,T}(sx)
    Smap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim1(P,T))
end
has_qu(::Type{Smap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Smap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Smap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Sfourier
struct Sfourier{P<:Flat,T<:Real} <: Field{P,T,S0}
    sk::Vector{Complex{T}}
    Sfourier{P,T}(sk::Vector) where {P<:Flat,T<:Real} = new{P,T}(complex.(sk))
    Sfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim1(P,T))
end
has_qu(::Type{Sfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Sfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


const SField{P,T} = Union{Sfourier{P,T}, Smap{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
Smap{P,T}(s::SField{P,T})     where {P<:Flat,T<:Real} = convert(Smap{P,T}, s) 
Sfourier{P,T}(s::SField{P,T}) where {P<:Flat,T<:Real} = convert(Sfourier{P,T}, s) 


############################################################
#  Specify the harmonic transform
############################################################

import FieldsBase: _dot, white_noise

function harmonic_transform(::Type{F}) where F<:SField{P,T} where {P<:Flat, T<:Real}
    return r𝕌𝔽1(P,T)
end

function white_noise(g::r𝕌𝔽1{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n)  #<----------- needs checking
end

# dot(map, map)
function _dot(a::X, b::X, ::Type{IsMap{true}}) where X<:SField{P} where P<:Flat
    _realdot(a,b) = sum(a.*b)
    sum(map(_realdot, data(a),data(b)))
end


# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:SField{P,T} where {T,P<:Flat{θ,m}} where {θ,m}
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
function getindex(f::SField{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :k)  ? Sfourier{P,T}(f).sk :
    (sym == :x)  ? Smap{P,T}(f).sx :
    error("index is not defined")
end


