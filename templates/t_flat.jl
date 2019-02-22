# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/t_flat.jl"))


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform


# Tmap
struct Tmap{P<:Flat,T<:Real} <: Field{P,T,S0}
    tx::Matrix{T}
    Tmap{P,T}(tx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx)
    Tmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim2(P,T))

end
has_qu(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct Tfourier{P<:Pix,T<:Real} <: Field{P,T,S0}
    tk::Matrix{Complex{T}}
    Tfourier{P,T}(tk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk))
    Tfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim2(P,T))
end
has_qu(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


const S0Field{P,T} = Union{Tfourier{P,T}, Tmap{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
Tfourier{P,T}(t::S0Field{P,T}) where {P<:Flat,T<:Real} = convert(Tfourier{P,T}, t) 
Tmap{P,T}(t::S0Field{P,T}) where {P<:Flat,T<:Real}     = convert(Tmap{P,T}, t) 


############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:S0Field{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end


############################################################
#  getindex (optional)
############################################################

import Base: getindex

function getindex(f::S0Field{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :tk)  ? Tfourier{P,T}(f).tk :
    (sym == :tx)  ? Tmap{P,T}(f).tx :
    error("index is not defined")
end

