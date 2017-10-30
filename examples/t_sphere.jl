# TODO: work this hypothetical spherical example to see how things fit together

Using LibHealpix

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# Tmap
struct Tmap{P<:Healpix,T<:Real} <: Field{P,T,S0}
    tn::RingHealpixMap{T}
    Tmap{P,T}(tn::Array) where {P<:Flat,T<:Real} = new{P,T}(RingHealpixMap{T}(tn))
end
has_qu(::Type{Tmap{P,T}}) where {P<:Healpix,T<:Real} = HasQU{false}
is_map(::Type{Tmap{P,T}}) where {P<:Healpix,T<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{P,T}}) where {P<:Healpix,T<:Real} = IsLenseBasis{true}


# Tlm
struct Tlm{P<:Pix,T<:Real} <: Field{P,T,S0}
    tlm::Alm{T}
    Tlm{P,T}(tlm::Array) where {P<:Flat,T<:Real} = new{P,T}(Alm{T}(tk))
end
has_qu(::Type{Tlm{P,T}}) where {P<:Healpix,T<:Real} = HasQU{false}
is_map(::Type{Tlm{P,T}}) where {P<:Healpix,T<:Real} = IsMap{false}


const MyField{P,T} = Union{Tlm{P,T}, Tmap{P,T}}


############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:MyField{P,T} where {P<:Flat, T<:Real}
    return rℍ(Px,Tx)
end
