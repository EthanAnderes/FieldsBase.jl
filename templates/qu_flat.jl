# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/qu_flat.jl"))


############################################################
#  Define the field types and their trait properties
############################################################


using FieldsBase
using FieldsBase: r𝔽
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# QUmap
struct QUmap{P<:Flat,T<:Real} <: Field{P,T,S2}
    qx::Matrix{T}
    ux::Matrix{T}
    QUmap{P,T}(qx::Matrix, ux::Matrix) where {P<:Flat,T<:Real} = new{P,T}(qx, ux)
    QUmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T))
end
has_qu(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# QUfourier
struct QUfourier{P<:Flat,T<:Real} <: Field{P,T,S2}
    qk::Matrix{Complex{T}}
    uk::Matrix{Complex{T}}
    QUfourier{P,T}(qk::Matrix, uk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(qk), complex.(uk))
    QUfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T))
end
has_qu(::Type{QUfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# EBmap
struct EBmap{P<:Flat, T<:Real} <: Field{P,T,S2}
    ex::Matrix{T}
    bx::Matrix{T}
    EBmap{P,T}(ex::Matrix, bx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(ex, bx)
    EBmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T))
end
has_qu(::Type{EBmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{EBmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# EBfourier
struct EBfourier{P<:Flat, T<:Real} <: Field{P,T,S2}
    ek::Matrix{Complex{T}}
    bk::Matrix{Complex{T}}
    EBfourier{P,T}(ek::Matrix, bk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(ek), complex.(bk))
    EBfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T))
end
has_qu(::Type{EBfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{EBfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# const S2Field{P,T} = Union{EBfourier{P,T}, EBmap{P,T}, QUfourier{P,T}, QUmap{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
QUmap{P,T}(p::Field{P,T,S2})     where {P<:Flat,T<:Real} = convert(QUmap{P,T}, p)
QUfourier{P,T}(p::Field{P,T,S2}) where {P<:Flat,T<:Real} = convert(QUfourier{P,T}, p)
EBmap{P,T}(p::Field{P,T,S2})     where {P<:Flat,T<:Real} = convert(EBmap{P,T}, p)
EBfourier{P,T}(p::Field{P,T,S2}) where {P<:Flat,T<:Real} = convert(EBfourier{P,T}, p)




############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:Field{P,T,S2} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end



############################################################
#  getindex (optional)
############################################################

import Base: getindex

function getindex(f::Field{P,T,S2}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :ek)  ? EBfourier{P,T}(f).ek :
    (sym == :bk)  ? EBfourier{P,T}(f).bk :
    (sym == :qk)  ? QUfourier{P,T}(f).qk :
    (sym == :uk)  ? QUfourier{P,T}(f).uk :
    (sym == :ex)  ? EBmap{P,T}(f).ex :
    (sym == :bx)  ? EBmap{P,T}(f).bx :
    (sym == :qx)  ? QUmap{P,T}(f).qx :
    (sym == :ux)  ? QUmap{P,T}(f).ux :
    error("index is not defined")
end


