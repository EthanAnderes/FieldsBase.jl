# To load use:
#
# FieldsBase_dir = dirname(dirname(FieldsBase.source_path))
# include(joinpath(FieldsBase_dir,"templates/tqu_flat.jl"))


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# TQUmap
struct TQUmap{P<:Flat,T<:Real} <: Field{P,T,S02}
    tx::Matrix{T}
    qx::Matrix{T}
    ux::Matrix{T}
    TQUmap{P,T}(tx::Matrix, qx::Matrix, ux::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx, qx, ux)
end
has_qu(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# TQUfourier
struct TQUfourier{P<:Flat,T<:Real} <: Field{P,T,S02}
    tk::Matrix{Complex{T}}
    qk::Matrix{Complex{T}}
    uk::Matrix{Complex{T}}
    TQUfourier{P,T}(tk::Matrix, qk::Matrix, uk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk), complex.(qk), complex.(uk))
end
has_qu(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# TEBmap
struct TEBmap{P<:Pix, T<:Real} <: Field{P,T,S02}
    tx::Matrix{T}
    ex::Matrix{T}
    bx::Matrix{T}
    TEBmap{P,T}(tx::Matrix, ex::Matrix, bx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx, ex, bx)
end
has_qu(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# TEBfourier
struct TEBfourier{P<:Pix, T<:Real} <: Field{P,T,S02}
    tk::Matrix{Complex{T}}
    ek::Matrix{Complex{T}}
    bk::Matrix{Complex{T}}
    TEBfourier{P,T}(tk::Matrix, ek::Matrix, bk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk), complex.(ek), complex.(bk))
end
has_qu(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


const S02Field{P,T} = Union{TEBfourier{P,T}, TEBmap{P,T}, TQUmap{P,T}, TQUfourier{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
TQUmap{P,T}(tp::S02Field{P,T})     where {P<:Flat,T<:Real} = convert(TQUmap{P,T}, tp) 
TQUfourier{P,T}(tp::S02Field{P,T}) where {P<:Flat,T<:Real} = convert(TQUfourier{P,T}, tp) 
TEBmap{P,T}(tp::S02Field{P,T})     where {P<:Flat,T<:Real} = convert(TEBmap{P,T}, tp) 
TEBfourier{P,T}(tp::S02Field{P,T}) where {P<:Flat,T<:Real} = convert(TEBfourier{P,T}, tp) 




############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:S02Field{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end






############################################################
#  getindex (optional)
############################################################

import Base: getindex

# NOTE: these are not type stable
function getindex(f::S02Field{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :ek)  ? TEBfourier{P,T}(f).ek :
    (sym == :bk)  ? TEBfourier{P,T}(f).bk :
    (sym == :tk)  ? TEBfourier{P,T}(f).tk :
    (sym == :qk)  ? TQUfourier{P,T}(f).qk :
    (sym == :uk)  ? TQUfourier{P,T}(f).uk :
    (sym == :ex)  ? TEBmap{P,T}(f).ex :
    (sym == :bx)  ? TEBmap{P,T}(f).bx :
    (sym == :qx)  ? TQUmap{P,T}(f).qx :
    (sym == :ux)  ? TQUmap{P,T}(f).ux :
    (sym == :tx)  ? TQUmap{P,T}(f).tx :
    error("index is not defined")
end

