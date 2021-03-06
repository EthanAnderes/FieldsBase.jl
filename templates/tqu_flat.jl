# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/tqu_flat.jl"))


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
using FieldsBase: r𝔽
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# TQUmap
struct TQUmap{P<:Flat,T<:Real} <: Field{P,T,S02}
    tx::Matrix{T}
    qx::Matrix{T}
    ux::Matrix{T}
    TQUmap{P,T}(tqux) where {P<:Flat,T<:Real} = new{P,T}(T.(tqux), T.(tqux), T.(tqux))
    TQUmap{P,T}(tx, qx, ux) where {P<:Flat,T<:Real} = new{P,T}(T.(tx), T.(qx), T.(ux))
    TQUmap{P,T}(tx::Matrix{T}, qx::Matrix{T}, ux::Matrix{T}) where {P<:Flat,T<:Real} = new{P,T}(tx, qx, ux)
    TQUmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T))
end
has_qu(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# TQUfourier
struct TQUfourier{P<:Flat,T<:Real} <: Field{P,T,S02}
    tk::Matrix{Complex{T}}
    qk::Matrix{Complex{T}}
    uk::Matrix{Complex{T}}
    TQUfourier{P,T}(tquk) where {P<:Flat,T<:Real} = new{P,T}(Complex{T}.(tquk), Complex{T}.(tquk), Complex{T}.(tquk))
    TQUfourier{P,T}(tk, qk, uk) where {P<:Flat,T<:Real} = new{P,T}(Complex{T}.(tk), Complex{T}.(qk), Complex{T}.(uk))
    TQUfourier{P,T}(tk::Matrix{CT}, qk::Matrix{CT}, uk::Matrix{CT}) where {P<:Flat,T<:Real,CT<:Complex{T}} = new{P,T}(tk, qk, uk)
    TQUfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T))
end
has_qu(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# TEBmap
struct TEBmap{P<:Flat, T<:Real} <: Field{P,T,S02}
    tx::Matrix{T}
    ex::Matrix{T}
    bx::Matrix{T}
    TEBmap{P,T}(tebx) where {P<:Flat,T<:Real} = new{P,T}(T.(tebx), T.(tebx), T.(tebx))
    TEBmap{P,T}(tx, ex, bx) where {P<:Flat,T<:Real} = new{P,T}(T.(tx), T.(ex), T.(bx))
    TEBmap{P,T}(tx::Matrix{T}, ex::Matrix{T}, bx::Matrix{T}) where {P<:Flat,T<:Real} = new{P,T}(tx, ex, bx)
    TEBmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T), FieldsBase.zero_map_dim2(P,T))
end
has_qu(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# TEBfourier
struct TEBfourier{P<:Flat, T<:Real} <: Field{P,T,S02}
    tk::Matrix{Complex{T}}
    ek::Matrix{Complex{T}}
    bk::Matrix{Complex{T}}
    TEBfourier{P,T}(tebk) where {P<:Flat,T<:Real} = new{P,T}(Complex{T}.(tebk), Complex{T}.(tebk), Complex{T}.(tebk))
    TEBfourier{P,T}(tk, ek, bk) where {P<:Flat,T<:Real} = new{P,T}(Complex{T}.(tk), Complex{T}.(ek), Complex{T}.(bk))
    TEBfourier{P,T}(tk::Matrix{CT}, ek::Matrix{CT}, bk::Matrix{CT}) where {P<:Flat,T<:Real,CT<:Complex{T}} = new{P,T}(tk, ek, bk)
    TEBfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T), FieldsBase.zero_fourier_dim2(P,T))
end
has_qu(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# const S02Field{P,T} = Union{TEBfourier{P,T}, TEBmap{P,T}, TQUmap{P,T}, TQUfourier{P,T}}


# This is needed for 0.7 since constructors do not fall back on convert
TQUmap{P,T}(tp::Field{P,T,S02})     where {P<:Flat,T<:Real} = convert(TQUmap{P,T}, tp) 
TQUfourier{P,T}(tp::Field{P,T,S02}) where {P<:Flat,T<:Real} = convert(TQUfourier{P,T}, tp) 
TEBmap{P,T}(tp::Field{P,T,S02})     where {P<:Flat,T<:Real} = convert(TEBmap{P,T}, tp) 
TEBfourier{P,T}(tp::Field{P,T,S02}) where {P<:Flat,T<:Real} = convert(TEBfourier{P,T}, tp) 




############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:Field{P,T,S02} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end






############################################################
#  getindex (optional)
############################################################

import Base: getindex

# NOTE: these are not type stable
function getindex(f::Field{P,T,S02}, sym::Symbol) where {P<:Flat, T<:Real}
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

