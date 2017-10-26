#=
	To load use

	```
	include(joinpath(Pkg.dir("FieldsBase"), "templates", "tqu_flat.jl"))
	```
=#


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# TQUmap
struct TQUmap{P<:Flat,T<:Real} <: Field{P,S02}
    tx::Matrix{T}
    qx::Matrix{T}
    ux::Matrix{T}
    TQUmap{P,T}(tx::Matrix, qx::Matrix, ux::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx, qx, ux)
end
has_qu(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{TQUmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# TQUfourier
struct TQUfourier{P<:Flat,T<:Real} <: Field{P,S02}
    tk::Matrix{Complex{T}}
    qk::Matrix{Complex{T}}
    uk::Matrix{Complex{T}}
    TQUfourier{P,T}(tk::Matrix, qk::Matrix, uk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk), complex.(qk), complex.(uk))
end
has_qu(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{TQUfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# TEBmap
struct TEBmap{P<:Pix, T<:Real} <: Field{P,S02}
    tx::Matrix{Complex{T}}
    ex::Matrix{Complex{T}}
    bx::Matrix{Complex{T}}
    TEBmap{P,T}(tx::Matrix, ex::Matrix, bx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx, ex, bx)
end
has_qu(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{TEBmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# TEBfourier
struct TEBfourier{P<:Pix, T<:Real} <: Field{P,S02}
    tk::Matrix{Complex{T}}
    ek::Matrix{Complex{T}}
    bk::Matrix{Complex{T}}
    TEBfourier{P,T}(tk::Matrix, ek::Matrix, bk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk), complex.(ek), complex.(bk))
end
has_qu(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{TEBfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


############################################################
#  Specify the harmonic transform
############################################################

const S02Field{P,T} = Union{TEBfourier{P,T}, TQUmap{P,T}}
function harmonic_transform(::Type{F}) where F<:S02Field{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end
