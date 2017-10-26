#=
	To load use

	```
	include(joinpath(Pkg.dir("FieldsBase"), "templates", "qu_flat.jl"))
	```
=#


############################################################
#  Define the field types and their trait properties
############################################################


using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# QUmap
struct QUmap{P<:Flat,T<:Real} <: Field{P,S2}
    qx::Matrix{T}
    ux::Matrix{T}
    QUmap{P,T}(qx::Matrix, ux::Matrix) where {P<:Flat,T<:Real} = new{P,T}(qx, ux)
end
has_qu(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{QUmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# QUfourier
struct QUfourier{P<:Pix,T<:Real} <: Field{P,S2}
    qk::Matrix{Complex{T}}
    uk::Matrix{Complex{T}}
    QUfourier{P,T}(qk::Matrix, uk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(qk), complex.(uk))
end
has_qu(::Type{QUfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# EBmap
struct EBmap{P<:Pix, T<:Real} <: Field{P,S2}
    ex::Matrix{T}
    bx::Matrix{T}
    EBmap{P,T}(ex::Matrix, bx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(ex, bx)
end
has_qu(::Type{EBmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{EBmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# EBfourier
struct EBfourier{P<:Pix, T<:Real} <: Field{P,S2}
    ek::Matrix{Complex{T}}
    bk::Matrix{Complex{T}}
    EBfourier{P,T}(ek::Matrix, bk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(ek), complex.(bk))
end
has_qu(::Type{EBfourier{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{EBfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const S2Field{P,T} = Union{EBfourier{P,T}, EBmap{P,T}, QUfourier{P,T}, QUmap{P,T}}

############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:S2Field{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end
