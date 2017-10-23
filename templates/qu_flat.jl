#=
	To load from another package use:

	```
	include(joinpath(Pkg.dir("FieldsBase"), "templates", "qu_flat.jl"))
	```
=#

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# QUmap
struct QUmap{Px<:Flat,Tx<:Real} <: Field{Px,S2}
    qx::Matrix{Tx}
    ux::Matrix{Tx}
    QUmap{Px,Tx}(qx::Matrix, ux::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(qx, ux)
end
has_qu(::Type{QUmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = HasQU{true}
is_map(::Type{QUmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{true}
is_lense_basis(::Type{QUmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsLenseBasis{true}


# QUfourier
struct QUfourier{Px<:Pix,Tx<:Real} <: Field{Px,S2}
    qk::Matrix{Complex{Tx}}
    uk::Matrix{Complex{Tx}}
    QUfourier{Px,Tx}(qk::Matrix, uk::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(complex.(qk), complex.(uk))
end
has_qu(::Type{QUfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real} = HasQU{true}
is_map(::Type{QUfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{false}


# EBmap
struct EBmap{Px<:Pix, Tx<:Real} <: Field{Px,S2}
    ex::Matrix{Tx}
    bx::Matrix{Tx}
    EBmap{Px,Tx}(ex::Matrix, bx::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(ex, bx)
end
has_qu(::Type{EBmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = HasQU{false}
is_map(::Type{EBmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{true}


# EBfourier
struct EBfourier{Px<:Pix, Tx<:Real} <: Field{Px,S2}
    ek::Matrix{Complex{Tx}}
    bk::Matrix{Complex{Tx}}
    EBfourier{Px,Tx}(ek::Matrix, bk::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(complex.(ek), complex.(bk))
end
has_qu(::Type{EBfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real}  = HasQU{false}
is_map(::Type{EBfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{false}


const S2Field{Px,Tx} = Union{EBfourier{Px,Tx}, EBmap{Px,Tx}, QUfourier{Px,Tx}, QUmap{Px,Tx}}
function harmonic_transform(::Type{F}) where F<:S2Field{Px,Tx} where {Px<:Flat, Tx<:Real}
    return r𝔽(Px,Tx)
end
