#=
	To load use

	```
	include(joinpath(Pkg.dir("FieldsBase"), "templates", "t_flat.jl"))
	```
=#


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# Tmap
struct Tmap{P<:Flat,T<:Real} <: Field{P,S0}
    tx::Matrix{T}
    Tmap{P,T}(tx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(tx)
end
has_qu(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct Tfourier{P<:Pix,T<:Real} <: Field{P,S0}
    tk::Matrix{Complex{T}}
    Tfourier{P,T}(tk::Matrix) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk))
end
has_qu(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const S0Field{P,T} = Union{Tfourier{P,T}, Tmap{P,T}}

############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:S0Field{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end
