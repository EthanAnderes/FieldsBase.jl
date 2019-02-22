# To load use:
#
# using FieldsBase
# include(joinpath(FieldsBase.module_dir,"templates/t_flat_1d.jl"))


############################################################
#  Define the field types and their trait properties
############################################################

using FieldsBase
import FieldsBase: has_qu, is_map, is_lense_basis, harmonic_transform

# Tmap
struct Tmap{P<:Flat,T<:Real} <: Field{P,T,S0}
    tx::Vector{T}
    Tmap{P,T}(tx::Vector) where {P<:Flat,T<:Real} = new{P,T}(tx)
    Tmap{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_map_dim1(P,T))
end
has_qu(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct Tfourier{P<:Pix,T<:Real} <: Field{P,T,S0}
    tk::Vector{Complex{T}}
    Tfourier{P,T}(tk::Vector) where {P<:Flat,T<:Real} = new{P,T}(complex.(tk))
    Tfourier{P,T}() where {P<:Flat,T<:Real} = new{P,T}(FieldsBase.zero_fourier_dim1(P,T))
end
has_qu(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{Tfourier{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}

const S0Field1d{P,T} = Union{Tfourier{P,T}, Tmap{P,T}}


############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:S0Field1d{P,T} where {P<:Flat, T<:Real}
    return r𝔽1(P,T)
end



############################################################
#  other methods that need extending
############################################################

import FieldsBase: _dot, white_noise
import Base: getindex

function white_noise(g::r𝔽1{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n) ./ sqrt(g.Ωx)
end

# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:S0Field1d{P,T} where {T,P<:Flat{θ,m}} where {θ,m}
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
    sum(map(_complexdot, data(a),data(b))) * FT.Ωk
end

function getindex(f::S0Field1d{P,T}, sym::Symbol) where {P<:Flat, T<:Real}
    (sym == :k)  ? Tfourier{P,T}(f).tk :
    (sym == :x)  ? Tmap{P,T}(f).tx :
    error("index is not defined")
end

