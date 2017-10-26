

if VERSION >= v"0.7.0-DEV.1"
    using FFTW
end
FFTW.set_num_threads(6)
BLAS.set_num_threads(6)

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


############################################################
#  The fields are ready to go ...
############################################################


nside  = 512
Θpix   = 2.0
P     = Flat{Θpix,nside}
T     = Float32
g     =  r𝔽(P,T);

# Test ....
tx, qx, ux = rand(T, nside, nside), rand(T, nside, nside), rand(T, nside, nside)
tk, ek, bk = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)

p1 = TQUmap{P,T}(tx, qx, ux)
p2 = TEBfourier{P,T}(tk, ek, bk)

@time convert(TEBfourier{P,T}, p1)
@time convert(TEBfourier{P,T}, p2)



p = convert(TEBfourier{P,T}, p1)
p = convert(TEBfourier{P,T}, p2)
p = convert(TQUmap{P,T}, p2)



using Base.Test

@inferred TEBfourier{P,T}(p1)
@inferred TEBfourier{P,T}(p2)
@inferred TQUmap{P,T}(p1)
@inferred TQUmap{P,T}(p2)

@inferred 2 * p1 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2
@inferred 2 * p2 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2


@test dot(p1, p2) == dot(p2, p1)
@test dot(p1, p1) > 0

@inferred dot(p1, p2)

wn1, wn2, wn3 = white_noise(P, T), white_noise(P, T), white_noise(P, T)
p = TQUmap{P,T}(wn1, wn2, wn3)
dot(p, p) / 3 / nside^2 # this should be near 1.0
dot(TEBfourier{P,T}(p), TEBfourier{P,T}(p)) / 3 / nside^2 # this should be near nside^2
