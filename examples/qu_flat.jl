
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


############################################################
#  Specify the harmonic transform
############################################################

const MyField{P,T} = Union{EBfourier{P,T}, EBmap{P,T}, QUfourier{P,T}, QUmap{P,T}}
function harmonic_transform(::Type{F}) where F<:MyField{P,T} where {P<:Flat, T<:Real}
    return r𝔽(P,T)
end




############################################################
#  The fields are ready to go ...
############################################################

nside  = 512
Θpix   = 2.0
P     = Flat{Θpix,nside}
T     = Float32
g      =  r𝔽(P,T);

qx, ux = rand(T, nside, nside), rand(T, nside, nside)
ex, bx = rand(T, nside, nside), rand(T, nside, nside)
qk, uk = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)
ek, bk = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)

p1 = QUmap{P,T}(qx, ux)
p2 = EBmap{P,T}(ex, bx)
p3 = QUfourier{P,T}(qk, uk)
p4 = EBfourier{P,T}(ek, bk)

p = convert(QUfourier{P,T}, p1)
p.qk - g * p1.qx

@time convert(QUfourier{P,T}, p4)
@time convert(QUfourier{P,T}, p1) 



p = convert(EBfourier{P,T}, p1)
p = convert(EBfourier{P,T}, p3)
p = convert(QUfourier{P,T}, p3)



using Base.Test

@inferred QUfourier{P,T}(p1)
@inferred QUfourier{P,T}(p2)
@inferred QUfourier{P,T}(p3)
@inferred QUfourier{P,T}(p4)

@inferred 2 * p1 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2
@inferred 2 * p1 - 5.0 * p3
@inferred 2 * p1 - 5.0 * p4
@inferred 2 * p2 - 5.0 * p1
@inferred 2 * p2 - 5.0 * p2
@inferred 2 * p2 - 5.0 * p3
@inferred 2 * p2 - 5.0 * p4
@inferred 2 * p3 - 5.0 * p1
@inferred 2 * p3 - 5.0 * p2
@inferred 2 * p3 - 5.0 * p3
@inferred 2 * p3 - 5.0 * p4
@inferred 2 * p4 - 5.0 * p1
@inferred 2 * p4 - 5.0 * p2
@inferred 2 * p4 - 5.0 * p3
@inferred 2 * p4 - 5.0 * p4



##### Testing dot
r1, r2, r3, r4 = randn(T,nside, nside), randn(T,nside, nside), randn(T,nside, nside), randn(T,nside, nside)
p1, p2 = QUmap{P,T}(r1, r2), QUmap{P,T}(r3, r4)

@test dot(p1, p2) == dot(p2, p1)
@test dot(p1, p2) == (dot(r1,r3) + dot(r2,r4))*r𝔽(P,T).Ωpix
@test dot(p1, p1) > 0

@inferred dot(p1, p2)
@inferred dot(EBfourier{P,T}(p1), EBfourier{P,T}(p2))
@inferred dot(QUfourier{P,T}(p1), p2)
@inferred dot(QUfourier{P,T}(p1), QUfourier{P,T}(p2))
@inferred dot(EBfourier{P,T}(p1), p2)
@inferred dot(p1, EBfourier{P,T}(p2))
@inferred dot(QUmap{P,T}(p1), p2)
@inferred dot(p1, QUmap{P,T}(p2))
@inferred dot(QUmap{P,T}(p1), QUmap{P,T}(p2))
@inferred dot(p1, QUmap{P,T}(QUfourier{P,T}(p2)))

@inferred dot(p1, EBmap{P,T}(p2))
@inferred dot(EBmap{P,T}(p1), EBmap{P,T}(p2))
@inferred dot(EBmap{P,T}(p1), p2)


wn1 = white_noise(P, T)
wn2 = white_noise(P, T)
p = QUmap{P,T}(wn1, wn2)
dot(p, p) # this should be near 2nside^2
dot(EBfourier{P,T}(p), EBfourier{P,T}(p)) # this should be near nside^2


##### Testing DiagOp
p1 = QUmap{P,T}(qx, ux)
p2 = EBmap{P,T}(ex, bx)
p3 = QUfourier{P,T}(qk, uk)
p4 = EBfourier{P,T}(ek, bk)

@inferred 𝕃(p1)*p1
@inferred 𝕃(p1)*p2
@inferred 𝕃(p1)*p3
@inferred 𝕃(p1)*p4

@inferred 𝕃(p2)*p1
@inferred 𝕃(p2)*p2
@inferred 𝕃(p2)*p3
@inferred 𝕃(p2)*p4

@inferred 𝕃(p3)*p1
@inferred 𝕃(p3)*p2
@inferred 𝕃(p3)*p3
@inferred 𝕃(p3)*p4

@inferred 𝕃(p4)*p1
@inferred 𝕃(p4)*p2
@inferred 𝕃(p4)*p3
@inferred 𝕃(p4)*p4

L1 = 𝕃(p1)^(-5.3)
L2 = 𝕃(p2)^(2)
L3 = 𝕃(p3)^(-.1)
L4 = 𝕃(p4)^(4)
L5 = 𝕃(p1)^(1)
L6 = 𝕃(p2)^(0)
L7 = 𝕃(p1)^(-1)
L8 = 𝕃(p4)^(1.5)
L9 = inv(𝕃(p1))
L10 = inv(𝕃(p2))
L11 = 𝕃(p2)^(-1)

L1*p4
L2*p3
L3*p2
L4*p1
L5*p2
L6*p2
L7*p2
L8*p2
L9*p2 - L7*p2
L9*p3 - L7*p3
L10*p4 - L11*p4
