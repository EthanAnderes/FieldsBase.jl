
if VERSION >= v"0.7.0-DEV.1"
    using FFTW
end
FFTW.set_num_threads(6)
BLAS.set_num_threads(6)
############################################################
#  load user defined Field types
############################################################
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


#NOTE: Need to choose a fourier or spherical harmonic transform
const MyField{Px,Tx} = Union{EBfourier{Px,Tx}, EBmap{Px,Tx}, QUfourier{Px,Tx}, QUmap{Px,Tx}}
function harmonic_transform(::Type{F}) where F<:MyField{Px,Tx} where {Px<:Flat, Tx<:Real}
    return r𝔽(Px,Tx)
end




############################################################
#  set grid geometry, etc.
############################################################

nside  = 512
Θpix   = 2.0
Px     = Flat{Θpix,nside}
#Tx     = Float64
Tx     = Float32
#odes   = Lense.default_ode_steps  # number of ODE steps
g      =  r𝔽(Px,Tx);
#kmag   = sqrt.(g.k[1].^2 .+ g.k[2].^2)
#Δ⁻¹    = I𝔽(-squash.(1./kmag.^2), Px, Tx)
#Δ      = I𝔽(-kmag.^2, Px, Tx)

# Test ....
qx, ux = rand(Tx, nside, nside), rand(Tx, nside, nside)
ex, bx = rand(Tx, nside, nside), rand(Tx, nside, nside)
qk, uk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)
ek, bk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)

p1 = QUmap{Px,Tx}(qx, ux)
p2 = EBmap{Px,Tx}(ex, bx)
p3 = QUfourier{Px,Tx}(qk, uk)
p4 = EBfourier{Px,Tx}(ek, bk)

p = convert(QUfourier{Px,Tx}, p1)
p.qk - g * p1.qx

# using BenchmarkTools
# @benchmark convert(QUfourier{Px,Tx}, $p4)
# @benchmark  convert(QUfourier{Px,Tx}, $p1)
@time convert(QUfourier{Px,Tx}, p4) # 0.0021
@time  convert(QUfourier{Px,Tx}, p1) # 0.000805 seconds



p = convert(EBfourier{Px,Tx}, p1)
p = convert(EBfourier{Px,Tx}, p3)
p = convert(QUfourier{Px,Tx}, p3)



using Base.Test

@inferred QUfourier{Px,Tx}(p1)
@inferred QUfourier{Px,Tx}(p2)
@inferred QUfourier{Px,Tx}(p3)
@inferred QUfourier{Px,Tx}(p4)

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
