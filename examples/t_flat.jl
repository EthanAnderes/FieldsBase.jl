
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

# Tmap
struct Tmap{Px<:Flat,Tx<:Real} <: Field{Px,Tx,S0}
    tx::Matrix{Tx}
    Tmap{Px,Tx}(tx::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(tx)
end
has_qu(::Type{Tmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = HasQU{false}
is_map(::Type{Tmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{true}
is_lense_basis(::Type{Tmap{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsLenseBasis{true}


# Tfourier
struct Tfourier{Px<:Pix,Tx<:Real} <: Field{Px,Tx,S0}
    tk::Matrix{Complex{Tx}}
    Tfourier{Px,Tx}(tk::Matrix) where {Px<:Flat,Tx<:Real} = new{Px,Tx}(complex.(tk))
end
has_qu(::Type{Tfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real} = HasQU{false}
is_map(::Type{Tfourier{Px,Tx}}) where {Px<:Flat,Tx<:Real} = IsMap{false}


const MyField{Px,Tx} = Union{Tfourier{Px,Tx}, Tmap{Px,Tx}}


############################################################
#  Specify the harmonic transform
############################################################

function harmonic_transform(::Type{F}) where F<:MyField{Px,Tx} where {Px<:Flat, Tx<:Real}
    return r𝔽(Px,Tx)
end


############################################################
#  The fields are ready to go ...
############################################################


using Base.Test
# using Test # for v0.7

nside  = 1024
Θpix   = 2.0
Px     = Flat{Θpix,nside}
Tx     = Float32
g      =  r𝔽(Px,Tx);

tx = rand(Tx, nside, nside)
tk = rand(Complex{Tx}, nside÷2+1, nside)

t1 = Tmap{Px,Tx}(tx)
t2 = Tfourier{Px,Tx}(tk)

t = convert(Tfourier{Px,Tx}, t1)
@test all(t.tk - g * t1.tx .== 0)

@test typeof(convert(Tfourier{Px,Tx}, t1)) == Tfourier{Px,Tx}
@test typeof(convert(Tfourier{Px,Tx}, t2)) == Tfourier{Px,Tx}
@test typeof(convert(Tmap{Px,Tx}, t1)) == Tmap{Px,Tx}
@test typeof(convert(Tmap{Px,Tx}, t2)) == Tmap{Px,Tx}

@inferred convert(Tfourier{Px,Tx}, t1)
@inferred convert(Tfourier{Px,Tx}, t2)
@inferred convert(Tmap{Px,Tx}, t1)
@inferred convert(Tmap{Px,Tx}, t2)

@inferred 2 * t1 - 5.0 * t1
@inferred 2 * t1 - 5.0 * t2
@inferred 2 * t2 - 5.0 * t1
@inferred 2 * t2 - 5.0 * t2

@time 2 * t1 - 5.0 * t2
function foo(t1, t2)
    a = 2 * t1 - 5 * t2
    b = 2 * t1 - 5 * a
    c = a-b
    return c
end
@time foo(t1, t2) 


##### Testing dot
t1x = rand(Tx, nside, nside)
t2x = rand(Tx, nside, nside)
t1 = Tmap{Px,Tx}(t1x)
t2 = Tmap{Px,Tx}(t2x)

@test dot(t1, t2) == dot(t2, t1)
@test dot(t1, t2) == sum(t1x.*t2x)*r𝔽(Px,Tx).Ωx
@test dot(t1, t1) > 0

@inferred dot(t1, t2)
@inferred dot(Tfourier{Px,Tx}(t1), Tfourier{Px,Tx}(t2))
@inferred dot(Tfourier{Px,Tx}(t1), t2)
@inferred dot(t1, Tfourier{Px,Tx}(t2))
@inferred dot(Tmap{Px,Tx}(t1), t2)
@inferred dot(t1, Tmap{Px,Tx}(t2))
@inferred dot(Tmap{Px,Tx}(t1), Tmap{Px,Tx}(t2))


wn1 = white_noise(g)
t = Tmap{Px,Tx}(wn1)
dot(t, t)/nside^2 # this should be near 1.0
dot(Tfourier{Px,Tx}(t), Tfourier{Px,Tx}(t))/nside^2 # this should be near nside^2


##### Testing DiagOp
t1 = Tmap{Px,Tx}(tx)
t2 = Tfourier{Px,Tx}(tk)

@inferred 𝕃(t1)*t1
@inferred 𝕃(t1)*t2
@inferred 𝕃(t2)*t1
@inferred 𝕃(t2)*t2


L1 = 𝕃(t1)^(-5.3)
L2 = 𝕃(t2)^(2)
L3 = 𝕃(t1)^(-.1)
L4 = 𝕃(t2)^(4)
L5 = 𝕃(t1)^(1)
L6 = 𝕃(t2)^(0)
L7 = 𝕃(t1)^(-1)
L8 = 𝕃(t2)^(1.5)
L9 = inv(𝕃(t1))
L10 = inv(𝕃(t2))
L11 = 𝕃(t2)^(-1)

@inferred L1*t1
@inferred L2*t1
@inferred L3*t2
@inferred L4*t2
@inferred L5*t1
@inferred L6*t1
@inferred L7*t2
@inferred L8*t2
@inferred L9*t2 - L7*t2
@inferred L9*t1 - L7*t1
@inferred L10*t1 - L11*t1
