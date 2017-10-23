using FieldsBase
using Base.Test
# using Test # for v0.7


########################################################
#     test predefined template qu_flat.jl
########################################################
# TODO figure out how to get Pkg.dir("FieldsBase") to point to local directory
# include(joinpath(Pkg.dir("FieldsBase"), "templates", "qu_flat.jl"))
include(joinpath("/Users/ethananderes/Dropbox/FieldsBase", "templates", "qu_flat.jl"))


nside  = 512
Θpix   = 2.0
Px     = Flat{Θpix,nside}
Tx     = Float32
g      =  r𝔽(Px,Tx);

qx, ux = rand(Tx, nside, nside), rand(Tx, nside, nside)
ex, bx = rand(Tx, nside, nside), rand(Tx, nside, nside)
qk, uk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)
ek, bk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)

p1 = QUmap{Px,Tx}(qx, ux)
p2 = EBmap{Px,Tx}(ex, bx)
p3 = QUfourier{Px,Tx}(qk, uk)
p4 = EBfourier{Px,Tx}(ek, bk)

p = convert(QUfourier{Px,Tx}, p1)
@test all(p.qk - g * p1.qx .== 0)

@inferred convert(EBfourier{Px,Tx}, p1)
@inferred convert(EBfourier{Px,Tx}, p3)
@inferred convert(QUfourier{Px,Tx}, p3)

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




##### Testing dot
r1, r2, r3, r4 = randn(Tx,nside, nside), randn(Tx,nside, nside), randn(Tx,nside, nside), randn(Tx,nside, nside)
p1, p2 = QUmap{Px,Tx}(r1, r2), QUmap{Px,Tx}(r3, r4)

@test dot(p1, p2) == dot(p2, p1)
@test dot(p1, p2) == (dot(r1,r3) + dot(r2,r4))*r𝔽(Px,Tx).Ωpix
@test dot(p1, p1) > 0

@inferred dot(p1, p2)
@inferred dot(EBfourier{Px,Tx}(p1), EBfourier{Px,Tx}(p2))
@inferred dot(QUfourier{Px,Tx}(p1), p2)
@inferred dot(QUfourier{Px,Tx}(p1), QUfourier{Px,Tx}(p2))
@inferred dot(EBfourier{Px,Tx}(p1), p2)
@inferred dot(p1, EBfourier{Px,Tx}(p2))
@inferred dot(QUmap{Px,Tx}(p1), p2)
@inferred dot(p1, QUmap{Px,Tx}(p2))
@inferred dot(QUmap{Px,Tx}(p1), QUmap{Px,Tx}(p2))
@inferred dot(p1, QUmap{Px,Tx}(QUfourier{Px,Tx}(p2)))

@inferred dot(p1, EBmap{Px,Tx}(p2))
@inferred dot(EBmap{Px,Tx}(p1), EBmap{Px,Tx}(p2))
@inferred dot(EBmap{Px,Tx}(p1), p2)


wn1 = white_noise(Px, Tx)
wn2 = white_noise(Px, Tx)
p = QUmap{Px,Tx}(wn1, wn2)
dot(p, p) # this should be near 2nside^2
dot(EBfourier{Px,Tx}(p), EBfourier{Px,Tx}(p)) # this should be near nside^2


##### Testing DiagOp
p1 = QUmap{Px,Tx}(qx, ux)
p2 = EBmap{Px,Tx}(ex, bx)
p3 = QUfourier{Px,Tx}(qk, uk)
p4 = EBfourier{Px,Tx}(ek, bk)

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
