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
