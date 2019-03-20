
using FieldsBase
using Test
import LinearAlgebra: I



########################################################
#     test predefined template qu_flat.jl
########################################################
include(joinpath(FieldsBase.module_dir, "templates/qu_flat.jl"))

nside  = 512
Θpix   = 2.0
Px     = Flat{Θpix,nside}
Tx     = Float64
g      =  r𝔽(Px,Tx)

qx, ux = rand(Tx, nside, nside), rand(Tx, nside, nside)
ex, bx = rand(Tx, nside, nside), rand(Tx, nside, nside)
qk, uk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)
ek, bk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)

ek, bk = g * (qx, ux)
qx, ux = g \ (ek, bk)

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


##### Testing Uniform scaling
for pp in [p1, p2, p3, p4]
  @inferred (2I) * pp
  @inferred pp * (2I) 
  for (d1, d2, d3) in zip(data((2I) * pp), data(pp * (2I)), data(2pp))
    @test all(d1 .== d2 .== d3)
  end
end 

##### Testing zero constructor 
@inferred QUmap{Px,Tx}()
@inferred EBmap{Px,Tx}()
@inferred QUfourier{Px,Tx}()
@inferred EBfourier{Px,Tx}()

p1 = QUmap{Px,Tx}()
p2 = EBmap{Px,Tx}()
p3 = QUfourier{Px,Tx}()
p4 = EBfourier{Px,Tx}()

for pp in (p1, p2, p3, p4)
  for ff in data(pp)
    @test all(ff .== 0)
  end
end

##### Testing dot
r1, r2, r3, r4 = randn(Tx,nside, nside), randn(Tx,nside, nside), randn(Tx,nside, nside), randn(Tx,nside, nside)
p1, p2 = QUmap{Px,Tx}(r1, r2), QUmap{Px,Tx}(r3, r4)

@test dot(p1, p2) == dot(p2, p1)
#@test dot(p1, p2) == (dot(r1,r3) + dot(r2,r4))*r𝔽(Px,Tx).Ωx
@test dot(p1, p2) == (sum(r1.*r3) + sum(r2.*r4))*r𝔽(Px,Tx).Ωx
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


wn1 = white_noise(g)
wn2 = white_noise(g)
p = QUmap{Px,Tx}(wn1, wn2)
dot(p, p)/2/nside^2 # this should be near 1
dot(EBfourier{Px,Tx}(p), EBfourier{Px,Tx}(p))/2/nside^2 # this should be near 1


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



############################################################
#  Test predefined template t_flat.jl
############################################################
include(joinpath(FieldsBase.module_dir, "templates/t_flat.jl"))


nside  = 512
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


##### Testing Uniform scaling
for pp in [t1, t2]
  @inferred (2I) * pp
  @inferred pp * (2I) 
  for (d1, d2, d3) in zip(data((2I) * pp), data(pp * (2I)), data(2pp))
    @test all(d1 .== d2 .== d3)
  end
end 

##### Testing zero constructor 
@inferred QUmap{Px,Tx}()
@inferred EBmap{Px,Tx}()

t1 = QUmap{Px,Tx}()
t2 = EBmap{Px,Tx}()

for pp in (t1, t2)
  for ff in data(pp)
    @test all(ff .== 0)
  end
end



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
dot(Tfourier{Px,Tx}(t), Tfourier{Px,Tx}(t))/nside^2 # this should be near 1


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





############################################################
#  Test predefined template tqu_flat.jl
############################################################
include(joinpath(FieldsBase.module_dir, "templates/tqu_flat.jl"))

nside  = 512
Θpix   = 2.0
P     = Flat{Θpix,nside}
#T     = Float64
T     = Float32
g     =  r𝔽(P,T);

# Test ....
tx, qx, ux = rand(T, nside, nside), rand(T, nside, nside), rand(T, nside, nside)
tx, ex, bx = rand(T, nside, nside), rand(T, nside, nside), rand(T, nside, nside)
tk, ek, bk = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)
tk, qk, uk = rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside), rand(Complex{T}, nside÷2+1, nside)

p1 = TQUmap{P,T}(tx, qx, ux)
p2 = TEBfourier{P,T}(tk, ek, bk)
p3 = TQUfourier{P,T}(tx, qk, uk)
p4 = TEBmap{P,T}(tx, ex, bx)


@time convert(TEBfourier{P,T}, p1)
@time convert(TEBfourier{P,T}, p2)



p = convert(TEBfourier{P,T}, p1)
p = convert(TEBfourier{P,T}, p2)
p = convert(TQUmap{P,T}, p2)


@inferred TEBfourier{P,T}(p1)
@inferred TEBfourier{P,T}(p2)
@inferred TQUmap{P,T}(p1)
@inferred TQUmap{P,T}(p2)

@inferred 2 * p1 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2
@inferred 2 * p2 - 5.0 * p1
@inferred 2 * p1 - 5.0 * p2


##### Testing Uniform scaling
for pp in [p1, p2, p3, p4]
  @inferred (2I) * pp
  @inferred pp * (2I) 
  for (d1, d2, d3) in zip(data((2I) * pp), data(pp * (2I)), data(2pp))
    @test all(d1 .== d2 .== d3)
  end
end 

##### Testing zero constructor 
@inferred TQUmap{P,T}()
@inferred TEBmap{P,T}()
@inferred TQUfourier{P,T}()
@inferred TEBfourier{P,T}()

tp1 = TQUmap{P,T}()
tp2 = TEBmap{P,T}()
tp3 = TQUfourier{P,T}()
tp4 = TEBfourier{P,T}()

for pp in (tp1, tp2, tp3, tp4)
  for ff in data(pp)
    @test all(ff .== 0)
  end
end




@test dot(p1, p2) == dot(p2, p1)
@test dot(p1, p1) > 0

@inferred dot(p1, p2)

wn1, wn2, wn3 = white_noise(g), white_noise(g), white_noise(g)
p = TQUmap{P,T}(wn1, wn2, wn3)
dot(p, p) / 3 / nside^2 # this should be near 1.0
dot(TEBfourier{P,T}(p), TEBfourier{P,T}(p)) / 3 / nside^2 # this should be near 1






#=

using BenchmarkTools

nside  = 512
Θpix   = 2.0
Px     = Flat{Θpix,nside}
Tx     = Float32
g      = r𝔽(Px,Tx);

qx, ux = rand(Tx, nside, nside), rand(Tx, nside, nside)
ex, bx = rand(Tx, nside, nside), rand(Tx, nside, nside)
qk, uk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)
ek, bk = rand(Complex{Tx}, nside÷2+1, nside), rand(Complex{Tx}, nside÷2+1, nside)

p1 = QUmap{Px,Tx}(qx, ux)
p2 = QUfourier{Px,Tx}(qk, uk)

@benchmark convert($(EBfourier{Px,Tx}), $p1)


# =============
@benchmark convert($(EBfourier{Px,Tx}), $p2)

# on v0.6.3 with fused loop
BenchmarkTools.Trial:
  memory estimate:  4.02 MiB
  allocs estimate:  10
  --------------
  minimum time:     717.703 μs (0.00% GC)
  median time:      903.073 μs (0.00% GC)
  mean time:        1.457 ms (11.99% GC)
  maximum time:     4.295 ms (23.42% GC)
  --------------
  samples:          3426
  evals/sample:     1


# on v0.6.3 without fused loop
BenchmarkTools.Trial:
  memory estimate:  4.02 MiB
  allocs estimate:  12
  --------------
  minimum time:     1.284 ms (0.00% GC)
  median time:      1.397 ms (0.00% GC)
  mean time:        1.950 ms (8.41% GC)
  maximum time:     4.456 ms (20.11% GC)
  --------------
  samples:          2560
  evals/sample:     1


# on v0.7 without fused loop
BenchmarkTools.Trial:
  memory estimate:  4.02 MiB
  allocs estimate:  24
  --------------
  minimum time:     798.389 μs (0.00% GC)
  median time:      975.551 μs (0.00% GC)
  mean time:        1.523 ms (11.92% GC)
  maximum time:     42.914 ms (92.72% GC)
  --------------
  samples:          3276
  evals/sample:     1


# on v0.7 with fused loop
julia>   @benchmark convert($(EBfourier{Px,Tx}), $p2)
BenchmarkTools.Trial: 
  memory estimate:  4.02 MiB
  allocs estimate:  9
  --------------
  minimum time:     549.295 μs (0.00% GC)
  median time:      788.903 μs (0.00% GC)
  mean time:        1.276 ms (14.40% GC)
  maximum time:     4.904 ms (28.45% GC)
  --------------
  samples:          3906
  evals/sample:     1



=#
