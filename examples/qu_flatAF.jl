# ArrayFire S0 and S2 fields

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

#  Define the field types and their trait properties

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################



using ArrayFire #; allowslow(AFArray, false) # dot doesn't work with this set to false
using FieldsBase

# -------------- define the field types ---------------------------
import FieldsBase: has_qu, is_map, is_lense_basis

# Tmap w/ArrayFire
struct TmapAF{P<:Flat,T<:Real} <: Field{P,T,S0}
    tx::AFArray{T,2}
    TmapAF{P,T}(tx::Matrix) where {P<:Flat,T<:Real} = new{P,T}(AFArray(T.(tx)))
    TmapAF{P,T}(tx::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tx)
end
has_qu(::Type{TmapAF{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{TmapAF{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{TmapAF{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# Tfourier
struct TfourierAF{P<:Pix,T<:Real} <: Field{P,T,S0}
    tk::AFArray{Complex{T},2}
    TfourierAF{P,T}(tk::Matrix) where {P<:Flat,T<:Real}  = new{P,T}(AFArray(complex.(T.(tk))))
    TfourierAF{P,T}(tk::AFArray) where {P<:Flat,T<:Real} = new{P,T}(tk)
end
has_qu(::Type{TfourierAF{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{TfourierAF{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}



# QUmap w/ArrayFire
struct QUmapAF{P<:Flat,T<:Real} <: Field{P,T,S2}
    qx::AFArray{T,2}
    ux::AFArray{T,2}
    QUmapAF{P,T}(qx::AFArray, ux::AFArray) where {P<:Flat,T<:Real} = new{P,T}(qx,ux)
    QUmapAF{P,T}(qx::Matrix,  ux::Matrix)  where {P<:Flat,T<:Real} = new{P,T}(AFArray(T.(qx)), AFArray(T.(ux)))
end
has_qu(::Type{QUmapAF{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUmapAF{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}
is_lense_basis(::Type{QUmapAF{P,T}}) where {P<:Flat,T<:Real} = IsLenseBasis{true}


# QUfourier w/ArrayFire
struct QUfourierAF{P<:Pix,T<:Real} <: Field{P,T,S2}
    qk::AFArray{Complex{T},2}
    uk::AFArray{Complex{T},2}
    QUfourierAF{P,T}(qk::AFArray, uk::AFArray) where {P<:Flat,T<:Real} = new{P,T}(qk, uk)
    QUfourierAF{P,T}(qk::Matrix,  uk::Matrix)  where {P<:Flat,T<:Real} = new{P,T}(AFArray(Complex{T}.(qk)), AFArray(Complex{T}.(uk)))
end
has_qu(::Type{QUfourierAF{P,T}}) where {P<:Flat,T<:Real} = HasQU{true}
is_map(::Type{QUfourierAF{P,T}}) where {P<:Flat,T<:Real} = IsMap{false}


# EBmap w/ArrayFire
struct EBmapAF{P<:Pix, T<:Real} <: Field{P,T,S2}
    ex::AFArray{T,2}
    bx::AFArray{T,2}
    EBmapAF{P,T}(ex::AFArray, bx::AFArray) where {P<:Flat,T<:Real} = new{P,T}(ex, bx)
    EBmapAF{P,T}(ex::Matrix,  bx::Matrix)  where {P<:Flat,T<:Real} = new{P,T}(AFArray(T.(ex)), AFArray(T.(bx)))
end
has_qu(::Type{EBmapAF{P,T}}) where {P<:Flat,T<:Real} = HasQU{false}
is_map(::Type{EBmapAF{P,T}}) where {P<:Flat,T<:Real} = IsMap{true}


# EBfourier w/ArrayFire
struct EBfourierAF{P<:Pix, T<:Real} <: Field{P,T,S2}
    ek::AFArray{Complex{T},2}
    bk::AFArray{Complex{T},2}
    EBfourierAF{P,T}(ek::AFArray, bk::AFArray) where {P<:Flat,T<:Real} = new{P,T}(ek, bk)
    EBfourierAF{P,T}(ek::Matrix, bk::Matrix) where {P<:Flat,T<:Real}   = new{P,T}(AFArray(Complex{T}.(ek)), AFArray(Complex{T}.(bk)))
end
has_qu(::Type{EBfourierAF{P,T}}) where {P<:Flat,T<:Real}  = HasQU{false}
is_map(::Type{EBfourierAF{P,T}}) where {P<:Flat,T<:Real}  = IsMap{false}



# -----------------Abstract Fields for dispatch ------------------------
const S0FieldAF{P,T} = Union{TmapAF{P,T}, TfourierAF{P,T}}
const S2FieldAF{P,T} = Union{EBfourierAF{P,T}, EBmapAF{P,T}, QUfourierAF{P,T}, QUmapAF{P,T}}
const FieldAF{P,T}  = Union{S0FieldAF{P,T}, S2FieldAF{P,T}}


# ---------------------harmonic_transform ----------------------
import FieldsBase: harmonic_transform
function harmonic_transform(::Type{F}) where F<:FieldAF{P,T} where {P<:Flat, T<:Real}
    return AFr𝔽(P,T)    #<--- return an instance 
end


#=
This is basically all that needs to be done for these new field types. 
The traits take care of all the conversion, promotion, field algebra, DiagLinOps. 
However, since we will be using a non-standard FFT we need to additionally 
define a specialized ArrayFire FFT object AFr𝔽{P,T}
=#


############################################################
############################################################
############################################################
############################################################
############################################################

# Specialized ArrayFire FFT since it is not pre-defined in FieldsBase. 

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################

#  ArrayFire FFT
struct AFr𝔽{P<:Flat,T<:Real} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωx::T
    period::T
    nyq::T
    k::Vector{AFArray{T,2}}
    x::Vector{AFArray{T,2}}
    sin2ϕk::AFArray{T,2}
    cos2ϕk::AFArray{T,2}
end

# real FFT generated function constructor
@generated function AFr𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωx     = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [AFArray(T.(reshape(k_side, 1, nside))), AFArray(T.(reshape(k_side[1:nside÷2+1], nside÷2+1, 1)))]
    x      = [AFArray(T.(reshape(x_side, 1, nside))), AFArray(T.(reshape(x_side, nside, 1)))]
    ϕk     = AFArray(T.( atan2.(Array(k[2]), Array(k[1]) )))
    AFr𝔽{P,T}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk))
end

import Base: *, \
(*)(g::AFr𝔽{P,T}, x) where {P<:Flat,T}  = rfft2(x, T(g.Ωx/2/π))
(\)(g::AFr𝔽{P,T}, x) where {P<:Flat{θ,n},T} where {θ,n} = irfft2(x, T(2*π/g.Ωx/n^2))




############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
#  The fields are ready to go ...
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################




using Base.Test #... or `using Test` in v0.7

nside  = 1024
Θpix   = 2.0
P      = Flat{Θpix,nside}
T      = Float32
g      =  AFr𝔽(P,T)

qx, ux = AFArray(rand(T, nside, nside)), AFArray(rand(T, nside, nside))
ex, bx = AFArray(rand(T, nside, nside)), AFArray(rand(T, nside, nside))
qk, uk = AFArray(rand(Complex{T}, nside÷2+1, nside)), AFArray(rand(Complex{T}, nside÷2+1, nside))
ek, bk = AFArray(rand(Complex{T}, nside÷2+1, nside)), AFArray(rand(Complex{T}, nside÷2+1, nside))

p1 = QUmapAF{P,T}(qx, ux)
p2 = EBmapAF{P,T}(ex, bx)
p3 = QUfourierAF{P,T}(qk, uk)
p4 = EBfourierAF{P,T}(ek, bk)

p = convert(QUfourierAF{P,T}, p1)
@test all(p.qk - g * p1.qx .== 0)

@inferred convert(EBfourierAF{P,T}, p1)
@inferred convert(EBfourierAF{P,T}, p3)
@inferred convert(QUfourierAF{P,T}, p3)

@inferred QUfourierAF{P,T}(p1)
@inferred QUfourierAF{P,T}(p2)
@inferred QUfourierAF{P,T}(p3)
@inferred QUfourierAF{P,T}(p4)

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
p1, p2 = QUmapAF{P,T}(r1, r2), QUmapAF{P,T}(r3, r4)

@test dot(p1, p2) == dot(p2, p1)
#@test dot(p1, p2) == (dot(r1,r3) + dot(r2,r4))*r𝔽(P,T).Ωx  #<-- for ArrayFire these sometimes fail due to a discrepancy withBLAS
@test dot(p1, p2) == (sum(r1.*r3) + sum(r2.*r4))*r𝔽(P,T).Ωx #<-- for ArrayFire these sometimes fail due to a discrepancy withBLAS
@test dot(p1, p1) > 0

@inferred dot(p1, p2)
@inferred dot(EBfourierAF{P,T}(p1), EBfourierAF{P,T}(p2))
@inferred dot(QUfourierAF{P,T}(p1), p2)
@inferred dot(QUfourierAF{P,T}(p1), QUfourierAF{P,T}(p2))
@inferred dot(EBfourierAF{P,T}(p1), p2)
@inferred dot(p1, EBfourierAF{P,T}(p2))
@inferred dot(QUmapAF{P,T}(p1), p2)
@inferred dot(p1, QUmapAF{P,T}(p2))
@inferred dot(QUmapAF{P,T}(p1), QUmapAF{P,T}(p2))
@inferred dot(p1, QUmapAF{P,T}(QUfourierAF{P,T}(p2)))

@inferred dot(p1, EBmapAF{P,T}(p2))
@inferred dot(EBmapAF{P,T}(p1), EBmapAF{P,T}(p2))
@inferred dot(EBmapAF{P,T}(p1), p2)


wn1 = white_noise(g)
wn2 = white_noise(g)
p = QUmapAF{P,T}(wn1, wn2)
dot(p, p)/2/nside^2 # this should be near 1
dot(EBfourierAF{P,T}(p), EBfourierAF{P,T}(p))/2/nside^2 # this should be near 1


##### Testing DiagOp
p1 = QUmapAF{P,T}(qx, ux)
p2 = EBmapAF{P,T}(ex, bx)
p3 = QUfourierAF{P,T}(qk, uk)
p4 = EBfourierAF{P,T}(ek, bk)

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




tx = AFArray(rand(T, nside, nside))
tk = AFArray(rand(Complex{T}, nside÷2+1, nside))

t1 = TmapAF{P,T}(tx)
t2 = TfourierAF{P,T}(tk)

t = convert(TfourierAF{P,T}, t1)
@test all(t.tk - g * t1.tx .== 0)

@test typeof(convert(TfourierAF{P,T}, t1)) == TfourierAF{P,T}
@test typeof(convert(TfourierAF{P,T}, t2)) == TfourierAF{P,T}
@test typeof(convert(TmapAF{P,T}, t1)) == TmapAF{P,T}
@test typeof(convert(TmapAF{P,T}, t2)) == TmapAF{P,T}

@inferred convert(TfourierAF{P,T}, t1)
@inferred convert(TfourierAF{P,T}, t2)
@inferred convert(TmapAF{P,T}, t1)
@inferred convert(TmapAF{P,T}, t2)

@inferred 2 * t1 - 5.0 * t1
@inferred 2 * t1 - 5.0 * t2
@inferred 2 * t2 - 5.0 * t1
@inferred 2 * t2 - 5.0 * t2

##### Testing dot
t1x = rand(T, nside, nside)
t2x = rand(T, nside, nside)
t1 = TmapAF{P,T}(t1x)
t2 = TmapAF{P,T}(t2x)

@test dot(t1, t2) == dot(t2, t1)
@test dot(t1, t2) == sum(t1x.*t2x)*AFr𝔽(P,T).Ωx
@test dot(t1, t1) > 0








############################################################
############################################################
############################################################
############################################################
############################################################
#  Benchmark Lensing ...
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


struct LenseFlowAF{P<:Flat, T<:Real} 
    ϕ::TfourierAF{P,T}
    tstart::T
    tstop::T
    nsteps::Int
    p1t::Vector{AFArray{T,2}}
    p2t::Vector{AFArray{T,2}}
    fulltspan::Array{T,1} # <-- note: not stored on the GPU
    ik1::AFArray{Complex{T},2}
    ik2::AFArray{Complex{T},2}
end

const default_ode_steps = 15

function (::Type{LenseFlowAF{P,T}})(ϕ::S0FieldAF; tstart = 0, tstop = 1, nsteps::Int = default_ode_steps) where {P<:Flat, T<:Real}
    p1t, p2t, fulltspan, ik1, ik2 = get_pt(ϕ, T(tstart), T(tstop), nsteps)
    return LenseFlowAF{P,T}(ϕ, tstart, tstop, nsteps, p1t, p2t, fulltspan, ik1, ik2)
end

function *(L::LenseFlowAF{P,T}, p::A) where A<:S2FieldAF{P,T} where {P<:Flat{θpix,nside}, T<:Real} where {θpix,nside}
    if L.tstart == L.tstop
        return p
    end
    qx, ux = p |> QUmapAF{P,T} |> data
    lnqx = lense!(qx, L)
    lnux = lense!(ux, L)
    return A(QUmapAF{P,T}(lnqx,lnux))
end

function \(L::LenseFlowAF{P,T}, p::S2FieldAF{P,T}) where {P<:Flat, T<:Real}
    return LenseFlowAF{P,T}(L.ϕ, L.tstop, L.tstart, L.nsteps,L.p1t, L.p2t, L.fulltspan, L.ik1, L.ik2) * p
end

function lense!(x_init::AFArray{T,2}, L) where T<:Real # <--- f(t, x)
    tspan = [L.tstart, L.tstop]
    nsteps = L.nsteps
    ϵ = T((tspan[end] - tspan[1]) / nsteps)
    newtspan = T[T(tspan[1] + n*ϵ) for n in 0:nsteps-1]
    x      = deepcopy(x_init)
    v1, v2 = zeros(AFArray{T,2},size(x_init)), zeros(AFArray{T,2},size(x_init))
    v3, v4 = zeros(AFArray{T,2},size(x_init)), zeros(AFArray{T,2},size(x_init))
    rk1, rk2, rk3, rk4 = T(1/6), T(2/6), T(2/6), T(1/6)
    for i = 1:nsteps
        𝒱_lens_flow_array_AF!(newtspan[i], x, v1, L)
        sync(v1)
        𝒱_lens_flow_array_AF!(newtspan[i] + (ϵ/2), x .+ (ϵ/2).*v1, v2, L)
        sync(v2)
        𝒱_lens_flow_array_AF!(newtspan[i] + (ϵ/2), x .+ (ϵ/2).*v2, v3, L)
        sync(v3)
        𝒱_lens_flow_array_AF!(newtspan[i] +   (ϵ), x .+   (ϵ).*v3, v4, L)
       x .+= ϵ.*rk1.*v1 .+ ϵ.*rk2.*v2 .+ ϵ.*rk3.*v3 .+ ϵ.*rk4.*v4
       sync(x)
    end
    return x
end

function 𝒱_lens_flow_array_AF!(tm, fx::AFArray{T,2}, fx_out::AFArray{T,2}, L::LenseFlowAF{P,T})  where {P<:Flat{θ,nside},T<:Real} where {θ,nside}
    g = AFr𝔽(P,T)
    tmi = findfirst(tm - T(1e-6) .<=  L.fulltspan .<= tm + T(1e-6) )
    fk_strg   = g * fx # A_mul_B!(fk_strg, prfft, fx)
    sync(fk_strg)
    ∂1fx_strg = g \ (L.ik1 .* fk_strg)
    ∂2fx_strg = g \ (L.ik2 .* fk_strg)
    sync(∂1fx_strg); sync(∂2fx_strg)
    fx_out .= L.p1t[tmi] .* ∂1fx_strg .+ L.p2t[tmi] .* ∂2fx_strg
    sync(fx_out)
    return nothing
end

function get_pt(ϕ::S0FieldAF{P,T}, tstart, tstop, nsteps) where {P<:Flat,T<:Real}
    ϵ = T(abs(tstop - tstart) / nsteps)
    fulltspan = T[T(min(tstart,tstop) + n*ϵ/2) for n in 0:2nsteps]
    g = AFr𝔽(P,T)
    zk1 = 0 .* g.k[1]
    zk2 = 0 .* g.k[2]
    ik1 = AFArray(im .* Array(g.k[1]) .+ Array(zk2)) # why do I need both Array converts for it to work? ... it gives an error otherwise
    ik2 = AFArray(im .* Array(g.k[2]) .+ Array(zk1))
    #---------------------------
    ϕk,   = ϕ |> TfourierAF{P,T} |> data
    ∂1ϕ   = g \ (ik1 .* ϕk) 
    ∂2ϕ   = g \ (ik2 .* ϕk) 
    ∂11ϕ  = g \ (ik1 .* ik1 .* ϕk) 
    ∂12ϕ  = g \ (ik1 .* ik2 .* ϕk) 
    ∂22ϕ  = g \ (ik2 .* ik2 .* ϕk) 
    #---------------------------
    p1t = Vector{AFArray{T,2}}(length(fulltspan))
    p2t = Vector{AFArray{T,2}}(length(fulltspan))
    M11  = similar(∂11ϕ) 
    M12  = similar(∂11ϕ) 
    M22  = similar(∂11ϕ) 
    detM = similar(∂11ϕ) 
    @inbounds for i in 1:length(fulltspan)
        tm = fulltspan[i]
        M11  .= 1 .+ tm .* ∂22ϕ
        M12  .=   .- tm .* ∂12ϕ
        M22  .= 1 .+ tm .* ∂11ϕ
        detM .= M11 .* M22 .- M12.*M12
        sync(M11); sync(M12); sync(M22)
        p1t[i]  = (M11 .* ∂1ϕ .+ M12 .* ∂2ϕ)./detM
        p2t[i]  = (M12 .* ∂1ϕ .+ M22 .* ∂2ϕ)./detM
    end
    return p1t, p2t, fulltspan, ik1, ik2
end





# ----------------------- now compare
import Lense 
using JLD2, FileIO, PyPlot
saved_cl_FILE  = joinpath(dirname(Lense.source_path), "saved_cls/cls_r_0.05.jld2") 

##############  set grid geometry #########
nside = 1024 # 
Θpix  = 2.0 # 2.0
P     = Flat{Θpix,nside}
T     = Float32
g     = Lense.r𝔽(P,T)
gAF   = AFr𝔽(P,T)

################ regular array fields #########################
cl      = Lense.unlensedCl(P, T, load(saved_cl_FILE,"cls"))
lncl    = Lense.lensedCl(P, T,   load(saved_cl_FILE,"cls"))
t, p, ϕ = Lense.simulate(cl)
Lϕ      = Lense.LenseFlow{P,T}(ϕ, tstart=0, tstop=1, nsteps=15)
Lϕ⁻¹    = Lense.LenseFlow{P,T}(ϕ, tstart=1, tstop=0, nsteps=15)


################ regular array fields #########################
pAF    = EBfourierAF{P,T}(AFArray(p[:ek]), AFArray(p[:bk]))
ϕAF    = TfourierAF{P,T}(AFArray(ϕ[:tk]))
LϕAF   = LenseFlowAF{P,T}(ϕAF, tstart=0, tstop=1, nsteps=15)
LϕAF⁻¹ = LenseFlowAF{P,T}(ϕAF, tstart=1, tstop=0, nsteps=15)



############### basic test 
@time qx    = Lϕ   * p   |> Lense.QUmap{P,T}   |> x->Lense.data(x)[1];
@time qxAF  = LϕAF * pAF |>       QUmapAF{P,T} |> x->data(x)[1] |> Array;

matshow(qxAF - qx)
matshow(qxAF)
matshow(qx)




