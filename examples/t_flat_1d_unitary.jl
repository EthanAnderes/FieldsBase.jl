
using PyPlot 
include("/Users/ethananderes/Dropbox/FieldsBase/templates/t_flat_1dimension_unitary.jl")

nside = 2^12
Θpix  = 1/nside
P     = Flat{Θpix,nside}
T     = Float64
g     =  r𝔽1d(P,T);
dm    = 1

#= Note:
For unitary fft, Σ = W* Λ W where Λ = W Σ W* and W W* = I
For stationary Σ(x,y) ≡ Σ(x-y) = (1/nside^(dm/2)) * W diagΛ so to get  1 = Σ(0) need mean(diagΛ) == 1
=#
function trang(x,ρ) 
    y = abs(x)
    return y > ρ ? 0 : 1 - y/ρ
end 
function normalize_Λ(Λk,nside,dm)
    ave_Λk = (2sum(Λk) - Λk[1] - Λk[end]) / (nside^dm)
    Λk   ./= ave_Λk # i.e. want ave to be 1
    return Λk
end


ρ = 0.01 * g.period
ν = 1.4


# -------- spectral densities ------
Maternk     = normalize_Λ(( 4ν/ρ^2 + abs2.(g.k) ) .^ (-ν - dm/2), nside, dm)
Triangk     = normalize_Λ(real.(g * map(x->trang(x,ρ), g.x)), nside, dm)
# plot( (1/nside^(dm/2)) * ( g \ Maternk)) # this is the auto-cov function

# ------ cov operators 
Σν = Maternk  |> Ffourier{P,T} |> 𝕃
Σt = Triangk  |> Ffourier{P,T} |> 𝕃
#= plot 
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
plot((Σν * t_impls)[:fx])
=#



# ---------- whitening operator
# TODO: Try writing Σ = Σ * Δ * Δ⁻¹ = exp(log(Σ * Δ)) * Δ⁻¹
# TODO: Try writing Σ⁻¹ = Δ * (Σ * Δ)⁻¹ 
w_pwr     = 0
w_rng     = 3*g.Δx
#Laplacek  = normalize_Λ(( 4*w_pwr/(w_rng)^2 + abs2.(g.k) ) .^ w_pwr, nside, dm)
Laplacek   = ( 4*w_pwr/(w_rng)^2 + abs2.(g.k) ) .^ w_pwr
# ave_Δ²k    = (2sum(Laplacek.^2) - Laplacek[1]^2 - Laplacek[end]^2) / (nside^dm)
# Laplacek ./= sqrt(ave_Δ²k)
# ave_Δ²k    = (2sum(Maternk .* Laplacek .^ 2) - Maternk[1]*Laplacek[1]^2 - Maternk[end]*Laplacek[end]^2) / (nside^dm)
# Laplacek ./= sqrt(ave_Δ²k)

Δ         = Laplacek |> Ffourier{P,T} |> 𝕃
#= plot 
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
plot((Δ * t_impls)[:fx])
=#


# ---- pre-whiten -----
Δ²Maternk  = Maternk .* Laplacek .^ 2
Δ²Triangk  = Triangk .* Laplacek .^ 2
Δ²Σν = Δ²Maternk  |> Ffourier{P,T} |> 𝕃
Δ²Σt = Δ²Triangk  |> Ffourier{P,T} |> 𝕃

#  ------ log operators ------ 
logΔ²Σν = log.(Δ²Maternk)  |> Ffourier{P,T} |> 𝕃
logΔ²Σt = log.(Δ²Triangk)  |> Ffourier{P,T} |> 𝕃


# ----- simulate ------
sqrtΣν = sqrt.(Maternk)  |> Ffourier{P,T} |> 𝕃 # note: no factor (2π)^(-dm/2)
sqrtΣt = sqrt.(Triangk)  |> Ffourier{P,T} |> 𝕃
sim_ν = sqrtΣν * Fmap{P,T}(white_noise(g))
sim_t = sqrtΣt * Fmap{P,T}(white_noise(g))
# plot(sim_ν[:fx])
# plot(sim_t[:fx])


########################################
#  sparse approx to log
########################################

aprx_logΔ²Σν = spzeros(nside, nside)
aprx_logΔ²Σt = spzeros(nside, nside)
sprs_sz = 50 # approx each col with this number of non-zeros
for ind = 1:nside
    t_impls = (t=zeros(T,nside); t[ind] = 1; t) |> Fmap{P,T}
    logΔ²Σν_col = (logΔ²Σν * t_impls)[:fx]
    logΔ²Σt_col = (logΔ²Σt * t_impls)[:fx]
    pν = sortperm(abs2.(logΔ²Σν_col), rev=true)
    pt = sortperm(abs2.(logΔ²Σt_col), rev=true)
    tmp_approx_ν = logΔ²Σν_col[pν[1:sprs_sz]]
    tmp_approx_t = logΔ²Σt_col[pt[1:sprs_sz]]
    aprx_logΔ²Σν[pν[1:sprs_sz], ind] = tmp_approx_ν
    aprx_logΔ²Σt[pt[1:sprs_sz], ind] = tmp_approx_t
end
# ----- plot --
ind = 2000
t_impls = (t=zeros(T,nside); t[ind] = 1; t) |> Fmap{P,T}
logΔ²Σν_col = (logΔ²Σν * t_impls)[:fx] 
logΔ²Σt_col = (logΔ²Σt * t_impls)[:fx] 
plot(aprx_logΔ²Σν[:,ind])
plot(logΔ²Σν_col)
plot(logΔ²Σν_col - aprx_logΔ²Σν[:,ind])


########################################
#  test the forward map
########################################
# ----- edge 
# t_impls = (t=zeros(T,nside); t[2000:2000+5*sprs_sz] = 1; t) |> Fmap{P,T}
# test_ν      = deepcopy(Δ \ t_impls)
# data_aprx_ν = deepcopy(Δ \ t_impls)[:fx] 
# test_t      = deepcopy(Δ \ t_impls)
# data_aprx_t = deepcopy(Δ \ t_impls)[:fx] 
# --- ipulse
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
test_ν      = deepcopy(Δ \ t_impls)
data_aprx_ν = deepcopy(Δ \ t_impls)[:fx] 
test_t      = deepcopy(Δ \ t_impls)
data_aprx_t = deepcopy(Δ \ t_impls)[:fx] 
# # # --- # simulation 
# test_ν      = deepcopy(Δ \ sim_ν)
# data_aprx_ν = deepcopy(Δ \ sim_ν)[:fx] 
# test_t      = deepcopy(Δ \ sim_t)
# data_aprx_t = deepcopy(Δ \ sim_t)[:fx] 
nsteps = 500; ϵ = 1/500
for n = 0:nsteps-1
    test_t += ϵ * (logΔ²Σt * test_t)
    test_ν += ϵ * (logΔ²Σν * test_ν)
    data_aprx_t += ϵ * (aprx_logΔ²Σt * data_aprx_t)
    data_aprx_ν += ϵ * (aprx_logΔ²Σν * data_aprx_ν)
end
data_aprx_t = data_aprx_t |> Fmap{P,T} |> x->deepcopy(Δ \ x)[:fx]
data_aprx_ν = data_aprx_ν |> Fmap{P,T} |> x->deepcopy(Δ \ x)[:fx]
data_t      = deepcopy(Δ \ test_t)[:fx]
data_ν      = deepcopy(Δ \ test_ν)[:fx] 

figure()
subplot(2,1,1)
plot(data_t)
subplot(2,1,2)
plot(data_aprx_t, alpha=0.2)

figure()
subplot(2,1,1)
plot(data_ν)
subplot(2,1,2)
plot(data_aprx_ν, alpha=0.2)

#=
plot(data_ν)
plot(data_ν-data_aprx_ν)
=#

#=
plot(data_ν[2000:2000+5*sprs_sz] ./ data_aprx_ν[2000:2000+5*sprs_sz])
=#

################################################
#  test the inverse map
################################################

# ----- edge 
# t_impls = (t=zeros(T,nside); t[2000:2000+2*sprs_sz] = 1; t) |> Fmap{P,T}
# test_ν      = deepcopy(Δ * t_impls)
# data_aprx_ν = deepcopy(Δ * t_impls)[:fx]
# test_t      = deepcopy(Δ * t_impls)
# data_aprx_t = deepcopy(Δ * t_impls)[:fx]
# ----- ipulse
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
test_ν      = deepcopy(Δ * t_impls)
data_aprx_ν = deepcopy(Δ * t_impls)[:fx]
test_t      = deepcopy(Δ * t_impls)
data_aprx_t = deepcopy(Δ * t_impls)[:fx]
# ---- simulation 
# test_ν      = deepcopy(Δ * sim_ν)
# data_aprx_ν = deepcopy( Δ * sim_ν )[:fx]
# test_t      = deepcopy(Δ * sim_t)
# data_aprx_t = deepcopy( Δ * sim_t )[:fx]

nsteps = 500; ϵ = 1/500
for n = nsteps:-1:1
    test_t -= ϵ * (logΔ²Σt * test_t)
    test_ν -= ϵ * (logΔ²Σν * test_ν)
    data_aprx_t -= ϵ * (aprx_logΔ²Σt * data_aprx_t)
    data_aprx_ν -= ϵ * (aprx_logΔ²Σν * data_aprx_ν)
end
data_aprx_t = data_aprx_t |> Fmap{P,T} |> x->deepcopy(Δ * x)[:fx]
data_aprx_ν = data_aprx_ν |> Fmap{P,T} |> x->deepcopy(Δ * x)[:fx]
data_t      = deepcopy(Δ * test_t)[:fx]
data_ν      = deepcopy(Δ * test_ν)[:fx]

figure()
subplot(2,1,1)
plot(data_t)
subplot(2,1,2)
plot(data_aprx_t, alpha=0.2)

figure()
subplot(2,1,1)
plot(data_ν)
subplot(2,1,2)
plot(data_aprx_ν, alpha=0.2)


#= 
semilogy(data_ν .|> abs)
semilogy(data_ν - data_aprx_ν .|> abs)
=#

#=
plot(data_ν)
plot(data_ν - data_aprx_ν)
=#




################################################
#  likelihood ratio test
################################################

data_aprx_ν = deepcopy(Δ * sim_ν)[:fx] 
data_aprx_t = deepcopy(Δ * sim_t)[:fx]
nsteps = 500; ϵ = 1/500
for n = nsteps:-1:1
    data_aprx_t -= ϵ * (aprx_logΔ²Σt * data_aprx_t)
    data_aprx_ν -= ϵ * (aprx_logΔ²Σν * data_aprx_ν)
end
data_aprx_t = data_aprx_t |> Fmap{P,T} |> x->deepcopy(Δ * x)[:fx]
data_aprx_ν = data_aprx_ν |> Fmap{P,T} |> x->deepcopy(Δ * x)[:fx]

# log det
aprx_logdet_ν = trace(aprx_logΔ²Σν) -  (2sum(log.(Laplacek.^2)) - log.(Laplacek.^2)[1] - log.(Laplacek.^2)[end])  
aprx_logdet_t = trace(aprx_logΔ²Σt) -  (2sum(log.(Laplacek.^2)) - log.(Laplacek.^2)[1] - log.(Laplacek.^2)[end])  
logdet_ν = (2sum(log.(Maternk)) - log.(Maternk)[1] - log.(Maternk)[end])  
logdet_t = (2sum(log.(Triangk)) - log.(Triangk)[1] - log.(Triangk)[end])  

# data term
aprx_ℓℓ_ν = - 0.5 * dot(sim_ν[:fx],  data_aprx_ν)
aprx_ℓℓ_t = - 0.5 * dot(sim_t[:fx],  data_aprx_t)
ℓℓ_ν      = - 0.5 * dot(sim_ν[:fx],  (Σν \ sim_ν)[:fx]) 
ℓℓ_t      = - 0.5 * dot(sim_t[:fx],  (Σt \ sim_t)[:fx])
#=
plot(data_aprx_ν)
plot((Σν \ sim_ν)[:fx])
plot(abs.((Σν \ sim_ν)[:fx]) ./ abs.(data_aprx_ν))
=#

# full likelihood
aprx_logPν = aprx_ℓℓ_ν  - 0.5 * aprx_logdet_ν
aprx_logPt = aprx_ℓℓ_t - 0.5 * aprx_logdet_t
logPν = - 0.5 * dot(sim_ν, Σν \ sim_ν) - 0.5 * logdet_ν
logPt = - 0.5 * dot(sim_t, Σt \ sim_t) - 0.5 * logdet_t


# Difference ... 
2*(aprx_logPν - logPν) / sqrt(2*nside)
2*(aprx_logPt - logPt) / sqrt(2*nside)






