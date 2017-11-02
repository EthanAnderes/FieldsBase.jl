
using PyPlot 
include("/Users/ethananderes/Dropbox/FieldsBase/templates/t_flat_1dimension_unitary.jl")

nside = 2^12
Θpix  = 1/nside
dm    = 1
P     = Flat{Θpix,nside}
T     = Float64
g     = r𝕌𝔽1(P,T);

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


ρ = 0.05 * g.period
ν = 1.4

# -------- spectral densities ------
Maternk     = normalize_Λ(( 4ν/ρ^2 + abs2.(g.k) ) .^ (- ν - dm/2), nside, dm)
Triangk     = normalize_Λ(real.(g * map(x->trang(x,ρ), g.x)), nside, dm)
# plot( (1/nside^(dm/2)) * ( g \ Maternk)) # this is the auto-cov function

# ------ cov operators 
Σν = Maternk  |> Ffourier{P,T} |> 𝕃
Σt = Triangk  |> Ffourier{P,T} |> 𝕃
#= plot 
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
plot((Σν * t_impls)[:fx])
=#

#  ------ log operators ------ 
logΣν = log.(Maternk)  |> Ffourier{P,T} |> 𝕃
logΣt = log.(Triangk)  |> Ffourier{P,T} |> 𝕃
#= plot 
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
plot((logΣν * t_impls)[:fx])
=#

# ----- simulate ------
sqrtΣν = sqrt.(Maternk)  |> Ffourier{P,T} |> 𝕃 # note: no factor (2π)^(-dm/2)
sqrtΣt = sqrt.(Triangk)  |> Ffourier{P,T} |> 𝕃
sim_ν = sqrtΣν * Fmap{P,T}(white_noise(g))
sim_t = sqrtΣt * Fmap{P,T}(white_noise(g))
#=
plot(sim_ν[:fx][1:1000])
plot(sim_t[:fx][1:1000])
=#

########################################
#  sparse approx to log
########################################

# --- truncation
# sprs_sz = 50 # approx each col with this number of non-zeros
# splogΣν = spzeros(nside, nside)
# splogΣt = spzeros(nside, nside)
# for ind = 1:nside
#     t_impls = (t=zeros(T,nside); t[ind] = 1; t) |> Fmap{P,T}
#     logΣν_col = (logΣν * t_impls)[:fx]
#     logΣt_col = (logΣt * t_impls)[:fx]
#     pν = sortperm(abs2.(logΣν_col), rev=true)
#     pt = sortperm(abs2.(logΣt_col), rev=true)
#     tmp_approx_ν = logΣν_col[pν[1:sprs_sz]]
#     tmp_approx_t = logΣt_col[pt[1:sprs_sz]]
#     splogΣν[pν[1:sprs_sz], ind] = tmp_approx_ν
#     splogΣt[pt[1:sprs_sz], ind] = tmp_approx_t
# end


# -------  local cov matrix log
window    = -50:50
sprs_sz   = length(window)
get_sparse_logΣ = function ()
	intrir_Σν = zeros(sprs_sz,sprs_sz)
	intrir_Σt = zeros(sprs_sz,sprs_sz)
	mid_pnt = nside÷2
	for k = 1:sprs_sz
		t_impls = (t=zeros(T,nside); t[mid_pnt + window[k]] = 1; t) |> Fmap{P,T}
		intrir_Σν[:,k] =  (Σν * t_impls)[:fx][mid_pnt .+ window]
		intrir_Σt[:,k] =  (Σt * t_impls)[:fx][mid_pnt .+ window]
	end
	intrir_logΣν = real.(logm(intrir_Σν))
	intrir_logΣt = real.(logm(intrir_Σt)) 
	
	splogΣν = spzeros(nside, nside)
	splogΣt = spzeros(nside, nside)
	for ind = 1:nside
		if 	all(1 .<= ind .+ window .<= nside)
			splogΣν[ind + window, ind] = intrir_logΣν[:,(sprs_sz-1)÷2+1]
			splogΣt[ind + window, ind] = intrir_logΣt[:,(sprs_sz-1)÷2+1]
		else 
			nearst_ind = sortperm(abs.(ind .- (1:nside)))[1:sprs_sz]
			bdry_Σν = zeros(sprs_sz, sprs_sz)
			bdry_Σt = zeros(sprs_sz, sprs_sz)
			for k = 1:sprs_sz
				t_impls = (t=zeros(T,nside); t[nearst_ind[k]] = 1; t) |> Fmap{P,T}
				bdry_Σν[:,k] =  (Σν * t_impls)[:fx][nearst_ind]
				bdry_Σt[:,k] =  (Σt * t_impls)[:fx][nearst_ind]
			end
			bdry_logΣν = real.(logm(bdry_Σν))
			bdry_logΣt = real.(logm(bdry_Σt)) 
			splogΣν[nearst_ind, ind] = bdry_logΣν[:,1]
			splogΣt[nearst_ind, ind] = bdry_logΣt[:,1]
		end
	end
	return splogΣν, splogΣt
end
@time splogΣν, splogΣt = get_sparse_logΣ()


# ----- plot --
#=
ind = 1000
t_impls = (t=zeros(T,nside); t[ind] = 1; t) |> Fmap{P,T}
logΣν_col = (logΣν * t_impls)[:fx] 
plot(splogΣν[:,ind])
plot(logΣν_col)
plot(logΣν_col - splogΣν[:,ind])
=#


#= plot 

  | [^AH12]
  |
  |  Awad H. Al-Mohy and Nicholas J. Higham, "Improved inverse scaling and squaring
  |  algorithms for the matrix logarithm", SIAM Journal on Scientific Computing, 34(4),
  |  2012, C153-C169. doi:10.1137/110852553

  | [^AHR13]
  |
  |  Awad H. Al-Mohy, Nicholas J. Higham and Samuel D. Relton, "Computing the Fréchet
  |  derivative of the matrix logarithm and estimating the condition number", SIAM
  |  Journal on Scientific Computing, 35(4), 2013, C394-C410. doi:10.1137/120885991

=#





########################################
#  test the forward map
########################################

# ----- edge 
t_impls  = (t=zeros(T,nside); t[2000:2000+5*sprs_sz] = 1; t) |> Fmap{P,T}
test_ν   = deepcopy(t_impls)
data_spν = deepcopy(t_impls)[:fx] 
test_t   = deepcopy(t_impls)
data_spt = deepcopy(t_impls)[:fx] 
# --- ipulse
# t_impls    = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
# test_ν     = deepcopy(t_impls)
# data_spν   = deepcopy(t_impls)[:fx] 
# test_t     = deepcopy(t_impls)
# data_spt   = deepcopy(t_impls)[:fx] 

nsteps     = 1000; ϵ = 1/1000
for n      = 0:nsteps-1
	test_t     += ϵ * (logΣt * test_t)
	test_ν     += ϵ * (logΣν * test_ν)
	data_spt   += ϵ * (splogΣt * data_spt)
	data_spν   += ϵ * (splogΣν * data_spν)
end
data_spt   = data_spt
data_spν   = data_spν
data_t     = deepcopy(test_t)[:fx]
data_ν     = deepcopy(test_ν)[:fx] 

figure()
subplot(2,1,1)
plot(data_t)
subplot(2,1,2)
plot(data_spt)

figure()
subplot(2,1,1)
plot(data_ν)
subplot(2,1,2)
plot(data_spν)

#=
plot(data_ν)
plot(data_ν-data_spν)
=#

################################################
#  test the inverse map
################################################

# ----- edge 
# t_impls  = (t=zeros(T,nside); t[2000:2000+2*sprs_sz] = 1; t) |> Fmap{P,T}
# test_ν   = deepcopy(t_impls)
# data_spν = deepcopy(t_impls)[:fx]
# test_t   = deepcopy(t_impls)
# data_spt = deepcopy(t_impls)[:fx]
# ----- ipulse
t_impls  = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
test_ν   = deepcopy(t_impls)
data_spν = deepcopy(t_impls)[:fx]
test_t   = deepcopy(t_impls)
data_spt = deepcopy(t_impls)[:fx]

nsteps = 1000; ϵ = 1/1000
for n = nsteps:-1:1
	test_t   -= ϵ * (logΣt * test_t)
	test_ν   -= ϵ * (logΣν * test_ν)
	data_spt -= ϵ * (splogΣt * data_spt)
	data_spν -= ϵ * (splogΣν * data_spν)
end
data_spt = data_spt
data_spν = data_spν
data_t   = deepcopy(test_t)[:fx]
data_ν   = deepcopy(test_ν)[:fx]

figure()
subplot(2,1,1)
plot(data_t)
subplot(2,1,2)
plot(data_spt)

figure()
subplot(2,1,1)
plot(data_ν)
subplot(2,1,2)
plot(data_spν)


#= 
semilogy(data_ν .|> abs)
semilogy(data_ν - data_spν .|> abs)
=#

#=
plot(data_ν)
plot(data_ν - data_spν)
=#




################################################
#  likelihood ratio test
################################################

data_spν = deepcopy(sim_ν)[:fx] 
data_spt = deepcopy(sim_t)[:fx]
nsteps = 1000; ϵ = 1/1000
for n = nsteps:-1:1
    data_spt -= ϵ * (splogΣt * data_spt)
    data_spν -= ϵ * (splogΣν * data_spν)
end
data_spt = data_spt
data_spν = data_spν

# log det
splogdet_ν = trace(splogΣν) #-  (2sum(log.(Laplacek.^2)) - log.(Laplacek.^2)[1] - log.(Laplacek.^2)[end])  
splogdet_t = trace(splogΣt) #-  (2sum(log.(Laplacek.^2)) - log.(Laplacek.^2)[1] - log.(Laplacek.^2)[end])  
logdet_ν = (2sum(log.(Maternk)) - log.(Maternk)[1] - log.(Maternk)[end])  
logdet_t = (2sum(log.(Triangk)) - log.(Triangk)[1] - log.(Triangk)[end])  

# data term
spℓℓ_ν = - 0.5 * dot(sim_ν[:fx],  data_spν)
spℓℓ_t = - 0.5 * dot(sim_t[:fx],  data_spt)
ℓℓ_ν      = - 0.5 * dot(sim_ν[:fx],  (Σν \ sim_ν)[:fx]) 
ℓℓ_t      = - 0.5 * dot(sim_t[:fx],  (Σt \ sim_t)[:fx])
#=
plot(data_spν)
plot((Σν \ sim_ν)[:fx])
plot(abs.((Σν \ sim_ν)[:fx]) ./ abs.(data_spν))
=#

# full likelihood
splogPν = spℓℓ_ν  - 0.5 * splogdet_ν
splogPt = spℓℓ_t - 0.5 * splogdet_t
logPν = - 0.5 * dot(sim_ν, Σν \ sim_ν) - 0.5 * logdet_ν
logPt = - 0.5 * dot(sim_t, Σt \ sim_t) - 0.5 * logdet_t


# Difference ... 
2*(splogPν - logPν) / sqrt(2*nside)
2*(splogPt - logPt) / sqrt(2*nside)















###################################################
# old 
######################################


# ---------- whitening operator
# TODO: Try writing Σ = Y⁻¹ * Y * Σ * Y⁻¹ * Y  =  Y⁻¹ * exp(log(Y * Σ * Y⁻¹)) * Y
# TODO: Try writing Σ⁻¹ = Δ * (Σ * Δ)⁻¹ 
#w_pwr     = 0
#w_rng     = 3*g.Δx
#Laplacek  = normalize_Λ(( 4*w_pwr/(w_rng)^2 + abs2.(g.k) ) .^ w_pwr, nside, dm)
#Laplacek   = ( 4*w_pwr/(w_rng)^2 + abs2.(g.k) ) .^ w_pwr
# ave_Δ²k    = (2sum(Laplacek.^2) - Laplacek[1]^2 - Laplacek[end]^2) / (nside^dm)
# Laplacek ./= sqrt(ave_Δ²k)
# ave_Δ²k    = (2sum(Maternk .* Laplacek .^ 2) - Maternk[1]*Laplacek[1]^2 - Maternk[end]*Laplacek[end]^2) / (nside^dm)
# Laplacek ./= sqrt(ave_Δ²k)

#Laplacek = ones(g.k)
#Δ         = Laplacek |> Ffourier{P,T} |> 𝕃
#= plot 
t_impls = (t=zeros(T,nside); t[1000] = 1; t) |> Fmap{P,T}
plot((Δ * t_impls)[:fx])
=#



# ---- pre-whiten -----
# Δ²Maternk  = Maternk .* Laplacek .^ 2
# Δ²Triangk  = Triangk .* Laplacek .^ 2
# Δ²Σν = Δ²Maternk  |> Ffourier{P,T} |> 𝕃
# Δ²Σt = Δ²Triangk  |> Ffourier{P,T} |> 𝕃
