# Some transform operators the user can choose from

fourier_transform(::Type{Any}) = error("define fourier_transform")

############################################
#  real FFT
#############################################

#  real FFT
struct r𝔽{P<:Flat,T<:Real,F}
    Δx::T
    Δk::T
    Ωk::T
    Ωpix::T
    period::T
    nyq::T
    k::Vector{Matrix{T}}
    x::Vector{Matrix{T}}
    sin2ϕk::Matrix{T}
    cos2ϕk::Matrix{T}
    FFT::F
end

# I want generated version of these....
# NOTE: these need to be defined for each fourier transform ...
function ebk_to_quk(ek, bk, ::Type{r𝔽{P,T}}) where {P<:Pix, T<:Real}
    g  = r𝔽(P,T)
    # ---------------
    @inbounds qk = .- ek .* g.cos2ϕk .+ bk .* g.sin2ϕk
    @inbounds uk = .- ek .* g.sin2ϕk .- bk .* g.cos2ϕk
    return qk, uk
    # -------------------
    # rw, cl = size(ek)
    # qk = Array{Complex{T},2}(rw,cl)
    # uk = Array{Complex{T},2}(rw,cl)
    # @inbounds for j = 1:cl
    #     @simd for i = 1:rw
    #         ϕk = atan2(g.k[2][i],g.k[1][j])
    #         cos2ϕk = cos(2ϕk)
    #         sin2ϕk = sin(2ϕk)
    #         qk[i,j] = - ek[i,j] * cos2ϕk + bk[i,j] * sin2ϕk
    #         uk[i,j] = - ek[i,j] * sin2ϕk - bk[i,j] * cos2ϕk
    #     end
    # end
    # return qk, uk
end
function quk_to_ebk(qk, uk, ::Type{r𝔽{P,T}}) where {P<:Pix, T<:Real}
    g  = r𝔽(P,T)
    # ----------------
    @inbounds ek = .- qk .* g.cos2ϕk .- uk .* g.sin2ϕk
    @inbounds bk =    qk .* g.sin2ϕk .- uk .* g.cos2ϕk
    return ek, bk
    # ----------------
    # rw, cl = size(qk)
    # ek = Array{Complex{T},2}(rw,cl)
    # bk = Array{Complex{T},2}(rw,cl)
    # @inbounds for j = 1:cl
    #     @simd for i = 1:rw
    #         ϕk = atan2(g.k[2][i],g.k[1][j])
    #         cos2ϕk = cos(2ϕk)
    #         sin2ϕk = sin(2ϕk)
    #         ek[i,j] = - qk[i,j] * cos2ϕk - uk[i,j] * sin2ϕk
    #         bk[i,j] =   qk[i,j] * sin2ϕk - uk[i,j] * cos2ϕk
    #     end
    # end
    # return ek, bk
end


# real FFT generated function constructor
@generated function r𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωpix   = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side[1:nside÷2+1], nside÷2+1, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    ϕk     = atan2.(k[2],k[1])
    FFT    =  Ωpix / (2π) * plan_rfft(rand(T,nside,nside); flags=PATIENT, timelimit=4)
    FFT    =  Ωpix / (2π) * plan_rfft(rand(T,nside,nside); flags=PATIENT, timelimit=4)
    r𝔽{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk), FFT)
end


r𝔽(::Type{P}) where P<:Flat = r𝔽(P,Float64)

(*)(::Type{r𝔽{P,T}}, x) where P<:Pix where T = r𝔽(P,T).FFT * x
(\)(::Type{r𝔽{P,T}}, x) where P<:Pix where T = r𝔽(P,T).FFT \ x

(*)(::Type{r𝔽{P}}, x)   where P<:Pix = r𝔽(P,Float64).FFT * x
(\)(::Type{r𝔽{P}}, x)   where P<:Pix = r𝔽(P,Float64).FFT \ x

(*)(g::r𝔽{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝔽{P,T}, x) where P<:Pix where T = g.FFT \ x




############################################
#  complex FFT
#############################################

#  complex FFT
struct 𝔽{P<:Flat,T<:Real,F}
    Δx::T
    Δk::T
    Ωk::T
    Ωpix::T
    period::T
    nyq::T
    k::Vector{Matrix{T}}
    x::Vector{Matrix{T}}
    sin2ϕk::Matrix{T}
    cos2ϕk::Matrix{T}
    FFT::F
end


function ebk_to_quk(ek, bk, ::Type{𝔽{P,T}}) where {P<:Pix, T<:Real}
    g  = r𝔽(P,T)
    qk = .- ek .* g.cos2ϕk .+ bk .* g.sin2ϕk
    uk = .- ek .* g.sin2ϕk .- bk .* g.cos2ϕk
    return qk, uk
end
function quk_to_ebk(qk, uk, ::Type{𝔽{P,T}}) where {P<:Pix, T<:Real}
    g  = r𝔽(P,T)
    ek = .- qk .* g.cos2ϕk .- uk .* g.sin2ϕk
    bk =    qk .* g.sin2ϕk .- uk .* g.cos2ϕk
    return ek, bk
end


# complex FFT generated function constructor
@generated function 𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωpix   = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side[1:nside÷2+1], nside÷2+1, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    ϕk     = atan2.(k[2],k[1])
    FFT    =  Ωpix / (2π) * plan_fft(rand(Complex{T},nside,nside); flags=PATIENT, timelimit=4)
    𝔽{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk), FFT)
end

𝔽(::Type{P}) where P<:Flat  = 𝔽(P,Float64)

(*)(::Type{𝔽{P,T}}, x) where P<:Pix where T = 𝔽(P,T).FFT * x
(\)(::Type{𝔽{P,T}}, x) where P<:Pix where T = 𝔽(P,T).FFT \ x

(*)(::Type{𝔽{P}}, x)   where P<:Pix = 𝔽(P,Float64).FFT * x
(\)(::Type{𝔽{P}}, x)   where P<:Pix = 𝔽(P,Float64).FFT \ x

(*)(g::𝔽{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::𝔽{P,T}, x) where P<:Pix where T = g.FFT \ x



############################################
#  Healpix transform
#############################################

 struct ℍ{P<:Healpix,T<:Real}
    Ωpix::T
    lmax::T
    l::Matrix{T}
    m::Matrix{T}
    φ::Matrix{T}  # azmuth
    Θ::Matrix{T}  # polar
end

# function ebk_to_quk(ek, bk, ::Type{ℍ{P,T}}) where {P<:Pix, T<:Real}
#     g  = r𝔽(P,T)
#     qk = .- ek .* g.cos2ϕk .+ bk .* g.sin2ϕk
#     uk = .- ek .* g.sin2ϕk .- bk .* g.cos2ϕk
#     return qk, uk
# end
# function quk_to_ebk(qk, uk, ::Type{ℍ{P,T}}) where {P<:Pix, T<:Real}
#     g  = r𝔽(P,T)
#     ek = .- qk .* g.cos2ϕk .- uk .* g.sin2ϕk
#     bk =    qk .* g.sin2ϕk .- uk .* g.cos2ϕk
#     return ek, bk
# end
