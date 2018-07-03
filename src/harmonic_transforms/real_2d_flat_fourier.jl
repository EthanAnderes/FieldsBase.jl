
############################################
#  real FFT
#############################################

#  real FFT
struct r𝔽{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωx::T
    period::T
    nyq::T
    k::Vector{Matrix{T}}
    x::Vector{Matrix{T}}
    sin2ϕk::Matrix{T}
    cos2ϕk::Matrix{T}
    FFT::F
end


# real FFT generated function constructor
@generated function r𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωx   = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side[1:nside÷2+1], nside÷2+1, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    ϕk     = atan.(k[2],k[1])
    FFT    =  Ωx / (2π) * plan_rfft(Array{T}(undef,nside,nside); flags=FFTW.MEASURE); #flags=FFTW.PATIENT, timelimit=4)
    r𝔽{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, sin.(2 .* ϕk), cos.(2 .* ϕk), FFT)
end


r𝔽(::Type{P}) where P<:Flat = r𝔽(P,Float64)

(*)(::Type{r𝔽{P,T}}, x) where P<:Pix where T = r𝔽(P,T).FFT * x
(\)(::Type{r𝔽{P,T}}, x) where P<:Pix where T = r𝔽(P,T).FFT \ x

(*)(::Type{r𝔽{P}}, x)   where P<:Pix = r𝔽(P,Float64).FFT * x
(\)(::Type{r𝔽{P}}, x)   where P<:Pix = r𝔽(P,Float64).FFT \ x

(*)(g::r𝔽{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝔽{P,T}, x) where P<:Pix where T = g.FFT \ x





