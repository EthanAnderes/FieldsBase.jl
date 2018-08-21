
############################################
#  real FFT
#############################################


#  1-d real FFT
struct r𝔽1{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωx::T
    period::T
    nyq::T
    k::Vector{T}
    x::Vector{T}
    FFT::F
end

@generated function r𝔽1(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk
    Ωx     = Δx
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = k_side[1:nside÷2+1]
    x      = x_side
    dm     = 1 #<-- dimension
    FFT    =  (Ωx * ((2π) ^ (-dm/2))) * plan_rfft(Array{T}(undef, nside))
    r𝔽1{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, FFT)
end

r𝔽1(::Type{P}) where P<:Flat = r𝔽1(P,Float64)

(*)(::Type{r𝔽1{P,T}}, x) where P<:Pix where T = r𝔽1(P,T).FFT * x
(\)(::Type{r𝔽1{P,T}}, x) where P<:Pix where T = r𝔽1(P,T).FFT \ x

(*)(::Type{r𝔽1{P}}, x)   where P<:Pix = r𝔽1(P,Float64).FFT * x
(\)(::Type{r𝔽1{P}}, x)   where P<:Pix = r𝔽1(P,Float64).FFT \ x

(*)(g::r𝔽1{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝔽1{P,T}, x) where P<:Pix where T = g.FFT \ x


