############################################
#  unitary real FFT in 1 dimension
#############################################


struct r𝕌𝔽1{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωpix::T
    period::T
    nyq::T
    k::Vector{T}
    x::Vector{T}
    FFT::F
end

@generated function r𝕌𝔽1(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    dm     = 1 #<-- dimension
    Δx     = Θpix
    period = Δx*nside
    Δk     = 1/period
    Ωk     = Δk^dm
    Ωpix   = Δx^dm
    nyq    = 1 / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = k_side[1:nside÷2+1]
    x      = x_side
    FFT    = (nside^(-dm/2)) * plan_rfft(rand(T,nside); flags=PATIENT, timelimit=4)  # unitary normization
    r𝕌𝔽1{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωpix, period, nyq, k, x, FFT)
end


r𝕌𝔽1(::Type{P}) where P<:Flat = r𝕌𝔽1(P,Float64)

(*)(::Type{r𝕌𝔽1{P,T}}, x) where P<:Pix where T = r𝕌𝔽1(P,T).FFT * x
(\)(::Type{r𝕌𝔽1{P,T}}, x) where P<:Pix where T = r𝕌𝔽1(P,T).FFT \ x

(*)(::Type{r𝕌𝔽1{P}}, x)   where P<:Pix = r𝕌𝔽1(P,Float64).FFT * x
(\)(::Type{r𝕌𝔽1{P}}, x)   where P<:Pix = r𝕌𝔽1(P,Float64).FFT \ x

(*)(g::r𝕌𝔽1{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝕌𝔽1{P,T}, x) where P<:Pix where T = g.FFT \ x


