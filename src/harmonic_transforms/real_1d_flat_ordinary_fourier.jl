############################################
#  ordinary real FFT in 1 dimension
# fk   = ∫ f(x) exp(-2πik⋅x) dx
# f(x) = sum_k fk exp(2πik⋅x)
#############################################


struct r𝕆𝔽1{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
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

@generated function r𝕆𝔽1(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    dm     = 1 #<-- dimension
    Δx     = Θpix
    period = Δx*nside
    Δk     = 1/period
    Ωk     = Δk^dm
    Ωx     = Δx^dm
    nyq    = 1 / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = k_side[1:nside÷2+1]
    x      = x_side
    #FFT    = (1/nside^dm) * plan_rfft(Array{T}(nside); flags=FFTW.PATIENT, timelimit=4)  # unitary normization
    FFT    = (1/nside^dm) * plan_rfft(Array{T}(undef, nside); flags=FFTW.MEASURE)  # unitary normization
    r𝕆𝔽1{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, FFT)
end


r𝕆𝔽1(::Type{P}) where P<:Flat = r𝕆𝔽1(P,Float64)

(*)(::Type{r𝕆𝔽1{P,T}}, x) where P<:Pix where T = r𝕆𝔽1(P,T).FFT * x
(\)(::Type{r𝕆𝔽1{P,T}}, x) where P<:Pix where T = r𝕆𝔽1(P,T).FFT \ x

(*)(::Type{r𝕆𝔽1{P}}, x)   where P<:Pix = r𝕆𝔽1(P,Float64).FFT * x
(\)(::Type{r𝕆𝔽1{P}}, x)   where P<:Pix = r𝕆𝔽1(P,Float64).FFT \ x

(*)(g::r𝕆𝔽1{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝕆𝔽1{P,T}, x) where P<:Pix where T = g.FFT \ x


