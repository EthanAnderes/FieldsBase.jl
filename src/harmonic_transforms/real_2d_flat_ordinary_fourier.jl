############################################
#  unitary real FFT in 2 dimension
# fk   = sum_x f(x) exp(-2πik⋅x)/n^(d/2)
# f(x) = sum_k fk exp(2πik⋅x)/n^(d/2)
#############################################


struct r𝕆𝔽2{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
    Δx::T
    Δk::T
    Ωk::T
    Ωx::T
    period::T
    nyq::T
    k::Vector{Matrix{T}}
    x::Vector{Matrix{T}}
    FFT::F
end


# real FFT generated function constructor
@generated function r𝕆𝔽2(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    dm     = 2 #<-- dimension
    Δx     = Θpix
    period = Δx*nside
    Δk     = 1/period
    Ωk     = Δk^dm
    Ωx     = Δx^dm
    nyq    = 1 / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side[1:nside÷2+1], nside÷2+1, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    FFT    = (1/nside^dm) * plan_rfft(Array{T}(undef, nside,nside); flags=FFTW.MEASURE) #; flags=FFTW.PATIENT, timelimit=4)
    r𝕆𝔽2{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, FFT)
end

r𝕆𝔽2(::Type{P}) where P<:Flat = r𝕆𝔽2(P,Float64)

(*)(::Type{r𝕆𝔽2{P,T}}, x) where P<:Pix where T = r𝕆𝔽2(P,T).FFT * x
(\)(::Type{r𝕆𝔽2{P,T}}, x) where P<:Pix where T = r𝕆𝔽2(P,T).FFT \ x

(*)(::Type{r𝕆𝔽2{P}}, x)   where P<:Pix = r𝕆𝔽2(P,Float64).FFT * x
(\)(::Type{r𝕆𝔽2{P}}, x)   where P<:Pix = r𝕆𝔽2(P,Float64).FFT \ x

(*)(g::r𝕆𝔽2{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝕆𝔽2{P,T}, x) where P<:Pix where T = g.FFT \ x




