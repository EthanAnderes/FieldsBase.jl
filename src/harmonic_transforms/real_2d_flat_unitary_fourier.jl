############################################
#  unitary real FFT in 2 dimension
# fk   = sum_x f(x) exp(-2πik⋅x)/n^(d/2)
# f(x) = sum_k fk exp(2πik⋅x)/n^(d/2)
#############################################


struct r𝕌𝔽2{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
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
@generated function r𝕌𝔽2(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
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
    
    #FFT    =  (nside^(-dm/2)) * plan_rfft(Array{T}(undef, nside,nside); flags=FFTW.PATIENT, timelimit=45)
    #r𝕌𝔽2{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, FFT)

    X   = zeros(T,nside,nside) 
    mlt = T(nside^(-dm/2))
    code_out = quote 
        FFT    = $mlt * plan_rfft($X; flags=FFTW.MEASURE, timelimit=30)
        r𝕌𝔽2{$P,$T,typeof(FFT)}($Δx, $Δk, $Ωk, $Ωx, $period, $nyq, $k, $x, FFT)
    end
    return code_out
end


r𝕌𝔽2(::Type{P}) where P<:Flat = r𝕌𝔽2(P,Float64)

(*)(::Type{r𝕌𝔽2{P,T}}, x) where P<:Pix where T = r𝕌𝔽2(P,T).FFT * x
(\)(::Type{r𝕌𝔽2{P,T}}, x) where P<:Pix where T = r𝕌𝔽2(P,T).FFT \ x

(*)(::Type{r𝕌𝔽2{P}}, x)   where P<:Pix = r𝕌𝔽2(P,Float64).FFT * x
(\)(::Type{r𝕌𝔽2{P}}, x)   where P<:Pix = r𝕌𝔽2(P,Float64).FFT \ x

(*)(g::r𝕌𝔽2{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::r𝕌𝔽2{P,T}, x) where P<:Pix where T = g.FFT \ x




