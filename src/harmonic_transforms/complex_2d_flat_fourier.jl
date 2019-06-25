



############################################
#  complex FFT
#############################################

#  complex FFT
struct 𝔽{P<:Flat,T<:Real,F} <: HarmonicTransform{P,T}
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


# complex FFT generated function constructor
@generated function 𝔽(::Type{P},::Type{T}) where T<:Real where P<:Flat{Θpix, nside}  where {Θpix, nside}
    Δx     = deg2rad(Θpix/60)
    period = Δx*nside
    Δk     = 2π/period
    Ωk     = Δk^2
    Ωx     = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side, nside, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    ϕk     = atan.(k[2],k[1])
    FFT    =  Ωx / (2π) * plan_fft(Array{Complex{T}}(undef, nside,nside); flags=FFTW.MEASURE) #; flags=FFTW.PATIENT, timelimit=4)
    #--- force the real hermitian symmitry for sin2ϕk ()
    sin2ϕk = sin.(2 .* ϕk)
    if iseven(nside)
        sin2ϕk[1, end:-1:(nside÷2+2)] .= sin2ϕk[1, 2:nside÷2]
        sin2ϕk[nside÷2+1, end:-1:(nside÷2+2)] .= sin2ϕk[nside÷2+1, 2:nside÷2] # needs testing
    end
    # ---------
    𝔽{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, sin2ϕk, cos.(2 .* ϕk), FFT)
end


# default T == Float64
𝔽(::Type{P}) where P<:Flat = 𝔽(P, Float64)

# forward transform for scalar fields
function *(g::𝔽{P,T}, tx::Matrix)::Matrix{Complex{T}} where {T<:Real, P<:Flat}
    g.FFT * tx
end

# inverse transform for scalar fields
function \(g::𝔽{P,T}, tk::Matrix{Complex{T}})::Matrix{Complex{T}} where {T<:Real, P<:Flat}
    g.FFT \ tk
end

(*)(::Type{𝔽{P,T}}, x) where P<:Pix where T = 𝔽(P,T).FFT * x
(\)(::Type{𝔽{P,T}}, x) where P<:Pix where T = 𝔽(P,T).FFT \ x

(*)(::Type{𝔽{P}}, x)   where P<:Pix = 𝔽(P,Float64).FFT * x
(\)(::Type{𝔽{P}}, x)   where P<:Pix = 𝔽(P,Float64).FFT \ x

(*)(g::𝔽{P,T}, x) where P<:Pix where T = g.FFT * x
(\)(g::𝔽{P,T}, x) where P<:Pix where T = g.FFT \ x




# allow the types to operate
(*)(::Type{𝔽{P,T}}, x) where P<:Flat where T = 𝔽(P,T) * x
(\)(::Type{𝔽{P,T}}, x) where P<:Flat where T = 𝔽(P,T) \ x
(*)(::Type{𝔽{P}}, x)   where P<:Flat = 𝔽(P) * x
(\)(::Type{𝔽{P}}, x)   where P<:Flat = 𝔽(P) \ x



