
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
    Ωx     = Δx^2
    nyq    = 2π / (2Δx)
    k_side = ifftshift(-nside÷2:(nside-1)÷2) * Δk
    x_side = ifftshift(-nside÷2:(nside-1)÷2) * Δx
    k      = [reshape(k_side, 1, nside), reshape(k_side[1:nside÷2+1], nside÷2+1, 1)]
    x      = [reshape(x_side, 1, nside), reshape(x_side, nside, 1)]
    ϕk     = atan.(k[2],k[1])
    FFT    =  T(Ωx / (2π)) * plan_rfft(Matrix{T}(undef,nside,nside); flags=FFTW.ESTIMATE)
    #--- force the real hermitian symmitry for sin2ϕk ()
    sin2ϕk, cos2ϕk = sin.(2 .* ϕk), cos.(2 .* ϕk)
    if iseven(nside)
        sin2ϕk[1, end:-1:(nside÷2+2)] .= sin2ϕk[1, 2:nside÷2]
        sin2ϕk[end, end:-1:(nside÷2+2)] .= sin2ϕk[end, 2:nside÷2]
    end
    # ---------
    r𝔽{P,T,typeof(FFT)}(Δx, Δk, Ωk, Ωx, period, nyq, k, x, sin2ϕk, cos2ϕk, FFT)
end

# default T == Float64
r𝔽(::Type{P}) where P<:Flat = r𝔽(P, Float64)


# forward transform for scalar fields
function *(g::r𝔽{P,T}, tx::Matrix)::Matrix{Complex{T}} where {T<:Real, P<:Flat}
    g.FFT * tx
end

# forward transform for S2 fields
function *(g::r𝔽{P,T}, qux::Tuple{Matrix{T},Matrix{T}})::Tuple{Matrix{Complex{T}},Matrix{Complex{T}}} where {T<:Real, P<:Flat}
    qx, ux = qux
    qk, uk = g * qx, g * ux
    ek = similar(qk)
    bk = similar(qk)
    @inbounds @simd for I in eachindex(qk)
        ek[I] =   qk[I] * g.cos2ϕk[I] + uk[I] * g.sin2ϕk[I]
        bk[I] = - qk[I] * g.sin2ϕk[I] + uk[I] * g.cos2ϕk[I]
    end
    return (ek, bk)
end


# inverse transform for scalar fields
function \(g::r𝔽{P,T}, tk::Matrix)::Matrix{T} where {T<:Real, P<:Flat}
    g.FFT \ tk
end

# inverse transform for S2 fields
function \(g::r𝔽{P,T}, ebk::Tuple{Matrix{Complex{T}},Matrix{Complex{T}}})::Tuple{Matrix{T},Matrix{T}} where {T<:Real, P<:Flat}
    ek, bk = ebk
    qk = similar(ek)
    uk = similar(bk)
    @inbounds @simd for I in eachindex(ek)
        qk[I] =   ek[I] * g.cos2ϕk[I] - bk[I] * g.sin2ϕk[I]
        uk[I] =   ek[I] * g.sin2ϕk[I] + bk[I] * g.cos2ϕk[I]
    end
    qx, ux = g \ qk, g \ uk
    return (qx, ux)
end 


# allow the types to operate
(*)(::Type{r𝔽{P,T}}, x) where P<:Flat where T = r𝔽(P,T) * x
(\)(::Type{r𝔽{P,T}}, x) where P<:Flat where T = r𝔽(P,T) \ x
(*)(::Type{r𝔽{P}}, x)   where P<:Flat = r𝔽(P) * x
(\)(::Type{r𝔽{P}}, x)   where P<:Flat = r𝔽(P) \ x







