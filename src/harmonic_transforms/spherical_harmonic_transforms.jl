

############################################
#  Healpix transform
#############################################
# real transform
 struct rℍ{P<:Healpix,T<:Real} <: HarmonicTransform{P,T}
    Ωx::T
    lmax::T
    l::Matrix{T}
    m::Matrix{T}
    φ::Matrix{T}  # azmuth
    Θ::Matrix{T}  # polar
end
