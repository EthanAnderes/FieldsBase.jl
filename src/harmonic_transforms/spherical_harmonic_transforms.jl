

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

# function harmonic_eb_to_qu(elm, blm, g::ℍ{P,T}) where {P<:Pix, T<:Real}
#     # TODO add the spherical conversion
#     return qlm, ulm
# end
# function harmonic_qu_to_eb(qlm, ulm, g::ℍ{P,T}) where {P<:Pix, T<:Real}
#     # TODO add the spherical conversion
#     return elm, blm
# end
