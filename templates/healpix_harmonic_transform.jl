
using FieldsBase
using PyCall

import Base: *, \

const UNSEEN = -1.6375e30


abstract type Healpix{nside,lmax} <: Pix end


struct rℍ{P<:Healpix,T<:Real} <: HarmonicTransform{P,T}
    nside::Int
    lmax::Int
    npix::Int
    Ωx::T
    l::Vector{Int}
    m::Vector{Int}
    θ::Vector{T}  # polar
    φ::Vector{T}  # azmuth
end


@generated function rℍ(::Type{P},::Type{T}) where {nside, lmax, T<:Real, P<:Healpix{nside,lmax}}
    hp    = pyimport("healpy") 
    θ, φ  = hp.pix2ang(nside, 0:(hp.nside2npix(nside)-1))
    Ωx    = abs2(hp.pixelfunc.nside2resol(nside, arcmin=false))
    npix  = hp.nside2npix(nside)
    l, m  = hp.sphtfunc.Alm.getlm(lmax)
    rℍ{P,T}(nside, lmax, npix, Ωx, l, m, θ, φ)
end

# default T == Float64
# rℍ(::Type{P}) where P<:Flat = rℍ(P, Float64)


# forward transform for scalar fields
function *(g::rℍ{P,T}, tx::Vector{T})::Vector{Complex{T}} where {nside, lmax, T<:Real, P<:Healpix{nside,lmax}}
    hps  = pyimport("healpy.sphtfunc") 
    hps.map2alm(tx, lmax=lmax, iter=0, pol=false) 
end

# forward transform for spin 2 fields
function *(g::rℍ{P,T}, qux::Tuple{Vector{T},Vector{T}})::Tuple{Vector{Complex{T}},Vector{Complex{T}}} where {nside, lmax, T<:Real, P<:Healpix{nside,lmax}}
    hps  = pyimport("healpy.sphtfunc") 
    ebℓm_vec = hps._sphtools.map2alm_spin_healpy(qux, 2, lmax=lmax) 
    return (ebℓm_vec[1], ebℓm_vec[2])
end

# inverse transform for scalar fields
function \(g::rℍ{P,T}, tℓm::Vector{Complex{T}})::Vector{T} where {nside, lmax, T<:Real, P<:Healpix{nside,lmax}}
    hps  = pyimport("healpy.sphtfunc") 
    hps.alm2map(tℓm, nside, lmax=lmax, pol=false, verbose=false)
end

# inverse transform for spin 2 fields
function \(g::rℍ{P,T}, ebℓm::Tuple{Vector{Complex{T}},Vector{Complex{T}}})::Tuple{Vector{T},Vector{T}} where {nside, lmax, T<:Real, P<:Healpix{nside,lmax}}
    hps  = pyimport("healpy.sphtfunc") 
    qux_vec = hps._sphtools.alm2map_spin_healpy(ebℓm, nside, 2, lmax)
    return (qux_vec[1], qux_vec[2])
end 



# allow the types to operate
(*)(::Type{rℍ{P,T}}, x) where {T<:Real, P<:Healpix} = rℍ(P,T) * x
(\)(::Type{rℍ{P,T}}, x) where {T<:Real, P<:Healpix} = rℍ(P,T) \ x
(*)(::Type{rℍ{P}}, x)   where P<:Healpix = rℍ(P) * x
(\)(::Type{rℍ{P}}, x)   where P<:Healpix = rℍ(P) \ x

