module FieldsBase

if VERSION >= v"0.7.0-DEV.1"
    using FFTW: ifftshift, plan_rfft, plan_fft, PATIENT
else
    using Base.DFT.FFTW: PATIENT
end
import Base: +, -, *, ^, \, getindex, promote_rule, convert, show, dot, inv

abstract type Pix end
abstract type Flat{Θpix,nside} <: Pix end
abstract type Healpix{nside}   <: Pix end
abstract type Spin end
abstract type S0 <: Spin end
abstract type S2 <: Spin end
abstract type S02 <: Spin end
abstract type Field{Px<:Pix, S<:Spin} end
export Pix, Flat, Healpix, Spin, S0, S2, S02, Field

include("convert_promote.jl")
export HasQU, has_qu, IsMap, is_map, IsLenseBasis, is_lense_basis # used for traits
export harmonic_transform

include("harmonic_transforms/fourier_transforms.jl")
export r𝔽, 𝔽

include("harmonic_transforms/spherical_harmonic_transforms.jl")
export ℍ # fourier and spherical transforms

include("field_ops.jl")

include("linear_ops.jl")
export LinOp, DiagOp, 𝕃

#TODO include("lense_transforms.jl")

include("util.jl")
export data, squash, white_noise


end # end Module
