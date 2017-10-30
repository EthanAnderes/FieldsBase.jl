module FieldsBase

if VERSION >= v"0.7.0-DEV.1"
    using FFTW: ifftshift, plan_rfft, plan_fft, PATIENT
else
    using Base.DFT.FFTW: PATIENT
end
import Base: +, -, *, ^, \, getindex, promote_rule, convert, show, dot, inv


####################################
abstract type Pix end
abstract type Flat{Θpix,nside} <: Pix end
abstract type Healpix{nside}   <: Pix end
abstract type Spin end
abstract type S0 <: Spin end
abstract type S2 <: Spin end
abstract type S02 <: Spin end
abstract type Field{P<:Pix, S<:Spin} end 
#TODO: Might be best to make this Field{P<:Pix, T<:Real, S<:Spin} 
# I think it would be more useful ... 
# test it out for complex fields in Pixel space and Healpix fields
abstract type HarmonicTransform{P<:Pix, T<:Real} end
export Pix, Flat, Healpix, Spin, S0, S2, S02, Field, HarmonicTransform


####################################
include("convert_promote.jl")
export HasQU, has_qu, IsMap, is_map, IsLenseBasis, is_lense_basis # used for traits
export harmonic_transform

####################################
include("harmonic_transforms/fourier_transforms.jl")
export r𝔽, 𝔽

####################################
include("harmonic_transforms/spherical_harmonic_transforms.jl")
export ℍ

####################################
include("field_ops.jl")

####################################
include("linear_ops.jl")
export LinOp, DiagOp, 𝕃

####################################
#TODO include("lense_transforms.jl")

####################################
include("util.jl")
export data, squash, white_noise



#TODO: code up squash that can be computed directly on the GPU.
#TODO: dot for Pix <: Healpix in field_ops.jl
#TODO: can we use splatting to get general spin SN fields for N>2
#TODO: get a spherical example up and running
#TODO: can we get rid of the "Method definition overwritten" when overloading harmonic_eb_to_qu, etc...
#TODO: can we add a basic pixel space lensing algo which is opt-in for users (and works on a GPU).

end # end Module
