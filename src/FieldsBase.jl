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
export data, HasQU, has_qu, IsMap, is_map, IsLenseBasis, is_lense_basis # used for traits

include("fourier_transforms.jl")
export fourier_transform, r𝔽, 𝔽, ℍ # <-- fourier and spherical transforms
#TODO make a general dimension r𝔽

include("field_ops.jl")
# export ...
#TODO

# include("linear_ops.jl")
# export ...
#TODO

# include("lense_transforms.jl")
# export ...
#TODO


end # end Module
