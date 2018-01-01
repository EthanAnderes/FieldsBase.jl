__precompile__() #getting warnings with v0.6 but works fine for v0.7

module FieldsBase

if VERSION >= v"0.7.0-DEV.1"
     using FFTW
end
FFTW.set_num_threads(Base.Threads.nthreads())
BLAS.set_num_threads(Base.Threads.nthreads())

import Base: +, -, *, ^, \, getindex, promote_rule, convert, show, dot, inv

const source_path = Base.source_path()

# abstract grid geometry and the corresponding harmonic transforms
abstract type Pix end
abstract type Flat{Θpix,nside} <: Pix end
abstract type Healpix{nside}   <: Pix end
abstract type HarmonicTransform{P<:Pix, T<:Real} end

# abstract field types
abstract type Spin end
abstract type S0 <: Spin end
abstract type S2 <: Spin end
abstract type S02 <: Spin end
abstract type Field{P<:Pix, T<:Real, S<:Spin} end 

# exported abstract types
export Pix, Flat, Healpix, Spin, S0, S2, S02, Field, HarmonicTransform

# traits used for convert and promote
include("convert_promote.jl")
export harmonic_transform, HasQU, has_qu, IsMap, is_map, IsLenseBasis, is_lense_basis 

# Harmonic transforms
include("harmonic_transforms/complex_2d_flat_fourier.jl")
include("harmonic_transforms/real_1d_flat_unitary_fourier.jl")
include("harmonic_transforms/real_2d_flat_unitary_fourier.jl")
include("harmonic_transforms/real_1d_flat_ordinary_fourier.jl")
include("harmonic_transforms/real_2d_flat_ordinary_fourier.jl")
include("harmonic_transforms/spherical_harmonic_transforms.jl")
include("harmonic_transforms/real_2d_flat_fourier.jl")
export r𝔽, r𝕆𝔽1, r𝕆𝔽2, r𝕌𝔽1, r𝕌𝔽2, rℍ

# field operations
include("field_ops.jl")

# linear ops
include("linear_ops.jl")
export LinOp, DiagOp, 𝕃

# misc 
include("util.jl")
export data, squash, white_noise


#TODO: dot for Pix <: Healpix in field_ops.jl
#TODO: get a spherical example up and running

end # end Module
