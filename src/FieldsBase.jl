module FieldsBase

using  FFTW
using  LinearAlgebra
import LinearAlgebra: dot
import Base: +, -, *, ^, \, getindex, promote_rule, convert, show, inv
export dot

const module_dir  = joinpath(@__DIR__, "..") |> normpath

# abstract grid geometry and the corresponding harmonic transforms
abstract type Pix end
abstract type Flat{Θpix,nside} <: Pix end
abstract type HarmonicTransform{P<:Pix, T<:Real} end

# abstract field types
abstract type Spin end
abstract type S0 <: Spin end
abstract type S2 <: Spin end
abstract type S02 <: Spin end
abstract type Field{P<:Pix, T<:Real, S<:Spin} end

# exported abstract types
export Pix, Flat, Spin, S0, S2, S02, Field, HarmonicTransform

# traits used for convert and promote
include("convert_promote.jl")
export harmonic_transform, HasQU, has_qu, IsMap, is_map, IsLenseBasis, is_lense_basis

# Harmonic transforms
include("harmonic_transforms/real_2d_flat_fourier.jl")
include("harmonic_transforms/real_2d_flat_unitary_fourier.jl")
include("harmonic_transforms/real_2d_flat_ordinary_fourier.jl")
include("harmonic_transforms/real_1d_flat_fourier.jl")
include("harmonic_transforms/real_1d_flat_unitary_fourier.jl")
include("harmonic_transforms/real_1d_flat_ordinary_fourier.jl")
include("harmonic_transforms/complex_2d_flat_fourier.jl")
# export r𝔽, 𝔽, r𝕌𝔽2, r𝕆𝔽2, r𝔽1, r𝕌𝔽1, r𝕆𝔽1

# field operations
include("field_ops.jl")

# linear ops
include("linear_ops.jl")
export LinOp, DiagOp # , 𝕃 # TODO get rid of 𝕃, in favor of DiagOp

# misc
include("util.jl")
export data, squash, white_noise


#TODO: dot for Pix <: Healpix in field_ops.jl
#TODO: get a spherical example up and running

end # end Module
