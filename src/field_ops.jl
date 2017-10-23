
################################################
# basic operations on fields
################################################

##### fields with scalars

(+)(f::T, n::Number) where T<:Field = T((data(f) .+ n)...)
(+)(n::Number, f::T) where T<:Field = T((data(f) .+ n)...)
(-)(a::T)            where T<:Field = T((.- data(a))...)
(-)(f::T, n::Number) where T<:Field = T((data(f) .- n)...)
(-)(n::Number, f::T) where T<:Field = T((n .- data(f))...)
(*)(f::T, n::Number) where T<:Field = T((n .* data(f))...)
(*)(n::Number, f::T) where T<:Field = T((n .* data(f))...)

##### op(fields, fields)

# operators which broadcast to the underlying data
for op in (:+, :-, :*)
    @eval ($op)(a::T, b::T) where {T<:Field} = T(map((a,b)->broadcast($op,a,b),data(a),data(b))...)
end

# operators for which we do automatic promotion
for op in (:+, :-, :dot)
    @eval ($op)(a::Field, b::Field) = ($op)(promote(a,b)...)
end

# dot for Pix <: Flat
dot(a::X, b::X) where X<:Field{P} where P<:Flat = _dot(a, b, is_map(X))

# dot(map, map)
function _dot(a::X, b::X, ::Type{IsMap{true}}) where X<:Field{P} where P<:Flat
    FFT = harmonic_transform(X)
    return sum(map(vecdot, data(a),data(b))) * FFT.Ωpix
end

# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:Field{P} where P<:Flat
    FFT = harmonic_transform(X)
    sum(map(_rfftdot, data(a),data(b))) * FFT.Ωk
end

# Dots two matrices assumed to be the output of a real-FFT, s.t. if
# A/B::Matrix{<:Real}, `real(vecdot(fft(A),fft(B))) == _rfftdot(rfft(A),rfft(B))`
function _rfftdot(a,b)
    n,m = size(a)
    @assert size(a)==size(b) && n==m÷2+1
    real(vecdot(a[1,:],b[1,:]) + 2vecdot(a[2:end-1,:],b[2:end-1,:]) + (iseven(m) ? 1 : 2) * vecdot(a[end,:],b[end,:]))
end

# TODO dot for Pix <: Healpix
