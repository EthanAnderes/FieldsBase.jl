
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
    return sum(map(_realdot, data(a),data(b))) * FFT.Ωx
end

# dot(fourier, fourier)
function _dot(a::X, b::X, ::Type{IsMap{false}}) where X<:Field{P} where P<:Flat
    FFT = harmonic_transform(X)
    sum(map(_complexdot, data(a),data(b))) * FFT.Ωk
end

# these work better for ArrayFire
_realdot(a,b) = sum(a.*b)

function _complexdot(a,b)
    n,m = size(a)
    @assert size(a)==size(b) && n==m÷2+1
    ra, ia = real(a), imag(a)
    rb, ib = real(b), imag(b)
    rtn  = 2*sum(ra.*rb)
    rtn += 2*sum(ia.*ib)
    rtn -= sum(ra[1,:].*rb[1,:])
    rtn -= sum(ia[1,:].*ib[1,:])
    if iseven(m)
           rtn -= sum(ra[end,:].*rb[end,:])
           rtn -= sum(ia[end,:].*ib[end,:])
    end
    return rtn
end
