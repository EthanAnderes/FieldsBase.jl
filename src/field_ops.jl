
################################################
# basic operations on fields
################################################

# operators which broadcast to the underlying data
for op in (:+, :-, :*)
    @eval ($op)(a::T, b::T) where {T<:Field} = T(map((a,b)->broadcast($op,a,b),data(a),data(b))...)
end

# operators for which we do automatic promotion
for op in (:+, :-, :dot)
    @eval ($op)(a::Field, b::Field) = ($op)(promote(a,b)...)
end

# fields with scalars
(+)(f::T, n::Number) where T<:Field = T((data(f) .+ n)...)
(+)(n::Number, f::T) where T<:Field = T((data(f) .+ n)...)
(-)(a::T)            where T<:Field = T((.- data(a))...)
(-)(f::T, n::Number) where T<:Field = T((data(f) .- n)...)
(-)(n::Number, f::T) where T<:Field = T((n .- data(f))...)
(*)(f::T, n::Number) where T<:Field = T((n .* data(f))...)
(*)(n::Number, f::T) where T<:Field = T((n .* data(f))...)
