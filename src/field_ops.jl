
################################################
# basic operations on fields
################################################


##### fields with scalars

(+)(f::F, n::Number) where F<:Field{P,T} where {P<:Pix,T} = F((data(f) .+ T(n))...)
(+)(n::Number, f::F) where F<:Field{P,T} where {P<:Pix,T} = F((data(f) .+ T(n))...)
(-)(a::F)            where F<:Field{P,T} where {P<:Pix,T} = F((.- data(a))...)
(-)(f::F, n::Number) where F<:Field{P,T} where {P<:Pix,T} = F((data(f) .- T(n))...)
(-)(n::Number, f::F) where F<:Field{P,T} where {P<:Pix,T} = F((T(n) .- data(f))...)
(*)(f::F, n::Number) where F<:Field{P,T} where {P<:Pix,T} = F((T(n) .* data(f))...)
(*)(n::Number, f::F) where F<:Field{P,T} where {P<:Pix,T} = F((T(n) .* data(f))...)

##### fields with UniformScaling

(*)(f::F, J::UniformScaling) where F<:Field{P,T} where {P<:Pix,T} = J.λ * f
(*)(J::UniformScaling, f::F) where F<:Field{P,T} where {P<:Pix,T} = J.λ * f

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
    # ----- this kills the repeated frequencies
    multip_ri = fill(true, n, m)   # both real and imaginary
    multip_ri[1,(n+1):m] .= false  # repeats
    if iseven(m)
        multip_ri[end,(n+1):m] .= false  # repeats
    end
    # ----- Now to a direct (x 2) dot product of the real and imag parts
    ra, ia = real(a), imag(a)
    rb, ib = real(b), imag(b)
    rtn  = 2*dot(ra, multip_ri .* rb) + 2*dot(ia, multip_ri .* ib)
    # ----- but we don't mult by 2 for the real terms... so we take away 1
    rtn -= ra[1,1]*rb[1,1]  
    rtn -= ra[1,n]*rb[1,n] 
    if iseven(m) 
        rtn -= ra[n,1]*rb[n,1] 
        rtn -= ra[n,n]*rb[n,n] 
    end 
    return rtn
end
