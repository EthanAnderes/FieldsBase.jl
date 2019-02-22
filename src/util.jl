


@generated data(x::Field) = :(tuple($((:(x.$f) for f=fieldnames(x))...)))

@inline squash(x::T) where T = ifelse(isfinite(x), x, T(0))

function white_noise(g::HarmonicTransform{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n,n) ./ sqrt(g.Ωx)
end

function zero_map_dim2(::Type{P}, ::Type{T}) where {θ,n,P<:Flat{θ,n},T<:Real}
	fill(T(0), n, n)
end
function zero_map_dim1(::Type{P}, ::Type{T}) where {θ,n,P<:Flat{θ,n},T<:Real}
	fill(T(0), n)
end

function zero_fourier_dim2(::Type{P}, ::Type{T}) where {θ,n,P<:Flat{θ,n},T<:Real}
	fill(Complex{T}(0), n÷2+1, n)
end
function zero_fourier_dim1(::Type{P}, ::Type{T}) where {θ,n,P<:Flat{θ,n},T<:Real}
	fill(Complex{T}(0), n÷2+1)
end




