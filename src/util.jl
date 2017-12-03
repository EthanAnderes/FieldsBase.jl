
@generated data(x::Field) = :(tuple($((:(x.$f) for f=fieldnames(x))...)))

squash(x::T) where T = ifelse(isfinite(x), x, zero(T))

function white_noise(g::HarmonicTransform{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n,n) ./ sqrt(g.Ωx)
end
