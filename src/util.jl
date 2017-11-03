

@generated data(x::Field) = :(tuple($((:(x.$f) for f=fieldnames(x))...)))

squash(x::T) where {T<:Number}  = isnan(x) ? zero(T) : isfinite(x) ? x : zero(T)

function white_noise(g::HarmonicTransform{P,T}) where {T,P<:Flat{θ,n}} where {θ,n}
    randn(T,n,n) ./ sqrt(g.Ωx)
end
# function white_noise(::Type{P}, ::Type{T}) where {Θ, n, P<:Flat{Θ,n}, T<:Real}
#     randn(T,n,n) ./ r𝔽(P,T).Δx #<-- e.g. to be used as input for Imap{Px,Tx}(...)
# end
