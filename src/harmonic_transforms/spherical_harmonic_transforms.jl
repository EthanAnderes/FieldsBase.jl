

############################################
#  Healpix transform
#############################################

 struct ℍ{P<:Healpix,T<:Real}
    Ωpix::T
    lmax::T
    l::Matrix{T}
    m::Matrix{T}
    φ::Matrix{T}  # azmuth
    Θ::Matrix{T}  # polar
end

# function ebk_to_quk(ek, bk, ::Type{ℍ{P,T}}) where {P<:Pix, T<:Real}
#     g  = r𝔽(P,T)
#     qk = .- ek .* g.cos2ϕk .+ bk .* g.sin2ϕk
#     uk = .- ek .* g.sin2ϕk .- bk .* g.cos2ϕk
#     return qk, uk
# end
# function quk_to_ebk(qk, uk, ::Type{ℍ{P,T}}) where {P<:Pix, T<:Real}
#     g  = r𝔽(P,T)
#     ek = .- qk .* g.cos2ϕk .- uk .* g.sin2ϕk
#     bk =    qk .* g.sin2ϕk .- uk .* g.cos2ϕk
#     return ek, bk
# end
