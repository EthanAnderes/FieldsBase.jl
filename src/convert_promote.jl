
harmonic_transform(::Type{Any}) = error("define harmonic_transform")
# this needs to be set by the user....

#############################################
##  THTT (Tim Holy Trait Trick)
#############################################
# NOTE: any user defined field needs to specify `has_qu(⋅), is_map(⋅), is_lense_basis(⋅)`

struct HasQU{Bool} end
# fallback
has_qu(::Type{X}) where X = error("no definition of has_qu") # --> HasQU{true/false}

struct IsMap{Bool} end
# fallback
is_map(::Type{X}) where X = error("no definition of is_map") # --> IsMap{true/false}

struct IsLenseBasis{Bool} end
# fallback
is_lense_basis(::Type{X}) where X = error("no definition of is_lense_basis") # --> IsLenseBasis{true/false}




#############################################
## Convert and Promote using Traits
#############################################

convert(::Type{X}, f::X) where X<:Field = X((d for d in data(f))...)
convert(::Type{X}, f::Y) where {X<:Field, Y<:Field} = _convert(X, has_qu(X), is_map(X), f, has_qu(Y), is_map(Y))

# convert(Xfourier,  f::Xmap)  for X ∈ {T, QU, EB, TQU, TEB}
function _convert(::Type{X}, ::Type{T}, ::Type{IsMap{false}}, f::Y, ::Type{T}, ::Type{IsMap{true}}) where {T<:HasQU, X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    X((harmonic_transform(X) * d for d in data(f))...)
end

# convert(Xmap,  f::Xfourier)  for X ∈ {T, QU, EB, TQU, TEB}
function _convert(::Type{X}, ::Type{T}, ::Type{IsMap{true}},  f::Y, ::Type{T}, ::Type{IsMap{false}}) where {T<:HasQU, X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    X((harmonic_transform(X) \ d for d in data(f))...)
end


############### convert(QUx, f::EBy) for x,y ∈ {fourier, map}

#convert(QUfourier,  f::EBfourier)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{false}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    ek, bk = data(f)
    X(ebk_to_quk(ek, bk, FT)...)
end

# convert(QUfourier,  f::EBmap)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{true}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    ex, bx = data(f)
    X(ebk_to_quk(FT * ex, FT * bx, FT)...)
end
#convert(QUmap,  f::EBfourier)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{false}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    ek, bk = data(f)
    qk, uk = ebk_to_quk(ek, bk, FT)
    X(FT \ qk, FT \ uk)
end
# convert(QUmap,  f::EBmap)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{true}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    ex, bx = data(f)
    qk, uk = ebk_to_quk(FT * ex, FT * bx, FT)
    X(FT \ qk, FT \ uk)
end

################# convert(EBx, f::QUy) for x,y ∈ {fourier, map}

# convert(EBfourier, f::QUfourier)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{false}},  f::Y, ::Type{HasQU{true}}, ::Type{IsMap{false}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    qk, uk = data(f)
    X(ebk_to_quk(qk, uk, FT)...)
end
# convert(EBfourier, f::QUmap)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{true}}, ::Type{IsMap{true}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    qx, ux = data(f)
    X(ebk_to_quk(FT * qx, FT * ux, FT)...)
end
#convert(EBmap, f::QUfourier)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{true}},  f::Y, ::Type{HasQU{true}}, ::Type{IsMap{false}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    qk, uk = data(f)
    ek, bk = ebk_to_quk(qk, uk, FT)
    X(FT \ ek, FT \ bk)
end
# convert(EBmap, f::QUmap)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{true}}, ::Type{IsMap{true}}) where {X<:Field{Px,S2},Y<:Field{Px,S2}} where Px<:Pix
    FT = harmonic_transform(X)
    qx, ux = data(f)
    ek, bk = ebk_to_quk(FT * qx, FT * ux, FT)
    X(FT \ ek, FT \ bk)
end

################# convert(IQUx, f::IEBy) for x,y ∈ {fourier, map}

# convert(IQUfourier,  f::IEBfourier)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{false}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tk, ek, bk = data(f)
    X(tk, ebk_to_quk(ek, bk, FT)...)
end
# convert(IQUfourier,  f::IEBmap)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{true}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tx, ex, bx = data(f)
    X(FT * tx, ebk_to_quk(FT * ex, FT * bx, FT)...)
end
#convert(IQUmap,  f::IEBfourier)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{false}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tk, ek, bk = data(f)
    qk, uk = ebk_to_quk(ek, bk, FT)
    X(FT \ tk, FT \ qk, FT \ uk)
end
# convert(IQUmap,  f::IEBmap)
function _convert(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{false}}, ::Type{IsMap{true}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tx, ex, bx = data(f)
    qk, uk = ebk_to_quk(FT \ ex, FT \ bx, FT)
    X(tx, FT \ qk, FT \ uk)
end

################# convert(IEBx, f::IQUy) for x,y ∈ {fourier, map}

# convert(IEBfourier, f::IQUfourier)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{false}},  f::Y, ::Type{HasQU{true}}, ::Type{IsMap{false}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tk, qk, uk = data(f)
    X(tk, ebk_to_quk(qk, uk, FT)...)
end
# convert(IEBfourier, f::IQUmap)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{false}}, f::Y, ::Type{HasQU{true}}, ::Type{IsMap{true}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tx, qx, ux = data(f)
    X(FT * tx, ebk_to_quk(FT * qx, FT * ux, FT)...)
end
#convert(IEBmap, f::IQUfourier)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{true}},  f::Y, ::Type{HasQU{true}}, ::Type{IsMap{false}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tk, qk, uk = data(f)
    ek, bk = ebk_to_quk(qk, uk, FT)
    X(FT \ tk, FT \ ek, FT \ bk)
end
# convert(IEBmap, f::IQUmap)
function _convert(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{true}}, f::Y, ::Type{HasQU{true}}, ::Type{IsMap{true}}) where {X<:Field{Px,S02},Y<:Field{Px,S02}} where Px<:Pix
    FT = harmonic_transform(X)
    tx, qx, ux = data(f)
    ek, bk = ebk_to_quk(FT * qx, FT * ux, FT)
    X(tx, FT \ ek, FT \ bk)
end



###################################################
##  Promotion with traits
###################################################

promote_rule(::Type{X}, ::Type{Y}) where {X<:Field, Y<:Field} = _promote_rule(X, has_qu(X), is_map(X), Y, has_qu(Y), is_map(Y))

#### Xfourier wins over Ymap for X,Y ∈ {T, QU, EB, TQU, TEB}
function _promote_rule(::Type{X}, ::Type{T1}, ::Type{IsMap{false}}, ::Type{Y}, ::Type{T2}, ::Type{IsMap{true}}) where {T1, T2, X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return X
end
function _promote_rule(::Type{X}, ::Type{T1}, ::Type{IsMap{true}}, ::Type{Y}, ::Type{T2}, ::Type{IsMap{false}}) where {T1, T2, X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return Y
end

#### QUmap wins over EBmap
function _promote_rule(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{true}}, ::Type{Y}, ::Type{HasQU{false}}, ::Type{IsMap{true}}) where {X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return X
end
function _promote_rule(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{true}}, ::Type{Y}, ::Type{HasQU{true}}, ::Type{IsMap{true}}) where {X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return Y
end

#### EBfourier wins over QUfourier
function _promote_rule(::Type{X}, ::Type{HasQU{true}}, ::Type{IsMap{false}}, ::Type{Y}, ::Type{HasQU{false}}, ::Type{IsMap{false}}) where {X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return Y
end
function _promote_rule(::Type{X}, ::Type{HasQU{false}}, ::Type{IsMap{false}}, ::Type{Y}, ::Type{HasQU{true}}, ::Type{IsMap{false}}) where {X<:Field{Px},Y<:Field{Px}} where Px<:Pix
    return X
end
