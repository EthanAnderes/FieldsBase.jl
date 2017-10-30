################################################
# linear operators
#############################################

abstract type LinOp{P<:Pix, T<:Real, S<:Spin}  end

struct DiagOp{F<:Field,P,T,S} <: LinOp{P,T,S}
    f::F
    DiagOp(f::F) where {P,T,S,F<:Field{P,T,S}} = new{F,P,T,S}(f)
end

const 𝕃 = DiagOp

(*)(O::𝕃{F}, f::Field) where {F<:Field} = O.f * F(f)
(\)(O::𝕃{F}, f::Field) where {F<:Field} = inv(O) * f

# define 𝕃^a
(^)(op::𝕃{F}, a::Number)  where F<:Field = 𝕃(F((i.^a for i in data(op.f))...))
(^)(op::𝕃{F}, a::Integer) where F<:Field = 𝕃(F((i.^a for i in data(op.f))...))

# inv(𝕃), includes a pre-squash
inv(op::𝕃{F}) where {F<:Field} = 𝕃(F((squash.(1 ./ i) for i in data(op.f))...))
