################################################
# linear operators
#############################################

abstract type LinOp{P<:Pix, S<:Spin}  end

struct DiagOp{F<:Field,P,S} <: LinOp{P,S}
    f::F
    DiagOp(f::F) where {P,S,F<:Field{P,S}} = new{F,P,S}(f)
end

const 𝕃 = DiagOp

*(O::ℒ{F}, f::Field) where {F} = O.f * F(f)


# # define ℒ^a
# (^)(op::ℒ{F}, a::Number)  where F<:Field = ℒ(F((i.^a for i in data(op.f))...))
# (^)(op::ℒ{F}, a::Integer) where F<:Field = ℒ(F((i.^a for i in data(op.f))...))
#
# # inv(ℒ), includes a pre-squash
# inv(op::ℒ{F}) where {F<:Field} = ℒ(F((squash.(1 ./ i) for i in data(op.f))...))
