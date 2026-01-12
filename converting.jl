#[ converting.jl ]#

# This file defines conversion functions between each type

#=
export index, x, y
=#

"""
###################################################################################################
                                        Operators
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################


# Operator()に通す -> 正規化される運命である ことの理解
# だから.opsを無理にいじることは推奨されていない(unsafe)

# 正規化される前を扱うコンストラクタ
function Operator()
    return Operator(Vector{AbstractOp}())
end
function Operator(op::AbstractOp)::Operator
    r = Operator()
    if op.cnt == 0
        return r
    end
    push!(r.ops,copy(op))   # copy内でOpDownはOpUpに直されている
    return r
end

# 最後にOperator()を通すことで正規化するコンストラクタ
function Operator(sym::AbstractString, n::Integer=1)::Operator
    if n == 0
        return Operator()
    end
    return Operator(str_to_op(sym,n))   # 正規化される
end
function Operator(v::Vector{<:Tuple{AbstractString, Integer}})::Operator
    regv = Vector{AbstractOp}(undef,lastindex(v))
    for vi in v
        regv[i] = str_to_op(vi) # UpDownの正規化がされていない
    end
    return Operator(regv)   # 正規化される
end
Operator(op::Operator)::Operator = op
Operator(op::Type{<:AbstractOp})::Operator = Operator((op)(1))
function Operator(op::Type{<:AbstractOp},cnt::Int,n::Int=1)::Operator
    if cnt == 0
        return Operator()
    end
    if op === OpDeriv
        return Operator(OpDeriv(cnt,n))
    else
        return Operator((op)(cnt))  # 正規化される
    end
end





"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################

# [============== about MonoIndex ==============]
IndexWord()::IndexWord = IndexWord(Tuple{}())
IndexWord(v::Vector)::IndexWord = IndexWord(Tuple{Vararg{Int}}(v))
function IndexWord(i::Index)::IndexWord
    if !is_monomial(i)
        throw(DomainError(i,"i must be monomial"))
    end
    return IndexWord(first(keys(i)))
end
function IndexWord(w::Hoffman)::IndexWord
    if !is_monomial(w)
        throw(DomainError(w,"w must be monomial"))
    end
    wo = first(keys(w))
    if get_index_orientation()
        if wo[1] != 2 # yに対応するだけ
            throw(DomainError(w,"w does not start with y"))
        end
        return IndexWord(idxprs(wo))
    else
        if wo[end] != 2
            throw(DomainError(w,"w does not end with y"))
        end
        return IndexWord(idxprs_r(wo))
    end
end
function IndexWord(w::HoffmanWord)::IndexWord
    if get_index_orientation()
        if w[1] != 2 # yに対応するだけ
            throw(DomainError(w,"w does not start with y"))
        end
        return IndexWord(idxprs(w))
    else
        if w[end] != 2
            throw(DomainError(w,"w does not end with y"))
        end
        return IndexWord(idxprs_r(w))
    end
end

IndexWord(a::Int...)::IndexWord = IndexWord(a)
IndexWord(m::IndexWord)::IndexWord = m

index(s...)::IndexWord = IndexWord(s...)    # = Index(s...)にするべき？
# TODO: MonoIndex for hrm shf mpl 


# [============== about Word ==============]
HoffmanWord()::HoffmanWord = HoffmanWord(Tuple{}())
HoffmanWord(s::Int...)::HoffmanWord = HoffmanWord(s)
const x = HoffmanWord(1)
const y = HoffmanWord(2)

HoffmanWord(v::Vector{Int})::HoffmanWord = HoffmanWord(Tuple{Vararg{Int}}(v))
function HoffmanWord(w::Hoffman)::HoffmanWord
    if !is_monomial(w)
        throw(DomainError(w,"w must be monomial"))
    end
    p = first(pairs(w))
    if !isone(p.second)
        throw(DomainError(w,"w's coefficient is not 1"))
    end
    return HoffmanWord(p.first)
end

function HoffmanWord(w::IndexWord)::HoffmanWord
    if get_index_orientation()
        return idxdprs(collect(w))
    else
        return idxdprs_r(collect(w))
    end
end
function HoffmanWord(i::Index)::HoffmanWord
    if !is_monomial(i)
        throw(DomainError(i,"i must be monomial"))
    end
    p = first(pairs(i))
    if !isone(p.second)
        throw(DomainError(i,"i's coefficient is not 1"))
    end
    return HoffmanWord(p.first)
end

HoffmanWord(w::HoffmanWord)::HoffmanWord = w
# TODO: Word for hrm shf mpl


# [============== about Index ==============]
Index()::Index = Index(Dict{IndexWord,Rational{BigInt}}())
function Index(idx::Int...)::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(idx),Rational(BigInt(1)) ) )
end
function Index(t::Tuple{Vararg{Int}})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(t),Rational(BigInt(1)) ) )
end
function Index(v::Vector{Int})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(v),Rational(BigInt(1)) ) )
end
function Index(c::NN, v::Vector{Int})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(v),Rational(BigInt(c)) ) )
end
function Index(v::Vector{Vector{Int}})::Index
    idx = Index()
    c = Rational(BigInt(1))
    for vi in v
        wvi = IndexWord(vi)
        idx[wvi] = getindex(idx,wvi) + c
    end
    return idx
end

# star-stuffle用 (Boolは+-がはいる)
function Index(v::Tuple{Vector{Vector{Int}}, Vector{Bool}})::Index
    idx = Index()
    c = Rational(BigInt(1))
    for i in 1:lastindex(v[1])
        wvi = IndexWord(v[1][i])
        if v[2][i]
            idx[wvi] = getindex(idx,wvi) - c
        else
            idx[wvi] = getindex(idx,wvi) + c
        end
    end
    return idx
end
function Index(m::IndexWord,coeff::NN = Rational(BigInt(1)))::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( m,Rational(BigInt(coeff)) ) )
end
function Index(w::HoffmanWord,coeff::NN = Rational(BigInt(1)))::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(w),Rational(BigInt(coeff)) ) )
end
function Index(w::Hoffman)::Index   # Hoffmanは正規化されている
    idx = Index()
    for (k,c) in w
        idx[IndexWord(k)] = c
    end
    return idx
end

Index(i::Index)::Index = i
# function Index(n::NN)::Index
#     idx = Index()
#     idx.terms[Word()] = n
#     return idx
# end
# TODO: Index for hrm shf mpl


# [============== about Hoffman ==============]
Hoffman()::Hoffman = Hoffman(Dict{HoffmanWord,Rational{BigInt}}())
function Hoffman(v::Vector{Vector{Int}})::Hoffman
    w = Hoffman()
    c = Rational(BigInt(1))
    for vi in v
        wvi = HoffmanWord(vi)
        w[wvi] = getindex(w,wvi) + c
    end
    return w
end
function Hoffman(m::IndexWord, coeff=Rational(BigInt(1)))::Hoffman
    return Hoffman( Dict{HoffmanWord,Rational{BigInt}}( HoffmanWord(m),Rational(BigInt(coeff)) ) )
end
function Hoffman(wm::HoffmanWord, coeff=Rational(BigInt(1))::NN)::Hoffman
    return Hoffman( Dict{HoffmanWord,Rational{BigInt}}( wm,Rational(BigInt(coeff)) ) )
end
function Hoffman(i::Index)::Hoffman # Indexは正規化されている
    w = Hoffman()
    for (k,c) in i
        w[HoffmanWord(k)] = c
    end
    return w
end

Hoffman(w::Hoffman)::Hoffman = w
# function Hoffman(n::NN)::Hoffman
#     w = Hoffman()
#     w.terms[Word()] = n
#     return w
# end
# TODO: Hoffman for hrm shf mpl


# [============== about Poly ==============]

poly_promote_rule(::Type{HoffmanWord}) = Hoffman
poly_promote_rule(::Type{Hoffman})     = Hoffman
poly_promote_rule(::Type{IndexWord})   = Index
poly_promote_rule(::Type{Index})       = Index
poly_promote_rule(::Type{NN})          = Rational{BigInt}
poly_promote_rule(::Type{<:Integer})   = Rational{BigInt}
poly_promote_rule(::Type{<:Rational})  = Rational{BigInt}
poly_promote_rule(::Type{T}) where T   = T

# 汎用
Poly{A}() where A = Poly{A}(Dict{Int,A}())
function Poly(w::A) where A
    pA = poly_promote_rule(A)
    r = Poly{pA}()
    r[0] = pA(w)
    return r
end

# 自己
function Poly(r::Poly{A})::Poly{A} where A
    return r
end

# 相互
function Hoffman(a::Poly{Index})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in a
        r[d] = Hoffman(c)
    end
    return r
end
Hoffman(a::Poly{Hoffman})::Poly{Hoffman} = a
function Index(a::Poly{Hoffman})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in a
        r[d] = Index(c)
    end
    return r
end
Index(a::Poly{Index})::Poly{Index} = a
function Hoffman(a::Poly{Rational{BigInt}})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in a
        r[d] = Hoffman(HoffmanWord(),c)
    end
    return r
end
function Index(a::Poly{Rational{BigInt}})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in a
        r[d] = Index(IndexWord(),c)
    end
    return r
end

# TODO: Poly for hrm shf mpl

