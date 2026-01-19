#[ arithmetic.jl ]#

# This file defines arithmetic functions

import Base: +, -, *, ^, //

#=
export shift_degree
=#

"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Arithmetic Functions ###############################################################

# isdefined table
#=

  Addition and Subtraction
--------------------------------------------------------------------------------------------------------------------------
|  Add/Sub       |  NN  |  HoffmanWord  |  Hoffman  |  IndexWord  |  Index  |  Poly NN   |  Poly Hoffman  |  Poly Index  |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  NN            |  NN  |  Hof          |  Hof      |  Idx        |  Idx    |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  HofW          |      |  Hof          |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  HoffmanWord   |      |               |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  IndexWord     |      |               |           |  Idx        |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Index         |      |               |           |             |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly NN       |      |               |           |             |         |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Hoffman  |      |               |           |             |         |            |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Index    |      |               |           |             |         |            |                |  Poly Idx    |
--------------------------------------------------------------------------------------------------------------------------

=#



############################## ADDITIVE INVERSE ##############################

# HoffmanWord
(-)(a::HoffmanWord)::Hoffman = Hoffman(a, -1)
(+)(a::HoffmanWord)::HoffmanWord = a

# Hoffman
function (-)(a::Hoffman)::Hoffman
    result = Hoffman()
    for (w, c) in a
        result.terms[w] = -c
    end
    return result
end
(+)(a::Hoffman)::Hoffman = a

# IndexWord
(-)(a::IndexWord)::IndexWord = Index(a, -1)
(+)(a::IndexWord)::IndexWord = a

# Index
function (-)(a::Index)::Index
    result = Index()
    for (w, c) in a
        result.terms[w] = -c
    end
    return result
end
(+)(a::Index)::Index = a

# Poly
function (-)(a::Poly{A})::Poly{A} where A
    r = copy(a)
    for (d, h) in a
        r.terms[d] = -h
    end
    return r
end
function (+)(a::Poly{A})::Poly{A} where A
    return a
end


############################## typelift ##############################

typelift(::Type{A}, a::A) where {A} = a

typelift(::Type{Rational{BigInt}},       c::NN)::NN                     = Rational(BigInt(c))
typelift(::Type{Hoffman},                c::NN)::Hoffman                = Hoffman( Dict{HoffmanWord,Rational{BigInt}}(HoffmanWord() => c) )
typelift(::Type{Index},                  c::NN)::Index                  = Index( Dict{IndexWord,Rational{BigInt}}(IndexWord() => c) )
typelift(::Type{Poly{Rational{BigInt}}}, c::NN)::Poly{Rational{BigInt}} = Poly( Dict{Int,Rational{BigInt}}(0 => c) )
typelift(::Type{Poly{Hoffman}},          c::NN)::Poly{Hoffman}          = Poly( Dict{Int,Hoffman}(0 => typelift(Hoffman,c)) )
typelift(::Type{Poly{Index}},            c::NN)::Poly{Index}            = Poly( Dict{Int,Index}(0 => typelift(Index,c)) )

typelift(::Type{Hoffman},      w::HoffmanWord)::Hoffman       = Hoffman( Dict{HoffmanWord,Rational{BigInt}}(w => Rational(BigInt(1))) )
typelift(::Type{Poly{Hoffman}},w::HoffmanWord)::Poly{Hoffman} = Poly( Dict{Int,Hoffman}(0 => typelift(Hoffman,w)) )

typelift(::Type{Poly{Hoffman}},w::Hoffman)::Poly{Hoffman} = Poly( Dict{Int,Hoffman}(0 => copy(w)) )

typelift(::Type{Index},        w::IndexWord)::Index       = Index( Dict{IndexWord,Rational{BigInt}}(w => Rational(BigInt(1))) )
typelift(::Type{Poly{Index}},  w::IndexWord)::Poly{Index} = Poly( Dict{Int,Index}(0 => typelift(Index,w)) )

typelift(::Type{Poly{Index}},  w::Index)::Poly{Index} = Poly( Dict{Int,Index}(0 => copy(w)) )

function typelift(::Type{Poly{Hoffman}},w::Poly{Rational{BigInt}})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in w
        r.terms[d] = typelift(Hoffman,c)
    end
    return r
end
function typelift(::Type{Poly{Index}},w::Poly{Rational{BigInt}})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in w
        r.terms[d] = typelift(Index,c)
    end
    return r
end


############################## ADD and SUBTRACT ##############################
const HNN = Union{HoffmanWord, NN}


########## HoffmanWord ##########

# HoffmanWord NN
+(a::HNN, b::HNN)::Hoffman = +(typelift(Hoffman,a),typelift(Hoffman,b))
-(a::HNN, b::HNN)::Hoffman = -(typelift(Hoffman,a),typelift(Hoffman,b))

########## Hoffman ##########

# Hoffman  NN,HoffmanWord
+(a::Hoffman, b::HNN)::Hoffman = +(a,typelift(Hoffman,b))
+(a::HNN, b::Hoffman)::Hoffman = +(typelift(Hoffman,a),b)
-(a::Hoffman, b::HNN)::Hoffman = -(a,typelift(Hoffman,b))
-(a::HNN, b::Hoffman)::Hoffman = -(typelift(Hoffman,a),b)

# Hoffman Hoffman
function +(a::Hoffman, b::Hoffman)::Hoffman
    result = copy(a)
    for (w, c) in b
        new = getindex(result,w) + c
        if iszero(new)
            delete!(result.terms, w)
        else
            result.terms[w] = new
        end
    end
    return result
end
function -(a::Hoffman, b::Hoffman)::Hoffman
    result = copy(a)
    for (w, c) in b
        new = getindex(result,w) - c
        if iszero(new)
            delete!(result.terms, w)
        else
            result.terms[w] = new
        end
    end
    return result
end

########## IndexWord ##########
const INN = Union{IndexWord, NN}

# IndexWord NN
+(a::INN, b::INN)::Index = +(typelift(Index,a),typelift(Index,b))
-(a::INN, b::INN)::Index = -(typelift(Index,a),typelift(Index,b))

########## Index ##########

# Index  NN,IndexWord
+(a::Index, b::INN)::Index = +(a,typelift(Index,b))
+(a::INN, b::Index)::Index = +(typelift(Index,a),b)
-(a::Index, b::INN)::Index = -(a,typelift(Index,b))
-(a::INN, b::Index)::Index = -(typelift(Index,a),b)

# Index Index
function +(a::Index, b::Index)::Index
    result = copy(a)
    for (w, c) in b
        new = getindex(result,w) + c
        if iszero(new)
            delete!(result.terms, w)
        else
            result.terms[w] = new
        end
    end
    return result
end
function -(a::Index, b::Index)::Index
    result = copy(a)
    for (w, c) in b
        new = getindex(result,w) - c
        if iszero(new)
            delete!(result.terms, w)
        else
            result.terms[w] = new
        end
    end
    return result
end

########## Poly ##########
const HHNN = Union{Hoffman, HoffmanWord, NN}
const IINN = Union{Index,   IndexWord,   NN}
const RNN  = Union{Rational{BigInt},     NN}

poly_typelift(::Type{<:RNN}, ::Type{<:RNN} ) = Rational{BigInt}
poly_typelift(::Type{<:HHNN},::Type{<:HHNN}) = Hoffman
poly_typelift(::Type{<:IINN},::Type{<:IINN}) = Index

function +(a::Poly{A},b::Poly{B}) where {A,B}
    C = poly_typelift(A,B)
    result = Poly{C}()
    for (deg,h) in a
        result.terms[deg] = typelift(C,h)
    end
    for (d, c) in b
        new = getindex(result,d) + c
        if iszero(new)
            delete!(result.terms, d)
        else
            result.terms[d] = new
        end
    end
    return result
end
function -(a::Poly{A},b::Poly{B}) where {A,B}
    C = poly_typelift(A,B)
    result = Poly{C}()
    for (deg,h) in a
        result.terms[deg] = typelift(C,h)
    end
    for (d, c) in b
        new = getindex(result,d) - c
        if iszero(new)
            delete!(result.terms, d)
        else
            result.terms[d] = new
        end
    end
    return result
end
function +(a::Poly{A},b::B) where {A,B}
    C = poly_typelift(A,B)
    result = Poly{C}()
    for (deg,h) in a
        result.terms[deg] = typelift(C,h)
    end
    new = getindex(result,0) + typelift(C,b)
    if iszero(new)
        delete!(result.terms,0)
    else
        result.terms[0] = new
    end
    return result
end
function -(a::Poly{A},b::B) where {A,B}
    C = poly_typelift(A,B)
    result = Poly{C}()
    for (deg,h) in a
        result.terms[deg] = typelift(C,h)
    end
    new = getindex(result,0) - typelift(C,b)
    if iszero(new)
        delete!(result.terms,0)
    else
        result.terms[0] = new
    end
    return result
end
+(a::A,b::Poly{B}) where {A,B} = +(b,a)
-(a::A,b::Poly{B}) where {A,B} = +(a,-b)


############################## MULTIPLICATION ##############################

# general (今後へ残しておく)
*(a::HHNN,b::HHNN) = *(typelift(Hoffman,a),typelift(Hoffman,b))
*(a::IINN,b::IINN) = *(typelift(Index  ,a),typelift(Index  ,b))

# specific
@inline function *(a::HoffmanWord, b::HoffmanWord)::HoffmanWord
    return HoffmanWord(a... , b...)
end
@inline function *(a::IndexWord, b::IndexWord)::IndexWord
    return IndexWord(a... , b...)
end
function *(a::NN, b::HoffmanWord)::Hoffman
    r = Hoffman()
    if a != 0
        r.terms[b] = a
    end
    return r
end
*(a::HoffmanWord,b::NN)::Hoffman = *(b,a)
function *(a::NN, b::IndexWord)::Index
    r = Index()
    if a != 0
        r.terms[b] = a
    end
    return r
end
*(a::IndexWord,b::NN)::Index = *(b,a)
function *(a::NN, b::Hoffman)::Hoffman
    r = Hoffman()
    for (w,c) in b
        r.terms[w] = c*a
    end
    return r
end
*(a::Hoffman, b::NN)::Hoffman = *(b,a)
function *(a::NN, b::Index)::Index
    r = Index()
    for (w,c) in b
        r.terms[w] = c*a
    end
    return r
end
*(a::Index, b::NN)::Index = *(b,a)
function *(a::HoffmanWord, b::Hoffman)::Hoffman
    r = Hoffman()
    for (w,c) in b
        r.terms[a*w] = c
    end
    return r
end
function *(a::Hoffman, b::HoffmanWord)::Hoffman
    r = Hoffman()
    for (w,c) in a
        r.terms[w*b] = c
    end
    return r
end
function *(a::IndexWord, b::Index)::Index
    r = Index()
    for (w,c) in b
        r.terms[a*w] = c
    end
    return r
end
function *(a::Index, b::IndexWord)::Index
    r = Index()
    for (w,c) in a
        r.terms[w*b] = c
    end
    return r
end

# base
function *(a::Hoffman, b::Hoffman)::Hoffman
    r = Hoffman()
    for (wa, ca) in a
        for (wb, cb) in b
            w = wa*wb
            r.terms[w] = getindex(r,w) + ca*cb
        end
    end
    filter!(p->!iszero(p.second),r.terms)
    return r
end
function *(a::Index, b::Index)::Index
    r = Index()
    for (wa, ca) in a
        for (wb, cb) in b
            w = wa*wb
            r.terms[w] = getindex(r,w) + ca*cb
        end
    end
    filter!(p->!iszero(p.second),r.terms)
    return r
end



########## Poly ##########


function *(a::Poly{A}, b::Poly{B}) where {A,B}
    C = poly_typelift(A,B)  # ここで積のmethodが無ければエラー

    # 単項最適化
    if is_monomial(a)
        (d,c) = first(a.terms)
        if isone(c)
            return shift_degree(typelift(Poly{C}, b), d)
        end
    end
    if is_monomial(b)
        (d,c) = first(b.terms)
        if isone(c)
            return shift_degree(typelift(Poly{C}, a), d)
        end
    end

    # 通常の積
    r = Poly{C}()
    for (da, ha) in a
        for (db, hb) in b
            d = da + db
            r.terms[d] = getindex(r,d) + ha*hb
        end
    end
    filter!(p->!iszero(p.second), r.terms)
    return r
end
function *(a::Poly{A}, b::B) where {A,B}
    C = poly_typelift(A,B)  # ここで積のmethodが無ければエラー
    r = Poly{C}()
    for (d,h) in a
        r.terms[d] = h*b
    end
    return r
end
*(a::A, b::Poly{B}) where {A,B} = *(b,a)


############################## POWER ##############################
# HoffmanWord
function ^(t::HoffmanWord, n::Integer)::HoffmanWord
    n <= 0 && return HoffmanWord()
    len = length(t)
    total = len * n
    v = Vector{eltype(t)}(undef, total)
    for i in 0:n-1
        copyto!(v, i*len + 1, t, 1, len)
    end
    return HoffmanWord(v)
end
# IndexWord
function ^(t::IndexWord, n::Integer)::IndexWord
    n <= 0 && return IndexWord()
    len = length(t)
    total = len * n
    v = Vector{eltype(t)}(undef, total)
    for i in 0:n-1
        copyto!(v, i*len + 1, t, 1, len)
    end
    return IndexWord(v)
end

# Hoffman
function ^(a::Hoffman, n::Integer)::Hoffman
    if n < 0
        throw(ArgumentError("Hoffman-type powers for negative exponents are not defined"))
    elseif n == 0
        return one(Hoffman)
    elseif n == 1
        return copy(a)
    elseif isone(a)
        return one(Hoffman)
    end

    result = one(Hoffman)
    base = copy(a)
    nn = n

    while nn > 0
        if (nn & 1) == 1
            result = result * base
        end
        nn >>= 1
        if nn > 0
            base = base * base
        end
    end

    return result
end

# Index
function ^(a::Index, n::Integer)::Index
    if n < 0
        throw(ArgumentError("Index-type powers for negative exponents are not defined"))
    elseif n == 0
        return one(Index)
    elseif n == 1
        return copy(a)
    elseif isone(a)
        return one(Index)
    end

    result = one(Index)
    base = copy(a)
    nn = n

    while nn > 0
        if (nn & 1) == 1
            result = result * base
        end
        nn >>= 1
        if nn > 0
            base = base * base
        end
    end

    return result
end

# Poly
function ^(a::Poly{A}, n::Integer)::Poly{A} where A
    if n < 0
        throw(DomainError(n, "negative power is not supported for Poly"))
    elseif n == 0
        return one(Poly{A})
    elseif n == 1
        return copy(a)
    end

    # 最適化: a が単項 T^k (係数が exactly one(Hoffman)) の場合
    if is_monomial(a)
        (deg, coeff) = first(a.terms)
        r = Poly{A}()
        r.terms[deg * Int(n)] = coeff ^ n
        return r
    end

    # 二乗累乗（binary exponentiation）
    base = copy(a)
    result = one(Poly{A})
    nn = Int(n)  # n は非負なので安全に Int に落とす

    while nn > 0
        if (nn & 1) == 1
            result = result * base   # 既に定義済みの *(Poly, Poly) を使用
        end
        nn >>= 1
        if nn > 0
            base = base * base
        end
    end

    return result
end

# degree shift
""" r , n -> r T^n """
function shift_degree(r::Poly{A},n::Int)::Poly{A} where A
    out = Poly{A}()
    for (d, h) in r.terms
        out.terms[d+n] = (h isa Hoffman || h isa Index ? copy(h) : h)
    end
    return out
end


############################# DIVISION for NN ##############################
# NN Word
function //(b::HoffmanWord, a::NN)::Hoffman
    r = Hoffman()
    r.terms[b] = 1//a
    return r
end

# NN Hoffman
function //(b::Hoffman, a::NN)::Hoffman
    r = Hoffman()
    for (w, c) in b.terms
        r.terms[w] = c//a
    end
    return r
end

# NN IndexWord
function //(b::IndexWord, a::NN)::Index
    r = Index()
    r.terms[b.word] = b.coeff//a
    return r
end

# NN Index
function //(b::Index, a::NN)::Index
    r = Index()
    for (w, c) in b.terms
        r.terms[w] = c//a
    end
    return r
end

# NN Poly
function //(b::Poly{A}, a::NN)::Poly{A} where A
    r = Poly{A}()
    for (d, h) in b.terms
        r.terms[d] = h//a
    end 
    return r
end