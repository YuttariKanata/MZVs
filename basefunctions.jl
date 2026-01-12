#[ basefunctions.jl ]#

# This file defines the basic functions
import Base: iszero, isone, zero, one, copy, ==, lastindex
import Base: getindex, length, iterate, firstindex, eltype, isempty, copy, collect, Tuple, vcat, hcat, reverse, ==, hash, isless, isequal, convert, in, iterate, haskey, keys, values, pairs

#=
export is_monomial, is_hoffman, is_index, is_monoindex,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr,is_zetaexpr
       isadmissible,
=#

"""
###################################################################################################
                                        Arithmatics
###################################################################################################
"""
#=
A
C(n+1,k+1) = (n+1)/(k+1)*C(n,k)
B
C(n+1,k) = (n+1)/(n-k+1)*C(n,k)
=#
#= (遅い)
function multinomial(v::Vector{Int})
    sorted_v = sort(v)
    binomm = Rational(BigInt(1))
    n = sorted_v[1]
    k = sorted_v[1]
    
    multinomm = Rational(BigInt(1))

    for i in 1:lastindex(v)-1
        for _ in 1:(sorted_v[i+1]-sorted_v[i])
            n += 1
            k += 1
            binomm *= n//k
        end
        for _ in 1:sorted_v[i]
            n += 1
            binomm *= n//(n-k)
        end
        multinomm *= binomm
    end

    return multinomm
end
=#
@inline function multinomial(v::Vector{Int})::BigInt
    num = factorial(BigInt(sum(v)))
    den = BigInt(1)
    for k in v
        den *= factorial(BigInt(k))
    end
    return div(num,den)
end
"""
v,n -> n/(v[1]!v[2]! ...)
"""
@inline function multinomial(v::Vector{Int},n::BigInt)::BigInt
    den = BigInt(1)
    for k in v
        den *= factorial(BigInt(k))
    end
    return div(n,den)
end


"""
###################################################################################################
                                        Operators
###################################################################################################
"""

###################################################################################################
############## Property Functions #################################################################


isone(op::Operator)::Bool = isempty(op.ops)
isone(op::Ts where Ts<:AbstractOp)::Bool  = (op.cnt == 0)

one(::Type{Operator})::Operator = Operator()
one(op::Type{<:AbstractOp})::AbstractOp = (op)(1)

==(a::Operator, b::Operator)::Bool = a.ops == b.ops
function ==(a::AbstractOp, b::AbstractOp)::Bool
    ta = typeof(a)
    tb = typeof(b)
    if ta === tb
        if ta === OpDeriv
            if a.cnt == b.cnt && a.n == b.n
                return true
            end
        else
            if a.cnt == b.cnt
                return true
            end
        end
        return false
    elseif (ta === OpUp && tb === OpDown) || (ta === OpDown && tb == OpUp)
        if a.cnt == -b.cnt
            return true
        end
    end
    return false
end
==(a::Operator, b::AbstractOp)::Bool = lastindex(a.ops) == 1 && a.ops[1] == b

lastindex(op::Operator)::Int = lastindex(op.ops)


###################################################################################################
############## Base Functions #####################################################################

# 文字->Operatorに変換するテーブル
const _String_to_Operator_Table = Dict{String, DataType}(
    "↑" => OpUp,
    "↓" => OpDown,
    "←" => OpLeft,
    "→" => OpRight,
    "-" => OpMinus,
    "τ" => OpTau,
    "†" => OpTau, # constで使おうとすると unknown unicode character '†'
    "𖼷" => OpTau, # ミャオ文字なので使える(泣)
    "η" => OpEta,
    "⋁" => OpEta,
    "φ" => OpPhi,
    "ϕ" => OpPhi,
    "∂" => OpDeriv,
)
# Operator->文字に変換するテーブル
const _Operator_to_String_Table = Dict{DataType, String}(
    OpUp => "↑",
    OpDown => "↓",
    OpLeft => "←",
    OpRight => "→",
    OpMinus => "-",
    OpTau => "τ",
    OpEta => "η",
    OpPhi => "φ",
    OpDeriv => "∂",
)

@inline function str_to_op(sym::AbstractString,n::Int = 1)::AbstractOp
    m = match(r"^∂(\d+)$", sym)
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(i,n)
    end
    return _String_to_Operator_Table[sym](n)
end
@inline function str_to_op(t::Tuple{AbstractString,Integer})::AbstractOp
    m = match(r"^∂(\d+)$", t[1])
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(i,t[2])
    end
    return _String_to_Operator_Table[t[1]](t[2])
end

# 単なる正規化
@inline function clean(vn::Vector{<:AbstractOp})::Vector{AbstractOp}
    
    vnl = lastindex(vn)
    if vnl == 0
        return Vector{AbstractOp}()
    end
    v = Vector{AbstractOp}(undef,vnl)
    for i in 1:vnl
        v[i] = copy(vn[i])  # copyでv内のOpDownはOpUpに直されている
    end

    r = Vector{AbstractOp}()
    sizehint!(r,vnl)

    # 常に op.ops内のUp,DownはすべてUpに直す
    for vi in v
        # vi:右側のAbstractOp (vn1,vn2,vn3,...)

        # bef:左側の最後    r1,r2,...,bef   <-append->   vn1,vn2,vn3,...
        bef = (isempty(r) ? nothing : r[end])

        tvi = typeof(vi)
        top = typeof(bef)

        if tvi !== top
            # 種類が違う場合は結合できない
            push!(r, vi)
        elseif tvi === OpTau
            # どちらもOpTauの場合 -> mod2 = xor(a,b) & 1
            r[end].cnt = xor(vi.cnt, bef.cnt) & 1
        elseif tvi === OpDeriv
            # どちらもOpDerivの場合 -> nが同じなら結合 そうでないなら追加
            if vi.n == bef.n
                r[end].cnt += vi.cnt
            else
                push!(r, vi)
            end
        else
            # それ以外の結合 -> cnt追加
            r[end].cnt += vi.cnt
        end

        # もしr[end]のcntが上の操作で0になったなら削除
        if !isempty(r) && r[end].cnt == 0
            pop!(r)
        end

    end

    return r
end
@inline function append_clean!(op::Operator, vn::Vector{<:AbstractOp})::Operator

    vnl = lastindex(vn)
    v = Vector{AbstractOp}(undef,vnl)
    for i in 1:vnl
        v[i] = copy(vn[i])  # copyでv内のOpDownはOpUpに直されている
    end

    # 常に op.ops内のUp,DownはすべてUpに直す
    for vi in v
        # vi:右側のAbstractOp (vn1,vn2,vn3,...)

        # bef:左側の最後    op1,op2,...,bef   <-append->   vn1,vn2,vn3,...
        bef = (isempty(op.ops) ? nothing : op.ops[end])

        tvi = typeof(vi)
        top = typeof(bef)

        if tvi !== top
            # 種類が違う場合は結合できない
            push!(op.ops, vi)
        elseif tvi === OpTau
            # どちらもOpTauの場合 -> mod2 = xor(a,b) & 1
            op.ops[end].cnt = xor(vi.cnt, bef.cnt) & 1
        elseif tvi === OpDeriv
            # どちらもOpDerivの場合 -> nが同じなら結合 そうでないなら追加
            if vi.n == bef.n
                op.ops[end].cnt += vi.cnt
            else
                push!(op.ops, vi)
            end
        else
            # それ以外の結合 -> cnt追加
            op.ops[end].cnt += vi.cnt
        end

        # もしop.ops[end]のcntが上の操作で0になったなら削除
        if !isempty(op.ops) && op.ops[end].cnt == 0
            pop!(op.ops)
        end

    end

    return op
end
@inline function append_clean!(a::Operator, bn::Operator)::Operator
    # 中央付近の結合だけ見ればいいという思想がある

    aftend = lastindex(bn.ops)
    b = copy(bn.ops)
    if isempty(b)
        return a
    end
    if isempty(a.ops)
        a.ops = b
        return a
    end

    # 中央境界
    #      1:a[1], 2:a[2], 3:a[3], ... , befidx:a[end]
    # aftidx:b[1], 2:b[2], 3:b[3], ...
    befidx = lastindex(a.ops)
    aftidx = 1

    bef = a.ops[befidx]     # aの最後
    aft = b[aftidx]         # bの最初

    aftend += 1     # +1の理由は後に出てくる

    while true
        # Operatorなので全てUpに正規化されている
        tb = typeof(bef)
        ta = typeof(aft)

        if tb !== ta
            # 型が異なる時
            break  # 中央で融合不可能なので終了 -> whileを抜ける
        elseif tb === OpTau
            # 両方OpTauの場合 -> mod2 = xor(a,b) & 1
            bef.cnt = xor(bef.cnt, aft.cnt) & 1
        elseif tb === OpDeriv
            # 両方OpDerivの場合 -> nが同じなら結合
            if bef.n == aft.n
                bef.cnt += aft.cnt
            end
        else
            # それ以外の場合は結合
            bef.cnt += aft.cnt
        end

        # もしbef.cnt==0なら消滅 → 再帰的処理へ
        if bef.cnt == 0
            # befは一つ減り aftは一つ後ろに行く
            befidx -= 1
            aftidx += 1

            if aftidx == aftend
                # aftidx += 1 してaftend(bnの最後のindexに+1したもの)と同じになる
                # -> aftidxは最後まで来たのでwhileを抜ける
                break
            end
            if befidx == 0
                # befidx -= 1 して0になる
                # -> befidxは最後まで来たのでwhileを抜ける
                break
            end
            # 何もないなら bef と aft を更新
            bef = a.ops[befidx]
            aft = b[aftidx]
            continue
        else
            # 結合成功かつ cnt != 0 (0にならなかったのでこれ以上は結合できない)
            # -> 右側先頭を削除し終了
            a.ops[befidx] = bef     #左側の最後は計算済みbef
            aftidx += 1             #右側の最初のindexを一つ後ろへ
            break
        end
    end

    # a.opsは 1:befidx まで(a.opsをわざわざpopしていない)
    resize!(a.ops,befidx)
    # 残りを連結
    # bはaftidxから最後まで(aftendが最後から1足されていることに注意)
    append!(a.ops,@view b[aftidx:aftend-1])
    return a
end

@inline function copy(op::Operator)::Operator
    r = Operator()
    if isone(op)
        return r
    end

    opl = lastindex(op.ops)
    b = Vector{AbstractOp}(undef,opl)   # 先に領域を確保してやると賢い
    for i in 1:opl
        x = op.ops[i]
        if x isa OpDeriv
            b[i] = OpDeriv(x.n,x.cnt)
        else
            b[i] = typeof(x)(x.cnt)
        end
    end
    r.ops = b
    return r
end
# これは正規化も兼ねている copy
@inline function copy(op::AbstractOp)::AbstractOp
    top = typeof(op)
    if top === OpDeriv
        return OpDeriv(op.n,op.cnt)
    elseif top === OpDown
        return OpUp(-op.cnt)
    elseif top === OpTau
        return OpTau(op.cnt & 1)
    else
        return (top)(op.cnt)
    end
end




"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Word compatible interface ##########################################################


# ===== Tuple 互換 =====
getindex(w::AbstractWord, i::Int) = w.t[i]
getindex(w::AbstractWord, r::UnitRange) = AbstractWord(w.t[r])
length(w::AbstractWord) = length(w.t)
iterate(w::AbstractWord, s...) = iterate(w.t, s...)
firstindex(w::AbstractWord) = firstindex(w.t)
lastindex(w::AbstractWord) = lastindex(w.t)
eltype(::Type{AbstractWord}) = Int
isempty(w::AbstractWord) = isempty(w.t)

# ===== 各種操作互換 =====
@inline copy(w::W) where {W<:AbstractWord} = W(w.t)                 # コピー
collect(w::AbstractWord) = collect(Int,w.t)                         # Vector化
Tuple(w::AbstractWord) = w.t                                        # タプル化
vcat(a::W, b::W) where {W<:AbstractWord} = W((a.t..., b.t...))      # 連結
reverse(w::W) where {W<:AbstractWord} = W(reverse(w.t))             # 反転
==(a::W, b::W) where {W<:AbstractWord} = a.t == b.t
hash(w::AbstractWord, h::UInt) = hash(w.t, h)
isless(a::W, b::W) where {W<:AbstractWord} = isless(a.t,b.t)
isequal(a::W, b::W) where {W<:AbstractWord} = isequal(a.t, b.t)
function convert(::Type{W}, t::Tuple{Vararg{Int}})::W where {W<:AbstractWord}
    return W(t)
end
# ===== スプラット互換 (a... が動くように) =====
in(item, w::AbstractWord) = in(item, w.t)
iterate(w::AbstractWord) = iterate(w.t)  # ←これが超重要！


###################################################################################################
############## Hoffman compatible interface ##########################################################


for fname in (:length, :isempty, :keys, :values, :pairs, :iterate)
    @eval begin
        $fname(x::Hoffman) = $fname(x.terms)
        $fname(x::Index)   = $fname(x.terms)
        $fname(x::Poly)    = $fname(x.terms)
    end
end
iterate(x::Hoffman,s...) = iterate(x.terms,s...)
iterate(x::Index,s...)   = iterate(x.terms,s...)
iterate(x::Poly,s...)    = iterate(x.terms,s...)
getindex(h::Hoffman, w::HoffmanWord) = get(h.terms, w, zero(Rational{BigInt}))
getindex(i::Index, w::IndexWord)     = get(i.terms, w, zero(Rational{BigInt}))
getindex(r::Poly{A}, d::Int) where A = get(r.terms, d, zero(A))
haskey(h::Hoffman, w::HoffmanWord) = haskey(h.terms, w)
haskey(i::Index, w::IndexWord)     = haskey(i.terms, w)
haskey(r::Poly, d::Int)            = haskey(r.terms, d)
==(a::Hoffman, b::Hoffman)         = a.terms == b.terms
==(a::Index, b::Index)             = a.terms == b.terms
==(a::Poly{A}, b::Poly{A}) where A = a.terms == b.terms
isequal(a::Hoffman, b::Hoffman)         = isequal(a.terms, b.terms)
isequal(a::Index, b::Index)             = isequal(a.terms, b.terms)
isequal(a::Poly{A}, b::Poly{A}) where A = isequal(a.terms, b.terms)
hash(h::Hoffman, u::UInt) = hash(h.terms, u)
hash(i::Index, u::UInt)   = hash(i.terms, u)
hash(r::Poly, u::UInt)    = hash(r.terms, u)
convert(::Type{Hoffman}, d::Dict{HoffmanWord, Rational{BigInt}})::Hoffman = Hoffman(d)
convert(::Type{Index}, d::Dict{IndexWord, Rational{BigInt}})::Index       = Index(d)
function convert(::Type{Poly{A}}, d::Dict{Int, A})::Poly{A} where A
    return Poly{A}(d)
end


###################################################################################################
############## Property Functions #################################################################


iszero(w::AbstractWord)::Bool     = false
iszero(x::Index)::Bool            = isempty(x.terms)
iszero(w::Hoffman)::Bool          = isempty(w.terms)
#iszero(hrm::HarmonicForm)::Bool   =  # TODO
#iszero(shf::ShuffleForm)::Bool    =  # TODO
#iszero(mpl::MPLCombination)::Bool =  # TODO
iszero(r::Poly)::Bool             = isempty(r.terms)

isone(w::AbstractWord)::Bool     = isempty(w.t)
isone(x::Hoffman)::Bool          = length(x) == 1 && haskey(x,HoffmanWord()) && isone(x[HoffmanWord()])
isone(x::Index)::Bool            = length(x) == 1 && haskey(x,IndexWord()) && isone(x[IndexWord()])
#isone(hrm::HarmonicForm)::Bool   =  # TODO
#isone(shf::ShuffleForm)::Bool    =  # TODO
#isone(mpl::MPLCombination)::Bool =  # TODO
isone(r::Poly)::Bool             = length(r) == 1 && haskey(r,0) && isone(r[0])

function isadmissible(idx::Index)::Bool
    if get_index_orientation()
        for w in keys(idx)
            if !isempty(w) && w[end] <= 1
                return false
            end
        end
        return true
    else
        for w in keys(idx)
            if !isempty(w) && w[1] <= 1
                return false
            end
        end
        return true
    end
end
function isadmissible(h::Hoffman)::Bool
    if get_index_orientation()
        for w in keys(h)
            if !isempty(w) && w[end] >= 2
                return false
            end
        end
        return true
    else
        for w in keys(h)
            if !isempty(w) && w[1] >= 2
                return false
            end
        end
        return true
    end
end

zero(::Type{Index})::Index               = Index()
zero(::Type{Hoffman})::Hoffman           = Hoffman()
zero(::Type{HarmonicForm})::HarmonicForm = HarmonicForm()
zero(::Type{ShuffleForm})::ShuffleForm   = ShuffleForm()
function zero(::Type{Poly{A}})::Poly{A} where A
    return Poly{A}()
end
one(::Type{HoffmanWord})::HoffmanWord = HoffmanWord()
one(::Type{IndexWord})::IndexWord     = IndexWord()
function one(::Type{Index})::Index
    idx = Index()
    idx[IndexWord()] = 1
    return  idx
end
function one(::Type{Hoffman})::Hoffman
    w = Hoffman()
    w[HoffmanWord()] = 1
    return  w
end
function one(::Type{Poly{A}})::Poly{A} where A
    r = Poly{A}()
    r[0] = one(A)
    return r
end
# TODO: one for hrm shf mpl

is_monomial(x::Index)::Bool                 = length(x.terms) == 1
is_monomial(w::Hoffman)::Bool               = length(w.terms) == 1
is_monomial(w::AbstractWord)::Bool          = true
#is_monomial(hrm::HarmonicForm)::Bool        = # TODO
#is_monomial(shf::ShuffleForm)::Bool         = # TODO
#is_monomial(mpl::MPLCombination)::Bool      = # TODO
is_monomial(r::Poly)::Bool                 = length(r.terms) == 1

is_hoffman(x::MPL)::Bool        = typeof(x) === Hoffman
is_index(x::MPL)::Bool          = typeof(x) === Index
is_monoindex(x::MPL)::Bool      = typeof(x) === MonoIndex
is_shuffleform(x::MPL)::Bool    = typeof(x) === ShuffleForm
is_harmonicform(x::MPL)::Bool   = typeof(x) === HarmonicForm
is_mplcombination(x::MPL)::Bool = typeof(x) === MPLCombination
is_shuffleexpr(x::MPL)::Bool    = typeof(x) <:  ShuffleExpr
is_harmonicexpr(x::MPL)::Bool   = typeof(x) <:  HarmonicExpr
is_zetaexpr(x::MPL)::Bool       = typeof(x) <:  ZetaExpr


###################################################################################################
############## Base Functions #####################################################################

# index compressor
""" example (2,2,1,1,2,1) -> [1,3,2] """
@inline function idxprs(w::HoffmanWord)::Vector{Int}
    n2 = count(==(2), w)
    v = Vector{Int}(undef,n2)
    idx = 0
    for x in w
        if x == 2
            idx += 1
            v[idx] = 1
        else
            v[idx] += 1
        end
    end
    return v
end
""" example (1,2,1,1,2,2) -> [2,3,1] """
@inline function idxprs_r(w::HoffmanWord)::Vector{Int}
    n2 = count(==(2), w)
    v = Vector{Int}(undef,n2)
    idx = 0
    cnt = 1
    for x in w
        if x == 1
            cnt += 1
        else
            idx += 1
            v[idx] = cnt
            cnt = 1
        end
    end
    return v
end
# index expander
""" example [1,3,2] -> (2,2,1,1,2,1) """
@inline function idxdprs(v::Vector{Int})::HoffmanWord
    # 総長さを一度に計算してから確保
    total_len = sum(v)
    w = ones(Int,total_len)

    pos = 1
    for t in v
        w[pos] = 2
        pos += t
    end
    return Word(w)
end
""" example [2,3,1] -> (1,2,1,1,2,2) """
@inline function idxdprs_r(v::Vector{Int})::HoffmanWord
    # 総長さを一度に計算してから確保
    total_len = sum(v)
    w = ones(Int,total_len)

    pos = 0
    for t in v
        pos += t
        w[pos] = 2
    end
    return Word(w)
end

@inline function copy(a::Hoffman)::Hoffman
    r = Hoffman()
    for (w,c) in a
        r[w] = c
    end
    return r
end
@inline function copy(a::Index)::Index
    r = Index()
    for (w,c) in a
        r[w] = c
    end
    return r
end
@inline function copy(p::Poly{A})::Poly{A} where A
    r = Poly{A}()
    for (deg,h) in p
        r[deg] = copy(h)
    end
    return r
end