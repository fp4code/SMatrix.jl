module SMatrix
import Base.convert

using LinearAlgebra

STerm = Union{Float32, Float64, Complex{Float32}, Complex{Float64},
              Array{Float32}, Array{Float64},
              Array{ Complex{Float32}}, Array{ Complex{Float64}}}

convert(::Type{STerm}, x::Int) = Float64(x)

mutable struct SMatI
    s11::STerm
    s12::STerm
    sa1::STerm
    sb2::STerm
    s21::STerm
    s22::STerm
    SMatI(s11,s12,sa1,sb2,s21,s22) = new(s11,s12,sa1,sb2,s21,s22)
    SMatI() = new()
end

mutable struct SMat
    s11::STerm
    s12::STerm
    s21::STerm
    s22::STerm
    SMat(s11,s12,s21,s22) = new(s11,s12,s21,s22)
    SMat(s::SMatI) = new(s.s11,s.s12,s.s21,s.s22)
    SMat() = new()
end

mutable struct SMatStack
    s11::Vector{STerm}
    s12::Vector{STerm}
    s21::Vector{STerm}
    s22::Vector{STerm}
    n::Int
    SMatStack(s::SMat) = new([s.s11],[s.s12],[s.s21],[s.s22],1)
    SMatStack(s11,s12,s21,s22) = new([s11],[s12],[s21],[s22],1)
    SMatStack(n::Int) = new(Vector{STerm}(n),Vector{STerm}(n),
                            Vector{STerm}(n),Vector{STerm}(n),n)
end

mutable struct Propagator
    p # plus direction
    m # minus direction
    Propagator(p::Number) = new(p,p)
    function Propagator(p::Vector)
        pd = Diagonal(p)
        new(pd,pd)
    end
    function Propagator(p::Vector,m::Vector)
        new(Diagonal(p),Diagonal(m))
    end
end

struct Stack
    s::Vector{Union{SMat,Array{SMat}}}
    p::Vector{Union{Propagator,Array{Propagator}}}
#    function Stack(s::Vector{Union{SMat,Array{SMat}}}, p::Vector{Union{Propagator,Array{Propagator}}})
    function Stack(s,p)
        if length(s) != length(p)+1
            throw(ArgumentError("size mismatch"))
        end
        new(s,p)
    end
    function Stack(stack::Vector{Any})
        s = stack[1:2:end]
        p = stack[2:2:end]
        if length(s) != length(p)+1
            throw(ArgumentError("stack should not begin or end with a P term"))
        end
        new(s,p)
    end
end

#=
useful to compute fields

1       s11
--------------
sa1[1]  ^^^^^^
vvvvvv  sb1[1]
--------------
sa1[2]  ^^^^^^
vvvvvv  sb1[2]
--------------
...
--------------
sa1[end]   ^^^
vvv   sb1[end]
--------------
s22

=#
mutable struct AmpStack
    s11::STerm
    sa1::Array{STerm}
    sb1::Array{STerm}
    s21::STerm
    AmpStack() = new()
    AmpStack(s11,sa1,ab1,s21) = new(s11,sa1,ab1,s21)
end

#function isgoodstask(stack::Stack)
#    TO BE DONE
#end

function add_layer(s_XY::SMat,pY::Propagator,s_YZ::SMat)
    # TODO: use code from add_layer_i
    ps11 = pY.m * s_YZ.s11
    ps12 = pY.m * s_YZ.s12
    ps22 = pY.p * s_XY.s22
    ps21 = pY.p * s_XY.s21
    # symmetric formulation, maybe not efficient...
    loop12 = ps11*ps22
    loop21 = ps22*ps11
    d1 = s_XY.s12 / (one(loop12) - loop12)
    d2 = s_YZ.s21 / (one(loop21) - loop21)
    #
    s11 = s_XY.s11 + d1 * ps11 * ps21
    s12 = d1 * ps12
    s22 = s_YZ.s22 + d2 * ps22 * ps12
    s21 = d2 * ps21
    #
    return SMat(s11,s12,s21,s22)
end

function add_layer_i(s_XY::SMat,pY::Propagator,s_YZ::SMat)
    a22pm = s_XY.s22 * pY.m
    b11pp = s_YZ.s11 * pY.p
    loopaba = a22pm*b11pp
    loopbab = b11pp*a22pm
    sa1 = (one(loopaba) - loopaba)\s_XY.s21
    sb2 = (one(loopbab) - loopbab)\s_YZ.s12
    ppsa1 = pY.p*sa1
    pmsb2 = pY.m*sb2
    s12 = s_XY.s12*pmsb2
    s21 = s_YZ.s21*ppsa1
    s11 = s_XY.s11 + s_XY.s12*pY.m*s_YZ.s11*ppsa1
    s22 = s_YZ.s22 + s_YZ.s21*pY.p*s_XY.s22*pmsb2
    return SMatI(s11,s12,sa1,sb2,s21,s22)
end

function internal_layer(s_XY::SMat,pY::Propagator,s_YZ::SMat)
    a22pm = s_XY.s22 * pY.m
    b11pp = s_YZ.s11 * pY.p
    loopaba = a22pm*b11pp
    sa = (one(loopaba) - loopaba)\s_XY.s21
    sb = pY.m * s_XY.s11 * pY.p * sa
    return (sa, sb)
end

function add_layer(a::SMatStack,p::Propagator,b::SMatStack)
    a11 = a.s11[1]
    a21 = a.s21[end]
    a22 = a.s22[end]
    a12 = a.s12[1]
    #
    b11 = b.s11[1]
    b21 = b.s21[end]
    b22 = b.s22[end]
    b12 = b.s12[1]
    #
    loopaba = a22*p.m*b11*p.p
    loopbab = b11*p.p*a22*p.m
    #
    a21l = (one(loopaba) - loopaba)\a21
    b12l = (one(loopbab) - loopbab)\b12
    #
    c = SMatStack(a.n + 1 + b.n)
    throw(AssertionError("Unfinished code"))
end


function forceone(i::Int, one::Bool)
    if one
        1
    else
        i
    end
end
    
function add_layer(a::Union{Array{SMat},Array{SMatStack}},
                   p::Array{Propagator},
                   b::Union{Array{SMat},Array{SMatStack}})
    s = max(size(a), size(p), size(b))
    if size(a) == (1,)
        ia = true
    elseif size(a) == s
        ia = false
    else
        throw(ArgumentError("Size mismatch for argument 1"))
    end
    if size(p) == (1,)
        ip = true
    elseif size(p) == s
        ip = false
    else
        throw(ArgumentError("Size mismatch for argument 2"))
    end
    if size(b) == (1,)
        ib = true
    elseif size(b) == s
        ib = false
    else
        throw(ArgumentError("Size mismatch for argument 3"))
    end
    sm = Vector{SMatrix.SMat}(undef, s)
    for i in 1:prod(s)
        sm[i] = add_layer(a[forceone(i, ia)],
                          p[forceone(i, ip)],
                          b[forceone(i, ib)])
    end
    return sm
end

add_layer(a::Array{SMat},p::Array{Propagator},b::SMat) = add_layer(a,p,[b])
add_layer(a::Array{SMat},p::Propagator,b::Array{SMat}) = add_layer(a,[p],b)
add_layer(a::Array{SMat},p::Propagator,b::SMat) = add_layer(a,[p],[b])
add_layer(a::SMat,p::Array{Propagator},b::Array{SMat}) = add_layer([a],p,b)
add_layer(a::SMat,p::Array{Propagator},b::SMat) = add_layer([a],p,[b])
add_layer(a::SMat,p::Propagator,b::Array{SMat}) = add_layer([a],[p],b)

function compute_stack_p(stack::Stack)
    n = length(stack.p)
    if length(stack.s) != n+1
        throw(ArgumentError("size mismatch"))
    end
    s = stack.s[1]
    for i in 1:n
        s = add_layer(s, stack.p[i], stack.s[i+1])
    end
    return s
end

function compute_stack_m(stack::Stack)
    n = length(stack.p)
    if length(stack.s) != n+1
        throw(ArgumentError("size mismatch"))
    end
    s = stack.s[n+1]
    for i in n:-1:1
        s = add_layer(stack.s[i], stack.p[i], s)
    end
    return s
end

function compute_stack_full(stack::Stack)
    n = length(stack.p)
    #
    sa21 = Vector{STerm}(undef, n+1)
    sa22 = Vector{STerm}(undef, n+1)
    s = stack.s[1]
    sa21[1] = s.s21
    sa22[1] = s.s22
    for i in 1:n
        s = add_layer(s, stack.p[i], stack.s[i+1])
        sa21[i+1] = s.s21
        sa22[i+1] = s.s22
    end
    #
    s11 = s.s11
    s21 = s.s21
    sb12 = Vector{STerm}(undef, n+1)
    sb11 = Vector{STerm}(undef, n+1)
    s = stack.s[n+1]
    sb12[n+1] = s.s12
    sb11[n+1] = s.s11
    for i in n:-1:1
        s = add_layer(stack.s[i], stack.p[i], s)
        sb12[i] = s.s12
        sb11[i] = s.s11
    end
    #
    sa1 = Vector{STerm}(undef, n)
    sb1 = Vector{STerm}(undef, n)
    for i in 1:n
        a22pm = sa22[i] * stack.p[i].m
        b11pp = sb11[i+1] * stack.p[i].p
        loopaba = a22pm*b11pp
        sa1[i] = (one(loopaba) - loopaba)\sa21[i]
        sb1[i] = b11pp * sa1[i]
    end
    # return (AmpStack(s11,sa1,sb1,s21),sa21,sa22,sb11,sb12)
    return AmpStack(s11,sa1,sb1,s21)
end


end # module
