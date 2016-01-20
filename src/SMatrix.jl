module SMatrix
import Base.convert

STerm = Union{Float32, Float64, Complex64, Complex128,
              Array{Float32}, Array{Float64},
              Array{Complex64}, Array{Complex128}}

convert(::Type{STerm}, x::Int) = Float64(x)

type SMatI
    s11::STerm
    s12::STerm
    sa1::STerm
    sb2::STerm
    s21::STerm
    s22::STerm
    SMatI(s11,s12,sa1,sb2,s21,s22) = new(s11,s12,sa1,sb2,s21,s22)
end

type SMat
    s11::STerm
    s12::STerm
    s21::STerm
    s22::STerm
    SMat(s11,s12,s21,s22) = new(s11,s12,s21,s22)
    SMat(s::SMatI) = new(s.s11,s.s12,s.s21,s.s22)
end

type SMatStack
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

type Propagator
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

Stack = Vector{Union{SMat,Propagator,Array{SMat},Array{Propagator}}}

#function isgoodstask(stack::Stack)
#    TO BE DONE
#end

eye(x::Number) = one(x)
eye(x) = Base.eye(x)

function add_layer(s_XY::SMat,pY::Propagator,s_YZ::SMat)
    ps11 = pY.m * s_YZ.s11
    ps12 = pY.m * s_YZ.s12
    ps22 = pY.p * s_XY.s22
    ps21 = pY.p * s_XY.s21
    # symmetric formulation, maybe not efficient...
    loop12 = ps11*ps22
    loop21 = ps22*ps11
    d1 = s_XY.s12 / (eye(loop12) - loop12)
    d2 = s_YZ.s21 / (eye(loop21) - loop21)
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
    sa1 = (eye(loopaba) - loopaba)\s_XY.s21
    sb2 = (eye(loopbab) - loopbab)\s_YZ.s12
    ppsa1 = pY.p*sa1
    pmsb2 = pY.m*sb2
    s12 = s_XY.s12*pmsb2
    s21 = s_YZ.s21*ppsa1
    s11 = s_XY.s11 + s_XY.s12*pY.m*s_YZ.s11*ppsa1
    s22 = s_YZ.s22 + s_YZ.s21*pY.p*s_XY.s22*pmsb2
    return SMatI(s11,s12,sa1,sb2,s21,s22)
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
    a21l = (eye(loopaba) - loopaba)\a21
    b12l = (eye(loopbab) - loopbab)\b12
    #
    c = SMatStack(a.n + 1 + b.n)
    throw(AssertionError("Unfinished code"))
end

    
function add_layer(a::Union{Array{SMat},Array{SMatStack}},
                   p::Array{Propagator},
                   b::Union{Array{SMat},Array{SMatStack}})
    s = max(size(a), size(p), size(b))
    if size(a) == (1,)
        ia(i::Int) = 1
    elseif size(a) == s
        ia(i::Int) = i
    else
        throw(ArgumentError("Size mismatch for argument 1"))
    end
    if size(p) == (1,)
        ip(i::Int) = 1
    elseif size(p) == s
        ip(i::Int) = i
    else
        throw(ArgumentError("Size mismatch for argument 2"))
    end
    if size(b) == (1,)
        ib(i::Int) = 1
    elseif size(b) == s
        ib(i::Int) = i
    else
        throw(ArgumentError("Size mismatch for argument 3"))
    end
    sm = Array(SMatrix.SMat, s)
    for i in 1:prod(s)
        sm[i] = add_layer(a[ia(i)],p[ip(i)],b[ib(i)])
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
    s = stack[1]
    for i in 2:2:length(stack)
        # TODO: add checks
        s = add_layer(s, stack[i], stack[i+1])
    end
    return s
end

function compute_stack_m(stack::Stack)
    s = stack[end]
    for i in length(stack)-1:-2:1
        # TODO: add checks
        s = add_layer(stack[i-1], stack[i], s)
    end
    return s
end



end # module
