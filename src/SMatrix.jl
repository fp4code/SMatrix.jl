module SMatrix

type SMat
    s11
    s12
    s21
    s22
    SMat(s11,s12,s21,s22) = new(s11,s12,s21,s22)
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

function add_layer(a::Array{SMat},p::Array{Propagator},b::Array{SMat})
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

end # module
