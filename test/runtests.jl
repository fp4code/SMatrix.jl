using SMatrix
using Base.Test

println("stupid scalar test")
sxy = SMatrix.SMat(0,1,1,0)
syz = SMatrix.SMat(0,1,1,0)
py = SMatrix.Propagator(1)
sxz = SMatrix.add_layer(sxy,py,syz)

@test sxy.s11 ≈ 0
@test sxy.s12 ≈ 1
@test sxy.s21 ≈ 1
@test sxy.s22 ≈ 0

println("anti-reflection coating on glass")
n1 = 1.0
n3 = 1.5
n2 = sqrt(n3)
t21 = 2*n1/(n1+n2)
t12 = 2*n2/(n1+n2)
t32 = 2*n2/(n2+n3)
t23 = 2*n3/(n2+n3)
r121 = t21-1
r212 = t12-1
r232 = t32-1
r323 = t23-1

sxy = SMatrix.SMat(r121,t12,t21,r212)
syz = SMatrix.SMat(r232,t23,t32,r323)
py = SMatrix.Propagator(1im)
sxz = SMatrix.add_layer(sxy,py,syz)

@test_approx_eq_eps sxz.s11 0.0 1e-14
@test_approx_eq_eps sxz.s11 0.0 1e-14
@test_approx_eq_eps abs(sxz.s21) sqrt(n1/n3) 1e-14
@test_approx_eq_eps abs(sxz.s12) sqrt(n3/n1) 1e-14

sxz = SMatrix.add_layer_i(sxy,py,syz)

@test_approx_eq_eps sxz.s11 0.0 1e-14
@test_approx_eq_eps sxz.s11 0.0 1e-14
@test_approx_eq_eps abs(sxz.s21) sqrt(n1/n3) 1e-14
@test_approx_eq_eps abs(sxz.s12) sqrt(n3/n1) 1e-14

#=
sxy = SMatrix.SMatStack(r121,t12,t21,r212)
syz = SMatrix.SMatStack(r232,t23,t32,r323)
py = SMatrix.Propagator(1im)
sxz = SMatrix.add_layer(sxy,py,syz)
=#

println("anti-reflection coating on glass, combination of scalar/array arguments")
na = 1.0
vnb = sqrt([1.4,1.5,1.6])
vnc = [1.4,1.5,1.6]
vwl = [0.4,0.5,0.6]
h = vwl[2]/(4*vnb[2])

scalarize(x::Array) = (length(x) == 1 ? x[1] : x)
scalarize(x) = x

vectorize(x::Array) = x
vectorize(x) = [x]

for nb in (vnb[2], vnb)
    for nc in (vnc[2], vnc)
        for wl in (vwl[2], vwl)
            ra = (na-nb)./(na+nb)
            rb = (nb-nc)./(nb+nc)
            sab = scalarize(map(x -> SMatrix.SMat(x, 1-x, 1+x, -x), ra))
            sbc = scalarize(map(x -> SMatrix.SMat(x, 1-x, 1+x, -x), rb))
            phases = exp(1im*2*pi*h*(nb./wl))
            pb  = scalarize(map(x -> SMatrix.Propagator(x), phases))
            sac = SMatrix.add_layer(sab, pb, sbc)
            R = scalarize(map(x -> abs(x.s11)^2, vectorize(sac)))
            if issubtype(typeof(R), Number)
                @test R<1e-10
            else
                @test R[1]>1e-5 && R[2]<1e-10 && R[3]>1e-5
            end
        end
    end
end

println("arbitrary 121 system")
a = SMatrix.SMat(0.2,[0.3 0.4],[0.5 0.6].',[0.3 -0.4;0.1 -0.03])
b = SMatrix.SMat([0.31 -0.42;0.13 -0.031im],[0.51 0.63].',[0.34 0.41],0.22)
pp = [0.1+1im,0.4]
pm = [0.2+1im,0.2]
p = SMatrix.Propagator(pp,pm)
ab = SMatrix.add_layer(a,p,b)
abi = SMatrix.add_layer_i(a,p,b)

pp = diagm(pp)
pm = diagm(pm)
@test ab.s11 ≈  a.s11 + a.s12*inv(eye(2) - pm*b.s11*pp*a.s22)*pm*b.s11*pp*a.s21
@test ab.s22 ≈  b.s22 + b.s21*inv(eye(2) - pp*a.s22*pm*b.s11)*pp*a.s22*pm*b.s12
@test ab.s21 ≈ b.s21*inv(eye(2) - pp*a.s22*pm*b.s11)*pp*a.s21
@test ab.s12 ≈ a.s12*inv(eye(2) - pm*b.s11*pp*a.s22)*pm*b.s12

@test ab.s11 ≈ abi.s11
@test ab.s22 ≈ abi.s22
@test ab.s21 ≈ abi.s21
@test ab.s12 ≈ abi.s12

println("multilayered mirror")

n0 = 1
na = 1.5
nb = 2
wl0 = 0.5
ha = wl0/(4*na)
hb = wl0/(4*nb)
vwl = 0.4:0.001:0.8                        # wavelengths
pa  = map(x -> SMatrix.Propagator(exp(1im*na*2*pi*ha/x)),
          vwl)                            # array of propagators
pb  = map(x -> SMatrix.Propagator(exp(1im*nb*2*pi*hb/x)),
          vwl)                            # array of propagators

r0a0 = (n0-na)./(n0+na)
rbab = (nb-na)./(nb+na)
s0a = SMatrix.SMat(r0a0, 1-r0a0, 1+r0a0, -r0a0)
sa0 = SMatrix.SMat(-r0a0, 1+r0a0, 1-r0a0, +r0a0)
sba = SMatrix.SMat(rbab, 1-rbab, 1+rbab, -rbab)
sab = SMatrix.SMat(-rbab, 1+rbab, 1-rbab, +rbab)
stack = SMatrix.Stack(Any[s0a,
                          pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba,
                          pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba,
                          pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba,pa,sab,pb,sba])
sp = SMatrix.compute_stack_p(stack)
sm = SMatrix.compute_stack_m(stack)
Rp = map(x->abs2(x.s11), sp)
Rm = map(x->abs2(x.s11), sm)
@test maximum(Rp - Rm) < 1e-14




