# SMatrix

Non commutative S-Matrix code, mainly for modal analysis of optical periodic structures.

## Installation

```
import Pkg;
Pkg.activate("mytestenv")
Pkg.add(url="git://github.com/fp4code/SMatrix.jl")
```

## Testing

```
Pkg.test("SMatrix")
```

## Usages

### scalar computation: anti-reflection coating glass

```
using SMatrix

na = 1.0                                   # air index
nb = sqrt(1.5)                             # AR index
nc = 1.5                                   # glass index
ra = (na-nb)/(na+nb)                       # amplitude of reflection air/AR
rb = (nb-nc)/(nb+nc)                       # amplitude of reflection AR/glass
s_ab = SMatrix.SMat(ra, 1-ra, 1+ra, -ra)   # air/AR S-matrix, 4 numbers
s_bc = SMatrix.SMat(rb, 1-rb, 1+rb, -rb)   # AR/glass S-matrix, 4 numbers
p_b  = SMatrix.Propagator(exp(1im*2*pi/4)) # AR propagator (1im indeed)
s_ac = SMatrix.add_layer(s_ab, p_b, s_bc)  # air/AR/glass S-matrix
R = abs2(s_ac.s11)                        # reflection is 0
T = abs2(s_ac.s21)*nc/na                  # transmission is 1
```

### add_layer arguments can be vectors (all having the same size, or 1)

```
h = 0.6/(4*nb)                             # AR thickness
vwl = 0.4:0.01:0.8                         # wavelengths
p_b  = map(x -> SMatrix.Propagator(exp(1im*nb*2*pi*h/x)),
       	   vwl)                            # array of propagators
s_ac = SMatrix.add_layer(s_ab, p_b, s_bc)  # s_ab and a_bc are fixed in this example
R = map(x->abs2(x.s11), s_ac)             # array of reflectivities
```

Now we can plot the air/AR/glass relectivity as a function of the wavelength:

```
# Pkg.add("PyPlot")
using PyPlot
plot(vwl, R)
```

### diffraction S-matrix combination

Here there is one mode in first medium, two modes in the layer and one mode in third medium.
Theses values are arbitrary, and certainly not reciprocal!

```
a = SMatrix.SMat([0.2], [0.3 0.4], [0.5; 0.6], [0.3 -0.4; 0.1 -0.03])
b = SMatrix.SMat([0.31 -0.42; 0.13 -0.031im], [0.51; 0.63], [0.34 0.41], [0.22])
p = SMatrix.Propagator([0.8+0.1im, 0.01])
ab = SMatrix.add_layer(a, p, b)              # 4 numbers ab.s11 ab.s12 ab.s21 ab.s22
```

### multilayered mirror

Not that compute_stack_p or compute_stack_m are not the more efficient ways
to compute an ababababab... multilayer.

```
using SMatrix

n0 = 1
na = 1.5
nb = 2
wl0 = 0.5
ha = wl0/(4*na)
hb = wl0/(4*nb)
vwl = 0.4:0.00001:0.8                        # wavelengths
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

using PyPlot
plot(vwl, Rp)
```

Using two such mirrors to create a Fabry-Perot etalon:

```
sa = sp
sb = map(x -> SMatrix.SMat(x.s22, x.s21, x.s12, x.s11), sa)

h = 10*wl0/(2*na)
p  = map(x -> SMatrix.Propagator(exp(1im*na*2*pi*h/x)),
         vwl)                            # array of propagators
s = SMatrix.compute_stack_p(SMatrix.Stack(Any[sa,p,sb]))
T = map(x->abs2(x.s21), s)
plot(vwl, T)
```

### field in a multilayered mirror

```
using SMatrix

n0 = 1
na = 1.5
nb = 2
wl0 = 0.5
ha = wl0/(4*na)
hb = wl0/(4*nb)
wl = wl0

r0a0 = (n0-na)./(n0+na)
rbab = (nb-na)./(nb+na)
s0a = SMatrix.SMat(r0a0, 1-r0a0, 1+r0a0, -r0a0)
sa0 = SMatrix.SMat(-r0a0, 1+r0a0, 1-r0a0, +r0a0)
sba = SMatrix.SMat(rbab, 1-rbab, 1+rbab, -rbab)
sab = SMatrix.SMat(-rbab, 1+rbab, 1-rbab, +rbab)

vi = [s0a,sab,sba,sab,sba,sab,sba,sab,sba,sab,sba,sab,sba,sab,sba,sab,sba];
vh = [ha,hb,ha,hb,ha,hb,ha,hb,ha,hb,ha,hb,ha,hb,ha,hb];
vn = [na,nb,na,nb,na,nb,na,nb,na,nb,na,nb,na,nb,na,nb];
vp = map(x -> SMatrix.Propagator(x), exp.(1im*vn*2*pi.*vh/wl));

stack = SMatrix.Stack(vi, vp)

sp = SMatrix.compute_stack_p(stack);
R = abs2(sp.s11)
si = SMatrix.compute_stack_full(stack);

vhs = cumsum([0;vh])
x = Vector{Float64}()
a = Vector{Complex{Float64}}()
for i in 1:length(vh)
    vdx = range(0, vh[i], length=11)[1:end-1]
    for xx in vdx
    	push!(x, vhs[i] + xx)
        pa  = exp(1im*vn[i]*2*pi*xx/wl)
        pb  = exp(1im*vn[i]*2*pi*(vh[i]-xx)/wl)
	a1 = pa*si.sa1[i]
	a2 = pb*si.sb1[i]
	push!(a, a1+a2)
    end
end

using PyPlot
plot(x, real(a),x, imag(a))
```
