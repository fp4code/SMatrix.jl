# SMatrix

Non commutative S-Matrix code, mainly for modal analysis of optical periodic structures.

## Installation

```
Pkg.clone("git://github.com/MiNaO/SMatrix.jl")
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
R = abs(s_ac.s11)^2                        # reflection is 0
T = abs(s_ac.s21)^2*nc/na                  # transmission is 1
```

### add_layer arguments can be vectors (all having the same size, or 1)

```
h = 0.6/(4*nb)                             # AR thickness
vwl = 0.4:0.01:0.8                         # wavelengths
p_b  = map(x -> SMatrix.Propagator(exp(1im*nb*2*pi*h/x)),
       	   vwl)                            # array of propagators
s_ac = SMatrix.add_layer(s_ab, p_b, s_bc)  # s_ab and a_bc are fixed in this example
R = map(x->abs(x.s11)^2, s_ac)             # array of reflectivities
```

Now we can plot the air/AR/glass relectivity as a fucntion of the wavelength:

```
using PyPlot
plot(vwl, R)
```

### Diffraction S-matrix combination

Here there is one mode in first medium, two modes in the layer and one mode in third medium.
Theses values are arbitrary, and certainly not reciprocal!

```
a = SMatrix.SMat(0.2,[0.3 0.4],[0.5 0.6].',[0.3 -0.4;0.1 -0.03])
b = SMatrix.SMat([0.31 -0.42;0.13 -0.031im],[0.51 0.63].',[0.34 0.41],0.22)
p = SMatrix.Propagator([0.8+0.1im, 0.01])
ab = SMatrix.add_layer(a,p,b)              # 4 numbers
```
