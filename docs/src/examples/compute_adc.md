# Compare ADCs

SpinDoctor comes with multiple approaches for computing the apparent diffusion coefficient
(ADC)
for a `ScalarGradient` $\vec{g}(t) = f(t) g \vec{d}$:

- The free diffusion coefficient $\frac{\vec{d}' D \vec{d}}{\vec{d}' \vec{d}}$, which
    represents unrestricted diffusion in the absence of boundaries
- Computing the short diffusion time approximation for the ADC
- Fitting the signal obtained by solving the BTPDE for different $b$-values
- Solving a homogenized model (HADC)
- Using the matrix formalism effective diffusion tensor

In this example we will compare the different approaches for a mesh.

```julia
using SpinDoctor
```
