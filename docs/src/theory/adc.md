# Apparent diffusion coefficient

Consider a unidirectional gradient ``\vec{g}(t) = g f(t) \vec{d}`` parametrized by an
amplitude ``g > 0``, time profile ``f : [0, T_\text{echo}] \to [-1, 1]`` and direction
``\vec{d}`` with ``\| \vec{d} \| = 1``.

Some commonly used time profiles (diffusion-encoding sequences) are:

- The pulsed-gradient spin echo (PGSE) sequence:
  ```math
  f(t) =
  \begin{cases}
  1, \quad & t_1 \leq t \leq t_1+\delta, \\
  -1, \quad & t_1 + \Delta < t \leq t_1 + \Delta + \delta, \\
  0, \quad & \text{otherwise};
  \end{cases}
  ```
- The double pulsed-gradient spin echo (DPGSE) sequence:
  ```math
  f(t) = f_{\delta, \Delta}(t) + f_{\delta, \Delta}\left(t - (\Delta + \delta +
  t_\text{pause}) \right),
  ```
  where ``f_{\delta, \Delta}`` is a normal PGSE sequence;
- The oscillating gradient spin echo (OGSE) sequence (here, cos-OGSE):
  ```math
  f(t) =
  \begin{cases}
      \cos\left(n \frac{2\pi}{\delta} (t - t_1)\right), \quad & t_1 < t \leq t_1
      + \delta, \\
      -\cos\left(n \frac{2 \pi}{\delta} (t - \Delta - t_1)\right), \quad & t_1 + \Delta < t
      \leq t_1 + \Delta + \delta, \\
      0, \quad & \text{otherwise},
  \end{cases}
  ```

where ``t_1 \geq 0`` and ``t_1 + \Delta \geq T_\text{echo} / 2``.

## B-value

In a DMRI experiment, ``f`` is usually fixed, while ``g`` or ``\vec{d}`` are varied. ``S``
is usually plotted against a quantity called the ``b``-value. The ``b``-value depends on
``f`` and ``g``:

```math
b(g, f) = \gamma^2 g^2 \int_0^{T_\text{echo}} \left( \int_0^t f(s) \, \mathrm{d} s \right)^2
\, \mathrm{d}t.
```

For PGSE, the b-value is:

```math
b(g, \delta, \Delta) = \gamma^2 g^2 \delta^2 \left( \Delta - \delta / 3 \right).
```

For the cosine OGSE with _integer_ number of periods ``n`` in each of the two
durations ``\delta``, the corresponding ``b``-value is:

```math
b(g, \delta) = \gamma^2 g^2 \frac{\delta^3}{4 n^2 \pi^2} = \gamma^2 g^2
\frac{\delta}{\omega^2}.
```

The reason for these definitions is that in a homogeneous medium, the signal attenuation is
``\exp(-\vec{d}^\mathsf{T} \mathbf{D}_0 \vec{d} b)``, where ``\mathbf{D}_0`` is the
intrinsic diffusion tensor.

## ADC

An important quantity that can be derived from the dMRI signal is the "Apparent Diffusion
Coefficient" (ADC), which gives an indication of the root mean squared distance travelled by
water molecules in the gradient direction ``\vec{d}``, averaged over all starting positions:

```math
D_\text{ADC} = \left. -\frac{\partial}{\partial b} \log{\frac{S(b)}{S(0)}}\right\vert_{b=0}.
```

## Fitting the ADC from the dMRI signal

We numerically compute ADC by a polynomial fit of

```math
\log\frac{S(b)}{S_\text{initial}} \approx c_0+c_1b+\dots+c_n b^n,
```

increasing ``n`` from 1 onwards until we get the value of ``c_1`` to be stable within a
numerical tolerance, where ``S_\text{initial} = \int_\Omega \rho(\vec{x}) \, \mathrm{d}
\Omega(\vec{x})``. The first coefficient is given by ``c_0 = \log
\frac{S(0)}{S_\text{initial}}``, which is equal to zero if ``T_{2,i} = \infty`` for all
``i``.

The ADC may be interpreted as a correction to the intrinsic diffusion coefficient, taking
into account the deviation from the free diffusion arising from interior interfaces,
non-isotropic diffusion, fluid movements etc. For suffiently small b-values, the signal
attenuation is given by

```math
\mathrm{e}^{-D_\text{ADC} b}.
```

Another quantity of interest is a slight generalization of the ADC -- an effective diffusion
tensor ``\mathbf{D}_\text{eff}``. The six coefficients of this symmetric positive tensor is
fitted to best approximate the following signal attenuation, for all directions ``\vec{d}``:

```math
\mathrm{e}^{-\vec{g}^\mathsf{T} \mathbf{D}_\text{eff} \vec{d} b}.
```

The resulting ADC in direction ``\vec{d}`` is then given by

```math
D_\text{ADC}(\vec{d}) = \vec{d}^\mathsf{T} \mathbf{D}_\text{eff} \vec{d}.
```

## HADC

In a previous work, a PDE model for the time-dependent ADC was obtained
starting from the Bloch-Torrey equation, using homogenization techniques. In the case of
negligible water exchange between compartments (low permeability), there is no coupling
between the compartments, at least to the quadratic order in ``g``, which is the ADC term.
The ADC in compartment ``\Omega`` is given by

```math
D_\text{HADC}(\vec{d}, f) = \vec{d}^\mathsf{T} \mathbf{D} \vec{d} -
\frac{\int_0^{T_\text{echo}} F(t) h(t) \, \mathrm{d}t}{\int_0^{T_\text{echo}} F(t)^2 \,
\mathrm{d} t},
```

where ``F(t) = \int_0^t f(s) \, \mathrm{d}s,`` and

```math
h(t) = \frac{1}{|\Omega|} \int_{\partial \Omega} \omega(\vec{x}, t) \, \vec{d} \cdot
\vec{n}(\vec{x}) \, \mathrm{d} \Gamma(\vec{x})
```

is a quantity related to the directional gradient of a function ``\omega`` that is the
solution of the homogeneous diffusion equation with Neumann boundary condition and zero
initial condition:

```math
\begin{alignedat}{2}
    \frac{\partial}{\partial t} \omega(\vec{x}, t) & = \nabla \cdot \left(\mathbf{D} \nabla
    \omega(\vec{x}, t)\right), \quad & & (\vec{x}, t) \in \Omega \times [0, T_\text{echo}],
    \\
    \omega(\vec{x}, 0) & = 0, \quad & & \vec{x} \in \Omega, \\
    \mathbf{D} \nabla \omega(\vec{x}, t) \cdot \vec{n} & = F(t) \mathbf{D} \vec{d} \cdot
    \vec{n}, \quad &  & \vec{x} \in \partial \Omega,
\end{alignedat}
```

``\vec{n}`` being the outward normal. The above set of equations comprise the homogenized
model that we call the HADC model.



## Short diffusion time approximation

A well-known formula for the ADC in the short diffusion time regime is the following short
time approximation (STA):

```math
D_\text{STA} = \left(1 - \frac{4\sqrt{D_0}}{3 \sqrt{\pi}}\sqrt{\Delta} \frac{{A}}{d\; V}
\right) D_0,
```
where ``\dfrac{A}{V}`` is the surface to volume ratio, ``d`` is the spatial dimension and
``D_0 = \vec{d}^\mathsf{T} \mathbf{D} \vec{d}`` is the intrinsic diffusion
coefficient in the gradient direction. In the above formula, the pulse duration ``\delta``
is assumed to be very small compared to ``\Delta``. A recent correction to the above formula
, taking into account the finite pulse duration ``\delta`` and the
gradient direction ``\vec{d}``, is the following:

```math
D_\text{STA} = \left(1 - \frac{4 \sqrt{D_0}}{3 \sqrt{\pi}} C_{\delta,\Delta}
\frac{A_{\vec{d}}}{V} \right) D_0,
```

where

```math
A_{\vec{d}} = \int_{\partial \Omega} \left(\vec{d} \cdot \vec{n}(\vec{x})\right)^2 \,
\mathrm{d} \Gamma(\vec{x}) = \vec{d}^\mathsf{T} \mathbf{N} \vec{d}
```

with ``\mathbf{N} = \int_{\partial \Omega} \vec{n}(\vec{x}) \vec{n}^\mathsf{T}\!(\vec{x}) \,
\mathrm{d} \Gamma(\vec{x})`` and

```math
C_{\delta, \Delta} = \dfrac{4}{35} \dfrac{\left(\Delta + \delta \right)^{7 / 2}+
\left(\Delta - \delta\right)^{7 / 2} - 2 \left(\delta^{ 7 / 2} + \Delta^{7 / 2}
\right)}{\delta^2 \left( \Delta - \delta / 3\right)} = \sqrt{\Delta} \left(1 + \dfrac{1}{3}
\dfrac{\delta}{\Delta} - \dfrac{8}{35} \left(\dfrac{\delta}{\Delta}\right)^{3 / 2} + \dots
\right).
```
When ``\delta \ll \Delta``, the value ``C_{\delta, \Delta}`` is approximately
``\sqrt{\Delta}``.
