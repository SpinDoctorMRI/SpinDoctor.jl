# Matrix formalism

# Eigenvalue decomposition of the generalized Laplace operator

Let ``\{(\phi, \lambda)\}`` be the ``L^2``-normalized eigenfunctions and eigenvalues of the
generalized Laplace operator ``\nabla \cdot \mathbf{D} \nabla``, defined almost everywhere on
the domain ``\Omega = \bigcup_{i = 1}^{N_\text{cmpt}} \Omega_i``. Denoting ``\phi^i`` the
restriction of ``\phi`` to compartment ``\Omega_i``, the eigenfunctions respect the
following set of equations, for ``(i, j) \in \{1, \dots, N_\text{cmpt}\}^2``:

```math
-\nabla \cdot \mathbf{D}_i \nabla \phi^i(\vec{x}) = \lambda \phi^i(\vec{x}), \quad \vec{x} \in
\Omega_i, \quad i \in \{1, \dots, N_\text{cmpt}\}.
```

where the boundary conditions are the same as for the BTPDE:

```math
\begin{alignedat}{3}
    \mathbf{D}_i \nabla \phi^i(\vec{x}) \cdot \vec{n}_i(\vec{x}) & = -\mathbf{D}_j \nabla
    \phi^j(\vec{x}) \cdot \vec{n}_j(\vec{x}), \quad & & \vec{x} \in \Gamma_{i j}, \quad (i,
    j) \in \{1, \dots, N_\text{cmpt}\}^2, & \\
    \mathbf{D}_i \nabla \phi^i(\vec{x}) \cdot \vec{n}_i(\vec{x}) & = \kappa_{i j} \left(c_{i j}
    \phi^j(\vec{x}) - c_{j i}\phi^i(\vec{x})\right), \quad & & \vec{x} \in \Gamma_{i j},
    \quad (i, j) \in \{1, \dots, N_\text{cmpt}\}^2, &\\
    \mathbf{D}_i \nabla \phi^i(\vec{x}) \cdot \vec{n}_i(\vec{x}) & = -\kappa_i \phi^i(\vec{x}),
    \quad & & \vec{x} \in \Gamma_i, \quad i \in \{1, \dots, N_\text{cmpt}\}. &
\end{alignedat}
```

Let the solutions ``(\phi, \lambda)`` to the above equations be denoted by ``\left\{(\phi_n,
\lambda_n)\right\}_{n \in \mathbb{N}^*}``. We assume the non-negative real-valued
eigenvalues are ordered in non-decreasing order:

```math
0 = \lambda_1 \leq \lambda_2 \leq \lambda_3 \leq \dots
```

If the domain ``\Omega`` consists of only one contiguous group of compartments connected
through a chain of permeable membranes, only the first eigenvalue will be zero, and the
corresponding eigenfunction will be the only constant function. If there are
``N_\text{group} \geq 2`` groups of connected compartments completely separated by interior
hard wall membranes, the first ``N_\text{group}`` eigenvalues will be zero:

```math
0 = \lambda_1 = \dots = \lambda_{N_\text{group}} < \lambda_{N_\text{group} + 1} \leq \dots,
```

and there will be ``N_\text{group}`` corresponding groupwise constant eigenfunctions. In the
latter case, the equations may be rewritten separately for each connected subdomain to
obtain a set of eigenvalues with a multiplicity of one, but the formulation is also valid in
the global form with multiple zero eigenvalues. These two formulations will lead to an
identical eigenfunction basis (up to a linear combination for the eigenvalues with
multiplicity higher than one), where the basis of each subdomain is a subset of the global
eigenfunction basis.

Let ``\mathbf{L}`` be the diagonal matrix containing the first ``N_\text{eig}`` Laplace
eigenvalues:

```math
\mathbf{L} = \operatorname{diag}(\lambda_1, \lambda_2, \dots, \lambda_{N_\text{eig}})\in
\R^{N_\text{eig} \times N_\text{eig}}.
```

Then the matrix ``\mathbf{L}`` represents the generalized Laplace operator ``-\nabla \cdot \mathbf{D}
\nabla`` in the truncated Laplace eigenfunction basis.


## Bloch-Torrey operator in the Laplace eigenfunction basis

Let ``\mathbf{A}(\vec{g})`` be the ``N_\text{eig} \times N_\text{eig}`` matrix defined by:

```math
\mathbf{A}(\vec{g}) = g_x \mathbf{A}^x + g_y \mathbf{A}^y + g_z \mathbf{A}^z,
```

where ``\vec{g} = (g_x, g_y, g_z)^\mathsf{T}`` is the gradient vector and ``\mathbf{A}^x``,
``\mathbf{A}^y``, and ``\mathbf{A}^z`` are three symmetric ``N_\text{eig} \times N_\text{eig}``
matrices whose entries are the first order moments in the coordinate directions of the
product of pairs of eigenfunctions:

```math
\begin{alignedat}{3}
    A_{mn}^x & = \int_\Omega x \, \phi_m(\vec{x}) \phi_n(\vec{x}) \, \mathrm{d}
    \Omega(\vec{x}), \quad & & (m, n) \in \{1, \dots, N_\text{eig}\}^2, \\
    A_{mn}^y & = \int_\Omega y \, \phi_m(\vec{x}) \phi_n(\vec{x}) \, \mathrm{d}
    \Omega(\vec{x}), \quad & & (m, n) \in \{1, \dots, N_\text{eig}\}^2, \\
    A_{mn}^z & = \int_\Omega z \, \phi_m(\vec{x}) \phi_n(\vec{x}) \, \mathrm{d}
    \Omega(\vec{x}), \quad & & (m, n) \in \{1, \dots, N_\text{eig}\}^2.
\end{alignedat}
```

Similarly, let ``\mathbf{T}`` be the ``N_\text{eig} \times N_\text{eig}`` Laplace relaxation matrix
defined by

```math
T_{mn} = \int_\Omega \frac{1}{T_2(\vec{x})} \phi_m(\vec{x}) \phi_n(\vec{x}) \, \mathrm{d}
\Omega(\vec{x}), \quad (m,n) \in \{1, \dots, N_\text{eig}\}^2.
```

Then the general time dependent Bloch-Torrey operator

```math
-\nabla \cdot \mathbf{D} \nabla + \frac{1}{T_2} + \underline{\mathrm{i}} \gamma \vec{g}(t)
\cdot \vec{x}, \quad t \in [0, T_\text{echo}]
```

in the truncated Laplace eigenfunction basis ``(\phi_j)_{j = 1, \dots, N_\text{eig}}`` is
given by the complex-valued matrix

```math
\mathbf{K}(\vec{g}(t)) = \mathbf{L} + \mathbf{T} + \underline{\mathrm{i}} \gamma
\mathbf{A}(\vec{g}(t)), \quad t \in [0, T_\text{echo}].
```

## Matrix Formalism approximation

Denoting ``\vec{\phi} = (\phi_1, \dots, \phi_{N_\text{eig}})^\mathsf{T}`` the vector of the
first ``N_\text{eig}`` Laplace eigenfunctions, we may project the solution to to the
Bloch-Torrey PDE onto the truncated Laplace eigenfunction basis:

```math
M^\text{MF}(\vec{x}, t) = \vec{\phi}^\mathsf{T}\!(\vec{x}) \vec{\nu}(t),
```

where the vector of coefficients ``\vec{\nu} = (\nu_1, \dots, \nu_{N_\text{eig}})^\mathsf{T}
: [0, T_\text{echo}] \to \mathbb{C}^{N_\text{eig}}`` is the solution to

```math
\frac{\mathrm{d} \vec{\nu}}{\mathrm{d} t} = -\mathbf{K}(\vec{g}(t)) \vec{\nu}(t), \quad t
\in [0, T_\text{echo}].
```

The initial coefficients are given by

```math
\vec{\nu}(0) = \int_\Omega \rho(\vec{x}) \vec{\phi}(\vec{x}) \, \mathrm{d} \Omega(\vec{x}).
```

By using a piece-wise constant approximation of the gradient ``\vec{g}``, we obtain

```math
\vec{\nu}(T_\text{echo}) \approx \left(\prod_{i = 1}^{N_\text{int}} \mathrm{e}^{-\delta_i
\mathbf{K}_i}\right) \vec{\nu}(0) = \mathrm{e}^{-\delta_{N_\text{int}}
\mathbf{K}_{N_\text{int}}} \dots \mathrm{e}^{-\delta_2 \mathbf{K}_2} \mathrm{e}^{-\delta_1
\mathbf{K}_1} \vec{\nu}(0),
```

where ``\{I_i\}_{i = 1, \dots, N_\text{int}}`` are intervals such that ``[0, T_\text{echo}]
= \bigcup_{i = 1}^{N_\text{int}} I_i``, ``\vec{g}(t) = \vec{g}_i`` for ``t \in I_i``,
``\delta_i = |I_i|``, and ``\mathbf{K}_i = \mathbf{K}(\vec{g}_i)``. The constants may be
computed through quadrature:

```math
\vec{g}_i = \frac{1}{\delta_i}\int_{I_i} \vec{g}(t) \, \mathrm{d} t \approx \frac{1}{2}
\left(\vec{g}(\min I_i) + \vec{g}(\max I_i)\right).
```

For the PGSE and double-PGSE sequences ``\vec{g}(t) = f(t) g \vec{d}``, where ``f(t) \in
\{-1, 0, 1\}`` for all ``t``, the coefficients of the final magnetization are given by

```math
\vec{\nu}(T_\text{echo}) = \mathrm{e}^{-\delta\mathbf{K}(g \vec{d})^*} \mathrm{e}^{-(\Delta -
\delta)(\mathbf{L} + \mathbf{T})} \mathrm{e}^{-\delta \mathbf{K}(g \vec{d}} \vec{\nu}(0)
```

and

```math
\vec{\nu}(T_\text{echo}) = \mathrm{e}^{-\delta \mathbf{K}(g \vec{d})^*} \mathrm{e}^{-(\Delta
- \delta)(\mathbf{L} + \mathbf{T})} \mathrm{e}^{-\delta \mathbf{K}(g \vec{d})}
\mathrm{e}^{-t_\text{pause}(\mathbf{L}
+ \mathbf{T})} \mathrm{e}^{-\delta \mathbf{K}(g \vec{d})^*} \mathrm{e}^{-(\Delta -
+ \delta)(\mathbf{L} +
\mathbf{T})} \mathrm{e}^{-\delta \mathbf{K}(g \vec{d})} \vec{\nu}(0)
```

respectively. These are the exact solutions, although there may still be truncation errors
from ``N_\text{eig}``.

The signal is given by

```math
S^\text{MF}(f, \vec{g}) = \vec{\nu}^\mathsf{T}\!(T_\text{echo}) \int_\Omega
\vec{\phi}(\vec{x}) \, \mathrm{d} \Omega(\vec{x}).
```
## Apparent diffusion coefficient

From the Matrix Formalism signal, the analytical expression of its ADC for a gradient
direction given by a unit vector ``\vec{d}`` and a sequence ``f`` is the following:

```math
D^\text{MF}(\vec{d}, f) = \frac{1}{|\Omega|} \sum_{n = 1}^{N_\text{eig}} (\vec{d} \cdot
\vec{a}_n)^2 \lambda_n \frac{\int_0^{T_\text{echo}} F(t) \int_0^t \mathrm{e}^{-\lambda_n(t -
s)} f(s) \, \mathrm{d} s \, \mathrm{d} t}{ \int_0^{T_\text{echo}} F^2(t) \, \mathrm{d} t},
```

where the coefficients of the three-row matrix ``\mathbf{a} = (\vec{a}_1^, \dots,
\vec{a}_{N_\text{eig}}) \in \mathbb{R}^{3 \times N_\text{eig}}`` are given by

```math
\begin{split}
    a_n^x & = \int_{\Omega} x \phi_n(\vec{x}) \, \mathrm{d} \Omega(\vec{x}), \\
    a_n^y & = \int_{\Omega} y \phi_n(\vec{x}) \, \mathrm{d} \Omega(\vec{x}), \\
    a_n^z & = \int_{\Omega} z \phi_n(\vec{x}) \, \mathrm{d} \Omega(\vec{x}),
\end{split}
\quad n \in \{1, \dots, N_\text{eig}\}.
```

To clarify the relationship between the ADC and the diffusion encoding direction
``\vec{d}``, we rewrite the Matrix Formalism ADC as

```math
D^\text{MF}(\vec{d}, f) = \vec{d}^\mathsf{T} \mathbf{D}^{\text{MF}}(f) \vec{d},
```

where the Matrix Formalism effective diffusion tensor is given by

```math
\mathbf{D}^{\text{MF}}(f) = \frac{1}{|\Omega|} \sum_{n = 1}^{N_\text{eig}} j_n(f) \vec{a}_n
\vec{a}_n^\mathsf{T} = \frac{1}{|\Omega|} \mathbf{a} \mathbf{j}(f) \mathbf{a}^\mathsf{T},
```

with ``\mathbf{j}(f) = \operatorname{diag}\left(j_1(f), \dots, j_{N_\text{eig}}(f)\right) \in
\mathbb{R}^{N_\text{eig} \times N_\text{eig}}`` depending on ``f``:

```math
j_n(f) = \lambda_n \frac{\int_0^{T_\text{echo}} F(t) \int_0^t e^{-\lambda_n(t - s)} f(s) \,
\mathrm{d} s \, \mathrm{d} t}{\int_0^{T_\text{echo}} F^2(t) \, \mathrm{d} t}, \quad n \in
\{1, \dots, N_\text{eig}\}.
```

In the case of a constant initial spin density ``\rho``, we also allow for computing the
Matrix Formalism Gaussian Approximation (MFGA) signal, given as

```math
S^{\text{MFGA}}(g, f, \vec{d}) = \rho |\Omega| \exp\left(-\vec{d}^\mathsf{T}
\mathbf{D}^{\text{MF}}(f) \vec{d} \; b(g, f) \right) }.
```

## Eigenfunction length scale and orientation

On a line segment of length ``L`` and diffusivity ``D_0``, the eigenvalues ``(\lambda_1,
\lambda_2, \dots)`` of the Laplace operator with Neumann boundary conditions are

```math
\lambda_n = \left(\frac{\pi (n - 1)}{L}\right)^2 D_0, \quad n = 1, 2, \dots
```

To make the link between the computed eigenvalue and the spatial scale of the eigenmode, we
will convert the computed ``\lambda_n`` into a length scale:

```math
\ell(\lambda) =
\begin{cases}
    + \infty, \quad                                & \lambda = 0, \\
    \pi \sqrt{\frac{\bar{\sigma}}{\lambda}}, \quad & \lambda > 0,
\end{cases}
```

and characterize the computed eigenmode by ``\ell(\lambda_n)`` instead of ``\lambda_n``. The
reference diffusivity ``\bar{\sigma}`` is taken as a volume weighted mean of the trace
average (spherical part) of ``(\mathbf{D}_i)_{1 \leq i \leq N_\text{cmpt}}``:

```math
\bar{\sigma} = \frac{1}{|\Omega|} \int_\Omega \frac{1}{3}\operatorname{tr} \mathbf{D}(\vec{x})
\, \mathrm{d} \Omega(\vec{x}).
```

To characterize the directional contribution of the eigenmode we use the fact that its
contribution to the ADC in the direction ``\vec{d}`` is ``j_n(f) (\vec{d} \cdot
\vec{a}_n)^2``. We thus call ``\vec{a}_n = (a_n^x, a_n^y, a_n^z)^\mathsf{T}`` the
"diffusion direction" of the ``n``th eigenmode. We remind that the three components of
``\vec{a}_n`` are the first moments in the 3 canonical basis directions of the associated
eigenfunction.
