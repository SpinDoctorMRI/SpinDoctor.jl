# Finite element discretization

In SpinDoctor, the finite element mesh generation is performed using an external package
called Tetgen. Each finite element mesh consists of

- a list of ``N_\text{node}`` nodes in three dimensions: ``(\mathbf{q}_1, \dots,
  \mathbf{q}_{N_\text{node}}) = (\mathbf{q}^x, \mathbf{q}^y, \mathbf{q}^z)^\mathsf{T} \in
  \R^{3 \times N_\text{node}}``;
- a list of ``N_\text{element}`` tetrahedral elements (``4 \times N_\text{element}`` indices
  referencing the nodes).

The list of nodes includes double nodes that are placed at the interfaces between
compartments connected by permeable membranes. This allows for representing discontinuous
magnetization fields ``M_i`` and ``M_j`` (or ``\phi^i`` and ``\phi^j``) on the same boundary
``\Gamma_{i j}``. To distinguish between the different compartments, let ``\{1, \dots,
N_\text{node}\} = \bigcup_{i = 1}^{N_\text{cmpt}} \mathcal{I}_i`` with ``\mathcal{I}_i \cap
\mathcal{I}_j = \emptyset`` for ``i \neq j``. The set ``\mathcal{I}_i`` contains the indices of
the nodes representing compartment ``\Omega_i``, including interface nodes.  In the adjacent
compartments, the corresponding interface nodes will have different indices, distinct from
``\mathcal{I}_i``.

In SpinDoctor, the finite element space is the space of compartment-wise continuous
piecewise linear functions on tetrahedral elements in three dimensions. This space has a set
of basis functions whose number is exactly the number of finite element nodes (including
double nodes), and that are defined on the entire domain ``\Omega``:

```math
\varphi_k : \Omega \to [0, 1], \quad k \in \{1, \dots, N_\text{node}\}.
```

Let the finite element nodes be denoted by ``\mathbf{q}_1, \dots, \mathbf{q}_{N_\text{node}}``. The
basis function ``\varphi_k``, ``k \in \mathcal{I}_i``, is a piece-wise linear function, non-zero
on the tetrahedra of ``\Omega_i`` that touch the node ``\mathbf{q}_k``, and zero on all other
tetrahedra (including tetrahedra of other compartments different than ``\Omega_i`` that do
touch ``\mathbf{q}_k``). At the interface ``\Gamma_{i j}`` between two compartments, the value of
``\varphi_k`` is set to be the value it has inside its own compartment, distinct from that of
the adjacent compartment. On a tetrahedron of ``\Omega_i`` that touches ``\mathbf{q}_k``,
``\varphi_k`` is equal to ``1`` on ``\mathbf{q}_k`` and it is equal to ``0`` on the other 3 vertices of
the tetrahedron. This completely describes the piece-wise linear function. The index sets
may then be defined by ``\mathcal{I}_i = \{k = 1, \dots, N_\text{node} \ | \ \operatorname{supp}
(\varphi_k) \subset \Omega_i\}``, the set of indices of the finite element nodal functions
whose supports lie entirely within ``\Omega_i``.

Any function ``u`` in the finite element space can be written as a linear combination of the
above basis functions:

```math
u(\vec{x}) = \sum_{k = 1}^{N_\text{node}} \alpha_k \varphi_k(\vec{x}) =
\bm{\alpha}^\mathsf{T} \bm{\varphi}(\vec{x}), 
```

where ``\bm{\alpha} = (\alpha_1, \dots, \alpha_{N_\text{node}})^\mathsf{T}`` is the vector of
coefficients and ``\bm{\varphi} = (\varphi_1, \dots, \varphi_{N_\text{node}})^\mathsf{T}`` is
the vector of finite element nodal basis functions.

To discretize the Bloch-Torrey and Laplace operators, we construct the following finite
element matrices: ``\mathbf{M},\mathbf{S},\mathbf{Q}\in\R^{N_\text{node}\times N_\text{node}}``, known
in the FEM literature as the mass, stiffness, and flux matrices, respectively. These
matrices matrices we need are defined as follows, for ``(k, l) \in \{1, \dots,
N_\text{node}\}``:

```math
M_{kl} = \int_\Omega \varphi_k(\vec{x}) \varphi_l(\vec{x}) \, \mathrm{d}
\Omega(\vec{x}),
```

```math
S_{kl} = \int_\Omega \mathbf{D}(\vec{x}) \nabla \varphi_k(\vec{x}) \cdot \nabla
\varphi_l(\vec{x}) \, \mathrm{d} \Omega(\vec{x}),
```

```math
Q_{kl} = \sum_{i = 1}^{N_\text{cmpt}} Q_{kl}^i + \sum_{i = 1}^{N_\text{cmpt}} \sum_{j =
1}^{N_\text{cmpt}} Q_{kl}^{i j},
```

the latter being defined as the sum of interface integrals:

```math
Q_{kl}^i = \kappa_i \int_{\Gamma_i} \varphi_k(\vec{x}) \varphi_l(\vec{x}) \, \mathrm{d}
\Gamma(\vec{x}), \quad i \in \{1, \dots, N_\text{cmpt}\},
```

```math
Q_{kl}^{i j} =
\begin{cases}
    \kappa_{i j} c_{j i} \int_{\Gamma_{i j}} \varphi_k(\vec{x}) \varphi_l(\vec{x}) \,
    \mathrm{d} \Gamma(\vec{x}), \quad  & (k, l) \in \mathcal{I}_i^2, \\
    -\kappa_{i j} c_{i j} \int_{\Gamma_{i j}} \varphi_k(\vec{x}) \varphi_l(\vec{x}) \,
    \mathrm{d} \Gamma(\vec{x}), \quad & (k, l) \in \mathcal{I}_i \times \mathcal{I}_j, \\
    0, \quad & \text{otherwise}.
\end{cases}
\quad (i, j) \in \{1, \dots, N_\text{cmpt}\}^2.
```

We remind the reader that ``i`` and ``j`` are compartment indices, while ``k`` and ``l`` are finite
element nodal indices. Note that the above sum formulation counts each interface twice, but
only considers the contribution to one side of the boundary at the time. It also correctly
accounts for corner nodes (if any) that belong to two different permeable boundaries at the
same time.

The weak form of the Bloch-Torrey PDE also requires computing the first order moments of the
product of pairs of finite element basis functions. We let these three matrices be denoted
by ``\mathbf{J}^x``, ``\mathbf{J}^y`` and ``\mathbf{J}^z``:

```math
J_{k l}^x = \int_\Omega x \, \varphi_k(\vec{x}) \varphi_l(\vec{x}) \, \mathrm{d}
\Omega(\vec{x}), \quad (k, l) \in \{1, \dots, N_\text{node}\}^2,
```

```math
J_{k l}^y = \int_\Omega y \, \varphi_k(\vec{x}) \varphi_l(\vec{x}) \, \mathrm{d}
\Omega(\vec{x}), \quad (k, l) \in \{1, \dots, N_\text{node}\}^2,
```

```math
J_{k l}^z = \int_\Omega z \, \varphi_k(\vec{x}) \varphi_l(\vec{x}) \, \mathrm{d}
\Omega(\vec{x}), \quad (k, l) \in \{1, \dots, N_\text{node}\}^2.
```

where ``\vec{x} = (x, y, z)^\mathsf{T}``. For a given gradient vector ``\vec{g} = (g_x, g_y,
g_z)^\mathsf{T}``, we define

```math
\mathbf{J}(\vec{g}) = g_x \mathbf{J}^x + g_y \mathbf{J}^y + g_z \mathbf{J}^z.
```

In addition, we define the finite element relaxation matrix ``\mathbf{R} \in
\mathbb{R}^{N_\text{node} \times N_\text{node}}`` given by

```math
R_{kl} = \int_\Omega \frac{1}{T_2(\vec{x})} \, \varphi_k(\vec{x}) \varphi_l(\vec{x}) \,
\mathrm{d} \Omega(\vec{x}), \quad (k, l) \in \{1, \dots, N_\text{node}\}^2.
```

In SpinDoctor, these matrices are assembled from local element matrices and the assembly
process is based on vectorized routines of, which replace expensive
loops over elements by operations with 3-dimensional arrays. All local element matrices in
the assembly of ``\mathbf{S}``, ``\mathbf{M}``, and ``\mathbf{Q}`` are evaluated at the same time and
stored in a full matrix of size ``4 \times 4 \times N_\text{element}``, where
``N_\text{element}`` denotes the number of tetrahedral elements.

The matrices ``\mathbf{J}^x``, ``\mathbf{J}^y``, and ``\mathbf{J}^z`` are assembled as coordinate weighted
mass matrices, where the three coordinate functions ``\vec{x} \mapsto x``, ``y``, ``z`` act as
nodal weights in the assembly process, given by ``\mathbf{q}^x``, ``\mathbf{q}^y``, ``\mathbf{q}^z``; the
vectors of ``x``, ``y``, and ``z`` coordinates of the finite element nodes.

With compartment-wise constant ``T_2``-relaxation times, \mathbf{R} is a block-diagonally scaled
version of the mass matrix \mathbf{M}, where the weights are the inverses of the relaxation
times.



## Finite element solution to the Bloch-Torrey PDE

The solution to the Bloch-Torrey partial differential equation may be
projected onto the finite element nodal basis ``(\varphi_k)_{1\leq k\leq N_\text{node}}``, in
which case it is given by

```math
M(\vec{x},t) = \sum_{k = 1}^{N_\text{node}}
\xi_k(t)\varphi_k(\vec{x}) = \bm{\xi}^\mathsf{T}\!(t) \bm{\varphi}(\vec{x}),
```

where ``\bm{\varphi} = (\varphi_1,\dots,\varphi_{N_\text{node}})^\mathsf{T} \in
\mathbb{R}^{N_\text{node}}`` is the vector of finite element basis functions and the
function ``\bm{\xi} = (\xi_1, \dots, \xi_{N_\text{node}})^\mathsf{T} : [0, T_\text{echo}]
\to \mathbb{C}^{N_\text{node}}`` is the solution to the following system of _ordinary_
differential equations (ODE):

```math
\mathbf{M} \frac{\mathrm{d} \bm{\xi}}{\mathrm{d} t} = -\left(\mathbf{S} + \mathbf{Q} + \mathbf{R} +
\underline{\mathrm{i}} \gamma f(t) \mathbf{J}(\vec{g})\right)\bm{\xi}(t).
```

## Finite element solution to the HADC model

Similarly, the solution to the HADC model, ``\omega``, may be obtained by
solving the equation

```math
\mathbf{M} \frac{\mathrm{d} \bm{\zeta}}{\mathrm{d} t} = - \mathbf{S} \bm{\zeta}(t) + F(t)
\mathbf{G} \mathbf{D} \vec{d},
```

where ``\bm{\zeta} = (\zeta_1, \dots, \zeta_{N_\text{node}}) : [0, T_\text{echo}] \to
\mathbb{R}^{N_\text{node}}`` is the unknown vector of coefficients of the solution
``\omega(\vec{x}, t) = \bm{\zeta}^\mathsf{T}(t) \bm{\varphi}(\vec{x})``, ``\bm{\varphi} =
(\varphi_1, \dots, \varphi_{N_\text{node}})`` is the vector of finite element basis
functions, and ``\mathbf{G} \in \mathbb{R}^{N_\text{node} \times 3}`` is given by

```math
\mathbf{G} = \sum_{k = 1}^{N_\text{node}} \int_{\partial \Omega} \varphi_k(\vec{x})
\bm{\varphi}(\vec{x}) \vec{n}^\mathsf{T}\!(\vec{x}) \, \mathrm{d}\Gamma(\vec{x}),
```

where ``\vec{n} = (n_x, n_y, n_z)^\mathsf{T}`` is the outwards unit normal. This three-column
matrix represents the components of the boundary integral of the quantity ``F(t)
\vec{d}^\mathsf{T} \mathbf{D} \vec{n}``, where the constant diffusivity and gradient
sequence dependent part ``F(t) \vec{d}^\mathsf{T} \mathbf{D}`` has been factored out.
They can thus be assembled independently of the gradient sequence, and be reused when
solving for multiple sequences or directions. They are assembled using the same routine as
for \mathbf{Q}.

The boundary integral ``h`` is then given by

```math
h(t) = \frac{1}{|\Omega|} \bm{\zeta}^\mathsf{T}\!(t) \mathbf{G} \vec{d},
```

where ``|\Omega| = \sum_{j,k = 1}^{N_\text{node}} M_{jk}`` is computed using the mass matrix.


## Finite element Matrix Formalism signal representation

Both of the above ODEs are of dimension ``N_\text{node}``, the number of finite element
nodes. The finite element discretized Matrix Formalism representation uses a different
approach, limiting the problem size to ``N_\text{eig}``, the number of Laplace
eigenfunctions (of which the choice is further explored in section. But first we need to
solve an eigenvalue problem involving matrices of size ``N_\text{node}\times
N_\text{node}``. The finite element discretization changes the continuous Laplace operator
eigenvalue problem into the following discrete, generalized _matrix_ eigenvalue problem:
find ``(\lambda, \mathbf{p}) \in \mathbb{R} \times \mathbb{R}^{N_\text{node}}`` such that

```math
\lambda \mathbf{M} \mathbf{p} = (\mathbf{S} + \mathbf{Q}) \mathbf{p},
```

of which we will retain the ``N_\text{eig}`` smallest eigenvalues and corresponding
eigenvectors ``\{(\lambda_n, \mathbf{p}_n)\}_{1 \leq n \leq N_\text{eig}}``, with ``N_\text{eig}
\leq N_\text{node}``. Note however that there are in total ``N_\text{node}`` solutions to the
problem. Moving back to the space of functions (the function space
``P_1``), the eigenfunction ``\phi_n(\vec{x})`` associated to the eigenvalue ``\lambda_n`` is then

```math
\phi_n(\vec{x}) = \sum_{k = 1}^{N_\text{node}} p_n^k \varphi_k(\vec{x}) =
\mathbf{p}_n^\mathsf{T} \bm{\varphi}(\vec{x}), \quad n \in \{1, \dots, N_\text{eig}\},
```

where the entries of the eigenvector ``\mathbf{p}_n`` are the coefficients of the eigenfunction
``\phi_n`` in the finite element basis. Using matrix notation, this conversion can also be
written ``\bm{\phi} = \mathbf{P}^\mathsf{T} \bm{\varphi}``, where ``\bm{\varphi} = (\varphi_1,
\dots, \varphi_{N_\text{node}})^\mathsf{T}``, ``\bm{\phi} = (\phi_1, \dots,
\phi_{N_\text{eig}})^\mathsf{T}``, and ``\mathbf{P} = (\mathbf{p}_1, \dots, \mathbf{p}_{N_\text{eig}})
\in \mathbb{R}^{N_\text{node} \times N_\text{eig}}``. The integrals of the finite element
discretized eigenfunctions are then given by ``\bm{\Phi} = \int_\Omega \bm{\phi}(\vec{x})
\, \mathrm{d} \Omega(\vec{x}) = \mathbf{P}^\mathsf{T} \mathbf{M} \mathbf{o}``, where ``\mathbf{o} = (1,
\dots, 1)^\mathsf{T} \in \mathbb{R}^{N_\text{node}}``. Similarly, the coefficients of the
initial spin density in the finite element discretized eigenfunction basis are given by
``\bm{\nu} = \int_\Omega \rho(\vec{x}) \bm{\phi}(\vec{x}) \, \mathrm{d} \Omega(\vec{x}) =
\mathbf{P}^\mathsf{T} \mathbf{M} \bm{\rho}``, where ``\bm{\rho} = (\rho_{i(k)})_{1 \leq k \leq
N_\text{node}} \in \mathbb{R}^{N_\text{node}}`` and ``i(k) \in \{1, \dots, N_\text{cmpt}\}`` is
such that ``k \in \mathcal{I}_{i(k)}``.

Then it is clear that the first order moments of the product of pairs of Laplace
eigenfunctions can be written as:

```math
\mathbf{A}^u = \mathbf{P}^\mathsf{T}\mathbf{J}^u\mathbf{P}, \quad u = x, y, z.
```

These three matrices are computed using a total of six matrix-matrix multiplications.
Similarly, the eigenfunction basis relaxation matrix \mathbf{T} may be obtained from the finite
element nodal basis relaxation matrix \mathbf{R} by

```math
\mathbf{T} = \mathbf{P}^\mathsf{T} \mathbf{R} \mathbf{P}.
```
