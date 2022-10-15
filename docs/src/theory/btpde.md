# Bloch-Torrey partial differential equation

In diffusion MRI, a time-varying magnetic field gradient ``\vec{g} : [0, T_\text{echo}] \to
\mathbb{R}^3`` is applied to the tissue to encode water diffusion. The resulting complex transverse
water proton magnetization in the rotating frame satisfies the Bloch-Torrey PDE:

```math
\frac{\partial}{\partial t} M_i(\vec{x},t) = -\left(-\nabla \cdot \mathbf{D}_i \nabla +
\frac{1}{T_{2,i}} + \underline{\mathrm{i}} \gamma \vec{g}(t) \cdot \vec{x}\right)
M_i(\vec{x}, t),
```

where

- ``\vec{x} \in \Omega_i``,
- ``t \in [0, T_\text{echo}]`` where ``T_\text{echo}`` is the echo time at which the signal is
  measured,
- ``\gamma = 2.67513\times 10^8 \, \mathrm{rad \, s^{-1} T^{-1}}`` is the gyromagnetic ratio
  of the water proton,
- ``\underline{\mathrm{i}}`` is the imaginary unit,
- ``\mathbf{D}_i`` is the intrinsic diffusion tensor in ``\Omega_i``,
- ``T_{2,i}`` is the ``T_2``-relaxation time in ``\Omega_i``, and
- ``M_i`` is the magnetization in ``\Omega_i``.

The initial conditions are given by

```math
M_i(\vec{x}, 0) = \rho_i \in \mathbb{C}, \quad i \in \{1, \dots, N_\text{cmpt}\}.
```

where ``\rho_i`` is the initial spin density in ``\Omega_i``. The outer boundary conditions are
given by

```math
\mathbf{D}_i \nabla M_i(\vec{x}, t) \cdot \vec{n}_i(\vec{x}) = -\kappa_i M_i(\vec{x}, t),
\quad \vec{x} \in \Gamma_i, \quad i \in \{1, \dots, N_\text{cmpt}\},
```

and the interface conditions by

```math
\mathbf{D}_i \nabla M_i(\vec{x}, t) \cdot \vec{n}_i(\vec{x}) = -\mathbf{D}_j \nabla
M_j(\vec{x}, t) \cdot \vec{n}_j(\vec{x}), \quad \vec{x} \in \Gamma_{i j}, \quad (i, j) \in
\{1, \dots, N_\text{cmpt}\}^2,
```

```math
\mathbf{D}_i \nabla M_i(\vec{x}, t) \cdot \vec{n}_i(\vec{x}) = \kappa_{i j} \left(c_{i j}
M_j(\vec{x}, t) - c_{j i} M_i(\vec{x}, t)\right), \quad \vec{x} \in \Gamma_{i j}, \quad (i,
j) \in \{1, \dots, N_\text{cmpt}\}^2.
```

where

- ``\vec{n}_i`` is the unit outward pointing normal vector of ``\Omega_i``
- ``\kappa_i`` is a boundary relaxation coefficient for ``\Omega_i``,
- ``\kappa_{i j} = \kappa_{j i}`` is the permeability coefficient on ``\Gamma_{i j}``,
- ``c_{i j}`` and ``c_{j i}`` account for the spin density equilibrium between the two
  compartments, with either
  - ``c_{i j} = c_{j i} = 1``, in which case a uniform spin density across compartments is
    favored in the absence of a
    gradient, or
  - ``c_{i j} = \frac{2 \rho_i}{\rho_i + \rho_j}`` and ``c_{j i} = \frac{2 \rho_j}{\rho_i +
    \rho_j}``, which ensures that the non-uniform intitial spin density is preserved in the
    absence of a gradient [Lee2021](@cite).

The diffusion MRI signal is measured at echo time ``t = T_\text{echo}``. This signal is
computed as the spatial integral of the final magnetization ``M(\cdot,
T_\text{echo})``:

```math
S(\vec{g}) = \int_\Omega M(\vec{x}, T_\text{echo}) \, \mathrm{d} \Omega(\vec{x}).
```
