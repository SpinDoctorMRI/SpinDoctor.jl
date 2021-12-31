# Definitions

## Abbreviations

- ADC, Apparent Diffusion Coefficient;
- BT, Bloch-Torrey;
- BTPDE, Bloch-Torrey PDE;
- DMRI, Diffusion MRI;
- ECS, Extra-Cellular Space;
- FEM, Finite Element Method;
- HADC, Homogenized ADC model;
- HARDI, High Angular Resolution Diffusion Imaging;
- MF, Matrix Formalism;
- MRI, Magnetic Resonance Imaging;
- ODE, Ordinary Differential Equation;
- OGSE, Oscillating Gradient Spin Echo sequence;
- PDE, Partial Differential Equation;
- PGSE, Pulsed-Gradient Spin Echo sequence;
- STA, Short Time Approximation.


## Geometrical configuration

In SpinDoctor, we consider a domain $\Omega = \bigcup_{i = 1}^{N_\text{cmpt}}
\Omega_i \subset \mathbb{R}^3$ consisting of $N_\text{cmpt}$ compartments $\{\Omega_i\}_{1
\leq i \leq N_\text{cmpt}}$. The permeable interface between two compartments is denoted by
$\Gamma_{i j} = \Omega_i \cap \Omega_j$ for $i \neq j$, $(i, j) \in \{1, \dots,
N_\text{cmpt}\}^2$. For $i = j$, we let $\Gamma_{ii} = \emptyset$ for the ease of notation.
Finally, let $\partial \Omega$ denote the outer boundary of the domain, and $\Gamma_i =
\Omega_i \cap \partial \Omega$ its restriction to $\Omega_i$. Note that for compartments
that do not touch, we have $\Gamma_{i j} = \emptyset$. Similarly, we have $\Gamma_i =
\emptyset$ for compartments that do not touch the outer boundary.
