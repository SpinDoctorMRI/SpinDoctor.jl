# Custom gradients

Here we will consider a custom gradient applie to a sphere geometry.

We start by loading SpinDoctor.

```julia
using LinearAlgebra
using SpinDoctor
```

The built in geometry recipes allow for making various cell configuration. We here consider
the case of twisted axons immersed in an extracellular space (ECS).

```julia
setup = SphereSetup(;
    name = "gaze-into-the-orb",
    ncell = 1,
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.3,
    bend = 0.0,
    twist = 9 / 4,
    include_in = false,
    in_ratio = 0.6,
    ecs_shape = :no_ecs,
    ecs_ratio = 0.5,
)
```

```julia
coeffs = coefficients(
    setup;
    D = (; in = 0.002 * I(3), out = 0.002 * I(3), ecs = 0.002 * I(3)),
    T₂ = (; in = Inf, out = Inf, ecs = Inf),
    ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
    κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    γ = 2.67513e-4,
)
```

```julia
mesh, = create_geometry(setup)
```

```julia
model = Model(; mesh, coeffs...)
matrices = assemble_matrices(model)
```

The Bloch-Torrey PDE takes a magnetic field gradient pulse sequence as an input. We may
define our custom gradient (given in T/m) and echo time `TE` (given in microseconds) as
follows:

```julia
TE = 5000.0
g⃗(t) = 0.5 * [sin(10π * t / TE), 0, 0]
gradient = GeneralGradient(g⃗, TE)
```

In order to follow the evolution of the solution during time stepping, we add a `Plotter` to
a list of callbacks.

```julia
callbacks = [Plotter{Float64}()]
```

We may then define the problem and solve for our gradient and callbacks.

```julia
btpde = GeneralBTPDE(; model, matrices)
ξ = solve(btpde, gradient; callbacks)
```
