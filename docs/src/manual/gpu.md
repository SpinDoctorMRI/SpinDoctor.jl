# GPU

The time stepping of the Bloch-Torrey equation and of the homogenized ADC model runs on
NVIDIA GPUs, with the sparse direct solvers of
[cuDSS](https://developer.nvidia.com/cudss), after loading CUDA and CUDSS:

```julia
using CUDA, CUDSS

solver = ThetaSolver(; timestep = 50.0, arraytype = CuArray)
u = solve_btpde(model, matrices, gradient, solver)
adcs = solve_hadc(model, matrices, gradient; arraytype = CuArray)
```

The matrices are assembled on the CPU and moved to the GPU, and the result is returned on
the CPU. Periodic boundary conditions are supported.

For packed myelinated axons with 73 000 mesh points, an RTX 4090 compared with 8 CPU
threads:

| Problem                          | CPU   | GPU    |
| :------------------------------- | ----: | -----: |
| BTPDE, PGSE, 500 time steps      | 61 s  | 3.3 s  |
| BTPDE, cos-OGSE, 220 time steps  | 161 s | 11 s   |
| HADC tensors                     | 8.1 s | 0.74 s |

The GPU tests are in `test/gpu`, with their own environment:

```bash
julia --project=test/gpu test/gpu/runtests.jl
```
