"""
Time stepping on NVIDIA GPUs: dense arrays on the GPU, and sparse matrices in the CSR format
with the sparse direct solvers of cuDSS.
"""
module SpinDoctorCUDAExt

using CUDA: CuArray
using CUDA.CUSPARSE: CuSparseMatrixCSR
using CUDSS
using LinearAlgebra: cholesky, lu
using SparseArrays: SparseMatrixCSC
using SpinDoctor

SpinDoctor.to_device(::Type{<:CuArray}, x::Array) = CuArray(x)
SpinDoctor.to_device(::Type{<:CuArray}, A::SparseMatrixCSC) = CuSparseMatrixCSR(A)
SpinDoctor.device_lu(::Type{<:CuArray}, A::SparseMatrixCSC) = lu(CuSparseMatrixCSR(A))
SpinDoctor.device_cholesky(::Type{<:CuArray}, A::SparseMatrixCSC) =
    cholesky(CuSparseMatrixCSR(A))

end
