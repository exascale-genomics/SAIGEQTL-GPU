/*! @file gpuSparseMult_Optimized.hpp
 *  @brief GPU-accelerated sparse matrix operations with persistent GPU memory
 */
#ifndef GPU_SPARSE_MULT_OPTIMIZED_HPP
#define GPU_SPARSE_MULT_OPTIMIZED_HPP


#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

namespace NullGENO {

// GPU-resident sparse matrix descriptor
// Keeps I_mat on GPU to avoid repeated transfers
struct GPUSparseMatrix {
    int m, n;                          // matrix dimensions
    int nnz;                           // number of non-zeros
    float* d_values;                   // values on GPU
    int* d_rowPtr;                     // row pointers on GPU
    int* d_colInd;                     // column indices on GPU
    cusparseSpMatDescr_t matA;        // cuSPARSE descriptor
    cublasHandle_t blasHandle;         // cuBLAS handle
    cusparseHandle_t sparseHandle;     // cuSPARSE handle
    
    GPUSparseMatrix() : m(0), n(0), nnz(0), d_values(nullptr), 
                        d_rowPtr(nullptr), d_colInd(nullptr), 
                        matA(nullptr), blasHandle(nullptr), sparseHandle(nullptr) {}
};

// Global GPU sparse matrix (persistent across PCG iterations)
extern GPUSparseMatrix g_gpu_I_matrix;

// Initialize GPU with I_mat (call once at start of variance ratio estimation)
void initGPUSparseMatrix(const arma::sp_fmat& I_mat);

// Cleanup GPU memory (call at end of variance ratio estimation)
void cleanupGPUSparseMatrix();

// GPU-accelerated I matrix group summation using cuSPARSE (PERSISTENT VERSION)
// Computes: resultVec = I_mat * I_mat.t() * bVec using cuSPARSE
// I_mat stays on GPU across multiple iterations
arma::fvec getprodImatImattbVec_GPU_Sparse(arma::fvec & bVec);

} // namespace NullGENO

#endif // GPU_SPARSE_MULT_OPTIMIZED_HPP
