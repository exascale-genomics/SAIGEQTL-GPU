#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <RcppArmadillo.h>
#include "Rcpp.h"

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
GPUSparseMatrix g_gpu_I_matrix;

// Initialize GPU with I_mat (call ONCE at start of variance ratio estimation)
void initGPUSparseMatrix(const arma::sp_fmat& I_mat) {
    std::cout << "[GPU-PERSISTENT] Initializing GPU sparse matrix (ONCE per variance ratio estimation)" << std::endl;
    std::cout << "[GPU-PERSISTENT] I_mat: " << I_mat.n_rows << " x " << I_mat.n_cols 
              << ", nnz=" << I_mat.n_nonzero << std::endl;
    
    g_gpu_I_matrix.m = I_mat.n_rows;
    g_gpu_I_matrix.n = I_mat.n_cols;
    g_gpu_I_matrix.nnz = I_mat.n_nonzero;
    
    if (g_gpu_I_matrix.nnz == 0) {
        std::cout << "[GPU-PERSISTENT] Empty matrix, skipping GPU initialization" << std::endl;
        return;
    }
    
    // Create handles
    cusparseCreate(&g_gpu_I_matrix.sparseHandle);
    cublasCreate(&g_gpu_I_matrix.blasHandle);
    
    // 1. Extract COO format from armadillo sparse matrix
    const float* values = I_mat.values;
    const arma::uword* row_indices = I_mat.row_indices;
    const arma::uword* col_ptrs = I_mat.col_ptrs;
    
    // Convert COO to CSR format
    std::vector<int> rowPtr(g_gpu_I_matrix.m + 1, 0);
    std::vector<int> colInd(g_gpu_I_matrix.nnz);
    std::vector<float> csrValues(g_gpu_I_matrix.nnz);
    
    int nnz_count = 0;
    for (int j = 0; j < g_gpu_I_matrix.n; j++) {
        for (arma::uword k = col_ptrs[j]; k < col_ptrs[j + 1]; k++) {
            int i = (int)row_indices[k];
            rowPtr[i + 1]++;
            colInd[nnz_count] = j;
            csrValues[nnz_count] = values[k];
            nnz_count++;
        }
    }
    
    // Cumulative sum for rowPtr
    for (int i = 1; i <= g_gpu_I_matrix.m; i++) {
        rowPtr[i] += rowPtr[i - 1];
    }
    
    // 2. Allocate GPU memory
    cudaMalloc((void**)&g_gpu_I_matrix.d_values, g_gpu_I_matrix.nnz * sizeof(float));
    cudaMalloc((void**)&g_gpu_I_matrix.d_rowPtr, (g_gpu_I_matrix.m + 1) * sizeof(int));
    cudaMalloc((void**)&g_gpu_I_matrix.d_colInd, g_gpu_I_matrix.nnz * sizeof(int));
    
    // 3. Copy data to GPU (ONE TIME!)
    cudaMemcpy(g_gpu_I_matrix.d_values, csrValues.data(), g_gpu_I_matrix.nnz * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(g_gpu_I_matrix.d_rowPtr, rowPtr.data(), (g_gpu_I_matrix.m + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(g_gpu_I_matrix.d_colInd, colInd.data(), g_gpu_I_matrix.nnz * sizeof(int), cudaMemcpyHostToDevice);
    
    // 4. Create cuSPARSE descriptor
    cusparseCreateCsr(&g_gpu_I_matrix.matA, g_gpu_I_matrix.m, g_gpu_I_matrix.n, g_gpu_I_matrix.nnz,
                      g_gpu_I_matrix.d_rowPtr, g_gpu_I_matrix.d_colInd, g_gpu_I_matrix.d_values,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
    
    std::cout << "[GPU-PERSISTENT] GPU sparse matrix initialized and transferred (will be reused)" << std::endl;
}

// Cleanup GPU memory (call at end of variance ratio estimation)
void cleanupGPUSparseMatrix() {
    std::cout << "[GPU-PERSISTENT] Cleaning up GPU sparse matrix" << std::endl;
    
    if (g_gpu_I_matrix.matA != nullptr) {
        cusparseDestroySpMat(g_gpu_I_matrix.matA);
    }
    if (g_gpu_I_matrix.blasHandle != nullptr) {
        cublasDestroy(g_gpu_I_matrix.blasHandle);
    }
    if (g_gpu_I_matrix.sparseHandle != nullptr) {
        cusparseDestroy(g_gpu_I_matrix.sparseHandle);
    }
    if (g_gpu_I_matrix.d_values != nullptr) {
        cudaFree(g_gpu_I_matrix.d_values);
    }
    if (g_gpu_I_matrix.d_rowPtr != nullptr) {
        cudaFree(g_gpu_I_matrix.d_rowPtr);
    }
    if (g_gpu_I_matrix.d_colInd != nullptr) {
        cudaFree(g_gpu_I_matrix.d_colInd);
    }
    
    g_gpu_I_matrix = GPUSparseMatrix();  // Reset
    std::cout << "[GPU-PERSISTENT] GPU cleanup complete" << std::endl;
}

// GPU-accelerated I matrix group summation using cuSPARSE (PERSISTENT VERSION)
// Computes: resultVec = I_mat * I_mat.t() * bVec using cuSPARSE
// ONLY transfers vectors, matrix stays on GPU!
arma::fvec getprodImatImattbVec_GPU_Sparse(arma::fvec & bVec) {
    int n = bVec.n_elem;
    
    std::cout << "[GPU-PERSISTENT] getprodImatImattbVec_GPU_Sparse: n=" << n << std::endl;
    
    if (g_gpu_I_matrix.nnz == 0) {
        std::cout << "[GPU-PERSISTENT] GPU matrix not initialized, returning zeros" << std::endl;
        return arma::fvec(n, arma::fill::zeros);
    }
    
    float* d_bVec, *d_temp, *d_result;
    
    // 1. Allocate GPU memory for VECTORS ONLY (not matrix!)
    cudaMalloc((void**)&d_bVec, n * sizeof(float));
    cudaMalloc((void**)&d_temp, g_gpu_I_matrix.n * sizeof(float));
    cudaMalloc((void**)&d_result, n * sizeof(float));
    
    std::cout << "[GPU-PERSISTENT] GPU vector memory allocated" << std::endl;
    
    // 2. Transfer input vector to GPU
    cudaMemcpy(d_bVec, bVec.memptr(), n * sizeof(float), cudaMemcpyHostToDevice);
    
    std::cout << "[GPU-PERSISTENT] Input vector transferred to GPU (only this is moved per iteration)" << std::endl;
    
    float alpha = 1.0f, beta = 0.0f;
    
    // Create dense vector descriptors
    cusparseDnVecDescr_t vecB, vecTemp, vecResult;
    cusparseCreateDnVec(&vecB, n, d_bVec, CUDA_R_32F);
    cusparseCreateDnVec(&vecTemp, g_gpu_I_matrix.n, d_temp, CUDA_R_32F);
    cusparseCreateDnVec(&vecResult, n, d_result, CUDA_R_32F);
    
    // 3. First multiplication: temp = I_mat.t() * bVec
    // This requires transpose, but cuSPARSE doesn't directly support transpose
    // We compute: I_mat.t() * bVec using CUSPARSE_OPERATION_TRANSPOSE
    std::cout << "[GPU-PERSISTENT] Starting SpMV (I_mat.t() * bVec)" << std::endl;
    
    size_t bufferSize = 0;
    cusparseSpMV_bufferSize(g_gpu_I_matrix.sparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
                            &alpha, g_gpu_I_matrix.matA, vecB, &beta, vecTemp,
                            CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
    
    void* d_buffer;
    cudaMalloc(&d_buffer, bufferSize);
    
    cusparseSpMV(g_gpu_I_matrix.sparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
                 &alpha, g_gpu_I_matrix.matA, vecB, &beta, vecTemp,
                 CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, d_buffer);
    
    std::cout << "[GPU-PERSISTENT] First SpMV complete" << std::endl;
    
    // 4. Second multiplication: result = I_mat * temp (using cuBLAS for dense-sparse)
    // Actually use SpMV again with non-transpose for I_mat * temp
    std::cout << "[GPU-PERSISTENT] Starting SpMV (I_mat * temp)" << std::endl;
    
    cusparseSpMV(g_gpu_I_matrix.sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                 &alpha, g_gpu_I_matrix.matA, vecTemp, &beta, vecResult,
                 CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, d_buffer);
    
    std::cout << "[GPU-PERSISTENT] Second SpMV complete" << std::endl;
    
    // 5. Copy result back to CPU
    arma::fvec resultVec(n);
    cudaMemcpy(resultVec.memptr(), d_result, n * sizeof(float), cudaMemcpyDeviceToHost);
    
    std::cout << "[GPU-PERSISTENT] Result transferred back to CPU" << std::endl;
    
    // 6. Cleanup temporary vectors (not the persistent matrix!)
    cusparseDestroyDnVec(vecB);
    cusparseDestroyDnVec(vecTemp);
    cusparseDestroyDnVec(vecResult);
    
    cudaFree(d_bVec);
    cudaFree(d_temp);
    cudaFree(d_result);
    cudaFree(d_buffer);
    
    std::cout << "[GPU-PERSISTENT] PCG iteration complete" << std::endl;
    
    return resultVec;
}

} // namespace NullGENO
