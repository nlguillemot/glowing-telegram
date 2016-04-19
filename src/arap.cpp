#include "arap.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// yanked in from the disgusting SVD headers
#include <algorithm>

// Used for debugging
#include <stdio.h>

// OpenMP for parallel loops
#include <omp.h>

// grabs the Intel one if available, otherwise grabs the stub included with this project
#include <mkl.h>

// CLAPACK imports, if Intel MKL is not available.
#ifndef INTEL_MKL_VERSION
extern "C" int spptrf_(char* uplo, int* n, float* ap, int* info);
extern "C" int spptrs_(char *uplo, int* n, int* nrhs, float* ap, float* b, int* ldb, int* info);
#endif

// --Incredibly-- over-engineered SVD implementation.
// Why is there no simple working svd3 library? This is ridiculous.
#define USE_SCALAR_IMPLEMENTATION
#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX
#include "Singular_Value_Decomposition_Preamble.hpp"

static void matrix_to_file(const char* filename, int nv, const float* F)
{
#ifdef _MSC_VER
    FILE* f = NULL;
    fopen_s(&f, filename, "wb");
#else
    FILE* f = fopen(filename, "wb");
#endif
    if (f)
    {
        int width = nv & 0xFFFF;
        int height = nv & 0xFFFF;

        // max image size = 2^16
        assert(width == nv);
        assert(height == nv);

        fputc(0, f);
        fputc(0, f);
        fputc(2, f); // uncompressed RGB
        fputc(0, f); fputc(0, f);
        fputc(0, f); fputc(0, f);
        fputc(0, f);
        fputc(0, f); fputc(0, f); // x origin
        fputc(0, f); fputc(0, f); // y origin
        fputc((width & 0x00FF), f);
        fputc((width & 0xFF00) / 256, f);
        fputc((height & 0x00FF), f);
        fputc((height & 0xFF00) / 256, f);
        fputc(32, f); // 32 bit bitmap
        fputc(1 << 5, f); // origin top left
        unsigned char* img = (unsigned char*)malloc(width * height * 4);
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                int vi = i;
                int vj = j;

                // enforce j <= i < n
                if (vi < vj)
                {
                    int tmp = vi;
                    vi = vj;
                    vj = tmp;
                }

                int lin = vi + vj * (2 * nv - vj - 1) / 2;
                int c = (int)ceilf(fabsf(F[lin]));
                if (c > 255)
                    c = 255;
                // img[(i * width + j) * 4 + 0] = (lin & 0xFF0000) >> 16; // B
                // img[(i * width + j) * 4 + 1] = (lin & 0xFF00) >> 8; // G
                // img[(i * width + j) * 4 + 2] = lin & 0xFF; // R
                img[(i * width + j) * 4 + 0] = c; // B
                img[(i * width + j) * 4 + 1] = c; // G
                img[(i * width + j) * 4 + 2] = c; // R
                                                  //if (F[lin] == 0.0f)
                                                  //{
                                                  //    img[(i * width + j) * 4 + 0] = 255; // B
                                                  //    img[(i * width + j) * 4 + 1] = 0; // G
                                                  //    img[(i * width + j) * 4 + 2] = 0; // R
                                                  //}
                img[(i * width + j) * 4 + 3] = 255; // A
            }
        }
        fwrite(img, width * height * 4, 1, f);
        free(img);
        fclose(f);
    }
}

void arap_factorize_system(
    int nv,
    const int* h_vfnpIDs,
    const float* e_ws, int ne,
    const int* cs, int nc,
    float* F)
{
    // number of packed weights
    int nw = nv * (nv + 1) / 2;

    // initialize packed weight matrix
    memset(F, 0, nw * sizeof(float));

    // set the weights for each edge in the laplacian matrix
    for (int h = 0; h < ne * 2; h += 2)
    {
        float wij = e_ws[h / 2];
        int vi = h_vfnpIDs[h *  4 + 0];
        int vj = h_vfnpIDs[(h + 1) * 4 + 0];

        // enforce j <= i < n
        if (vi < vj)
        {
            int tmp = vi;
            vi = vj;
            vj = tmp;
        }

        // convert to linear index and write into packed matrix
        int lin = vi + vj * (2 * nv - vj - 1) / 2;
        F[lin] = -wij;

        // diagonals ensure each row sum to zero (laplacian matrix)
        // affects two diagonals since each edge affects two rows (two vertices)
        int diag1 = vi + vi * (2 * nv - vi - 1) / 2;
        F[diag1] += wij;
        int diag2 = vj + vj * (2 * nv - vj - 1) / 2;
        F[diag2] += wij;
    }

    // implement constrained vertices by setting their row/col's weights to identity
    //for (int c = 0; c < nc; c++)
    //{
    //    // vertex to constrain == row/col to set to identity
    //    int vi = cs[c];

    //    // set the columns in the constrained row to identity
    //    for (int vj = 0; vj < vi; vj++)
    //    {
    //        int lin = vi + vj * (2 * nv - vj - 1) / 2;
    //        F[lin] = 0.0f;
    //    }

    //    // set the diagonal of the constrained row to identity
    //    int diag = vi + vi * (2 * nv - vi - 1) / 2;
    //    F[diag] = 1.0f;

    //    // set the rows of the constrained column to identity
    //    // This is actually wrong, mathematically speaking.
    //    // Only the row of vi should be set to zero.
    //    // However, that creates a non-symmetric matrix, which is less efficient.
    //    // In order to make up for this hack, the rhs of the equation will be offsetted before solving.
    //    if (vi < nv - 1)
    //        memset(&F[diag + 1], 0, sizeof(float) * (nv - vi - 1));
    //}
    
    // Output image of matrix (for debugging)
    matrix_to_file("arap_factorize_system_1.tga", nv, F);

    // cholesky factorize
    char uplo = 'L';
    int n = nv;
    float* ap = F;
    int info;
    spptrf_(&uplo, &n, ap, &info);
    assert(info == 0);

    // Output image of matrix (for debugging)
    matrix_to_file("arap_factorize_system_2.tga", nv, F);
}

static void update_rotations(
    const float* p_bind_XYZs, const float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    float* Ris)
{
    // Compute new rigid rotation for every vertex
#pragma omp parallel for
    for (int i = 0; i < nv; i++)
    {
        // build the covariance matrix of all outgoing edges in the 1-ring
        float Si[9] = {};
        int curr_hID = v_hIDs[i];
        do {
            int j = h_vfnpIDs[4 * curr_hID + 0];
            float wij = e_ws[curr_hID / 2];
            
            // wij * eij = wij * (pi - pj)
            const float* pi_bind = &p_bind_XYZs[i * 3];
            const float* pj_bind = &p_bind_XYZs[j * 3];
            float weij[3] = {
                wij * (pi_bind[0] - pj_bind[0]),
                wij * (pi_bind[1] - pj_bind[1]),
                wij * (pi_bind[2] - pj_bind[2]),
            };
            
            // eij' = (pi' - pj')
            const float* pi_guess = &p_guess_XYZs[i * 3];
            const float* pj_guess = &p_guess_XYZs[j * 3];
            float eij_guess[3] = {
                pi_guess[0] - pj_guess[0],
                pi_guess[1] - pj_guess[1],
                pi_guess[2] - pj_guess[2],
            };

            // Si += weij * eij'^T
            Si[0] += weij[0] * eij_guess[0];
            Si[1] += weij[1] * eij_guess[0];
            Si[2] += weij[2] * eij_guess[0];
            Si[3] += weij[0] * eij_guess[1];
            Si[4] += weij[1] * eij_guess[1];
            Si[5] += weij[2] * eij_guess[1];
            Si[6] += weij[0] * eij_guess[2];
            Si[7] += weij[1] * eij_guess[2];
            Si[8] += weij[2] * eij_guess[2];
            
            curr_hID = h_vfnpIDs[4 * (curr_hID ^ 1) + 2];
        } while (curr_hID != v_hIDs[i]);

        // compute SVD of covariance matrix, while attempting to contain how disgusting this svd library is.
        float test_A11, test_A21, test_A31, test_A12, test_A22, test_A32, test_A13, test_A23, test_A33;
        test_A11 = Si[0];
        test_A21 = Si[1];
        test_A31 = Si[2];
        test_A12 = Si[3];
        test_A22 = Si[4];
        test_A32 = Si[5];
        test_A13 = Si[6];
        test_A23 = Si[7];
        test_A33 = Si[8];

#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"

        ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=test_A11;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=test_A21;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=test_A31;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=test_A12;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=test_A22;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=test_A32;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=test_A13;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=test_A23;) 
        ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=test_A33;) 

#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"

        float Ui[9], Vi[9];
        Ui[0] = Su11.f;
        Ui[1] = Su21.f;
        Ui[2] = Su31.f;
        Ui[3] = Su12.f;
        Ui[4] = Su22.f;
        Ui[5] = Su32.f;
        Ui[6] = Su13.f;
        Ui[7] = Su23.f;
        Ui[8] = Su33.f;

        Vi[0] = Sv11.f;
        Vi[1] = Sv21.f;
        Vi[2] = Sv31.f;
        Vi[3] = Sv12.f;
        Vi[4] = Sv22.f;
        Vi[5] = Sv32.f;
        Vi[6] = Sv13.f;
        Vi[7] = Sv23.f;
        Vi[8] = Sv33.f;

        // Ri = Vi * Ui^T
        float* Ri = &Ris[i * 9];
        Ri[0] = Vi[0] * Ui[0] + Vi[3] * Ui[3] + Vi[6] * Ui[6];
        Ri[1] = Vi[1] * Ui[0] + Vi[4] * Ui[3] + Vi[7] * Ui[6];
        Ri[2] = Vi[2] * Ui[0] + Vi[5] * Ui[3] + Vi[8] * Ui[6];
        Ri[3] = Vi[0] * Ui[1] + Vi[3] * Ui[4] + Vi[6] * Ui[7];
        Ri[4] = Vi[1] * Ui[1] + Vi[4] * Ui[4] + Vi[7] * Ui[7];
        Ri[5] = Vi[2] * Ui[1] + Vi[5] * Ui[4] + Vi[8] * Ui[7];
        Ri[6] = Vi[0] * Ui[2] + Vi[3] * Ui[5] + Vi[6] * Ui[8];
        Ri[7] = Vi[1] * Ui[2] + Vi[4] * Ui[5] + Vi[7] * Ui[8];
        Ri[8] = Vi[2] * Ui[2] + Vi[5] * Ui[5] + Vi[8] * Ui[8];

        // check for and handle reflection case
        float det = Ri[0] * Ri[4] * Ri[8] + Ri[3] * Ri[7] * Ri[2] + Ri[6] * Ri[1] * Ri[5] - Ri[2] * Ri[4] * Ri[6] - Ri[5] * Ri[7] * Ri[0] - Ri[8] * Ri[1] * Ri[3];
        if (det < 0.0f)
        {
            Ri[6] *= -1.0f;
            Ri[7] *= -1.0f;
            Ri[8] *= -1.0f;
        }
    }
}

static void update_positions(
    const float* p_bind_XYZs, float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    const float* F, float* b,
    const int* cs, int nc,
    const float* Ris)
{
    // column major, three columns of size n
    float* bx = &b[nv * 0];
    float* by = &b[nv * 1];
    float* bz = &b[nv * 2];

    // initialize rhs
#pragma omp parallel for
    for (int i = 0; i < nv; i++)
    {
        // summation of neighbors rhs
        float bx_tmp = 0.0f;
        float by_tmp = 0.0f;
        float bz_tmp = 0.0f;

        int curr_hID = v_hIDs[i];
        do {
            int j = h_vfnpIDs[4 * curr_hID + 0];
            float wij = e_ws[curr_hID / 2];

            const float* Ri = &Ris[i * 9];
            const float* Rj = &Ris[j * 9];
            const float* pi_bind = &p_bind_XYZs[i * 3];
            const float* pj_bind = &p_bind_XYZs[j * 3];

            // pi - pj
            float pi_pj[3] = {
                pi_bind[0] - pj_bind[0],
                pi_bind[1] - pj_bind[1],
                pi_bind[2] - pj_bind[2]
            };

            // Ri + Rj
            float Ri_Rj[9] = {
                Ri[0] + Rj[0], Ri[1] + Rj[1], Ri[2] + Rj[2],
                Ri[3] + Rj[3], Ri[4] + Rj[4], Ri[5] + Rj[5],
                Ri[6] + Rj[6], Ri[7] + Rj[7], Ri[8] + Rj[8]
            };

            // wij 1/2 (Ri + Rj) (pi - pj)
            float wij_2 = wij / 2.0f;
            float rhs[3] = {
                wij_2 * (Ri_Rj[0] * pi_pj[0] + Ri_Rj[3] * pi_pj[1] + Ri_Rj[6] * pi_pj[2]),
                wij_2 * (Ri_Rj[1] * pi_pj[0] + Ri_Rj[4] * pi_pj[1] + Ri_Rj[7] * pi_pj[2]),
                wij_2 * (Ri_Rj[2] * pi_pj[0] + Ri_Rj[5] * pi_pj[1] + Ri_Rj[8] * pi_pj[2])
            };

            // all neighbors are summed
            bx_tmp += rhs[0];
            by_tmp += rhs[1];
            bz_tmp += rhs[2];

            curr_hID = h_vfnpIDs[4 * (curr_hID ^ 1) + 2];
        } while (curr_hID != v_hIDs[i]);

        // writeback to matrix
        bx[i] = bx_tmp;
        by[i] = by_tmp;
        bz[i] = bz_tmp;
    }

    // fudge the rhs, as discussed in arap_factorize_system
    // again, this is because the columns of the constrained vertices are set to zero in the system matrix.
    // to counter this, the rhs is added to with the values that were originally erased.
    // this is done in order to make the matrix symmetric, since the constraint rows are not symmetric otherwise.
    //for (int c = 0; c < nc; c++)
    //{
    //    // vertex to constrain
    //    int i = cs[c];
    //    const float* pi_guess = &p_guess_XYZs[i * 3];

    //    int curr_hID = v_hIDs[i];
    //    do {
    //        int j = h_vfnpIDs[4 * curr_hID + 0];
    //        float wij = e_ws[curr_hID / 2];

    //        bx[j] += wij * pi_guess[0];
    //        by[j] += wij * pi_guess[1];
    //        bz[j] += wij * pi_guess[2];

    //    } while (curr_hID != v_hIDs[i]);
    //}
    
    // enforce the constraints on the vertices themselves
    // done in a second pass to handle constrained vertices that are neighbors to constrained vertices.
    //for (int c = 0; c < nc; c++)
    //{
    //    // vertex to constrain
    //    int i = cs[c];
    //    const float* pi_guess = &p_guess_XYZs[i * 3];

    //    // enforce constraint
    //    bx[i] = pi_guess[0];
    //    by[i] = pi_guess[1];
    //    bz[i] = pi_guess[2];
    //}

    // for constraint "unit test"
    //float constrained_test[3];
    //if (nc > 0)
    //{
    //    constrained_test[0] = p_guess_XYZs[cs[0] * 3 + 0];
    //    constrained_test[1] = p_guess_XYZs[cs[0] * 3 + 1];
    //    constrained_test[2] = p_guess_XYZs[cs[0] * 3 + 2];
    //}

    // cholesky solve
    char uplo = 'L';
    int n = nv;
    int nrhs = 3;
    float* ap = (float*)F;
    int ldb = nv;
    int info;
    spptrs_(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    assert(info == 0);

    // Update new guess
    for (int i = 0; i < nv; i++)
    {
        p_guess_XYZs[3 * i + 0] = bx[i];
        p_guess_XYZs[3 * i + 1] = by[i];
        p_guess_XYZs[3 * i + 2] = bz[i];
    }

    //// constraint "unit test"
    //if (nc > 0)
    //{
    //    assert(fabsf(constrained_test[0] - p_guess_XYZs[cs[0] * 3 + 0]) < 0.001f);
    //    assert(fabsf(constrained_test[1] - p_guess_XYZs[cs[0] * 3 + 1]) < 0.001f);
    //    assert(fabsf(constrained_test[2] - p_guess_XYZs[cs[0] * 3 + 2]) < 0.001f);
    //}
}

void arap(
    const float* p_bind_XYZs, float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    const float* F,
    const int* cs, int nc,
    int ni)
{
    // rotation matrices computed at each iteration
    float* Ris = (float*)malloc(sizeof(float) * 9 * nv);

    // storage for solver
    float* b = (float*)malloc(sizeof(float) * 3 * nv);

    // iteratively refine guess by optimizing rotation and position
    for (int iter = 0; iter < ni; iter++)
    {
        update_rotations(
            p_bind_XYZs, p_guess_XYZs, nv,
            v_hIDs,
            h_vfnpIDs, e_ws,
            Ris);

        update_positions(
            p_bind_XYZs, p_guess_XYZs, nv,
            v_hIDs,
            h_vfnpIDs, e_ws,
            F, b,
            cs, nc,
            Ris);
    }
    
    free(b);
    free(Ris);
}
