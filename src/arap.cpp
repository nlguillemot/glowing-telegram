#include "arap.h"

// only svd() is used from this header
// could probably inline it all to keep this file standalone,
// or better yet come up with a solution that doesn't require a full svd.
#include "svd3.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

// grabs the Intel one if available, otherwise grabs the stub included with this project
#include <mkl.h>

// CLAPACK imports, if Intel MKL is not available.
#ifndef INTEL_MKL_VERSION
extern "C" int spptrf_(char* uplo, int* n, float* ap, int* info);
extern "C" int spptrs_(char *uplo, int* n, int* nrhs, float* ap, float* b, int* ldb, int* info);
#endif

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
    for (int c = 0; c < nc; c++)
    {
        // vertex to constrain == row to set to identity
        int vi = cs[c];

        // set the columns in the constrained row to identity
        for (int vj = 0; vj < vi; vj++)
        {
            int lin = vi + vj * (2 * nv - vj - 1) / 2;
            F[lin] = 0.0f;
        }

        // set the diagonal of the constrained row to identity
        int diag = vi + vi * (2 * nv - vi - 1) / 2;
        F[diag] = 1.0f;

        // set the rows of the constrained column to identity
        if (vi < nv - 1)
            memset(&F[diag + 1], 0, sizeof(float) * (nv - vi - 1));
    }

    // cholesky factorize
    char uplo = 'L';
    int n = nv;
    float* ap = F;
    int info;
    spptrf_(&uplo, &n, ap, &info);
    assert(info == 0);
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

        // compute SVD of covariance matrix
        float Ui[9], Sigma_i[9], Vi[9];
        svd(Si[0], Si[3], Si[6], Si[1], Si[4], Si[7], Si[2], Si[5], Si[8],
            Ui[0], Ui[3], Ui[6], Ui[1], Ui[4], Ui[7], Ui[2], Ui[5], Ui[8],
            Sigma_i[0], Sigma_i[3], Sigma_i[6], Sigma_i[1], Sigma_i[4], Sigma_i[7], Sigma_i[2], Sigma_i[5], Sigma_i[8],
            Vi[0], Vi[3], Vi[6], Vi[1], Vi[4], Vi[7], Vi[2], Vi[5], Vi[8]);
        
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

    // overwrite rhs with vertex constraints
    for (int c = 0; c < nc; c++)
    {
        // vertex to constrain
        int i = cs[c];
        const float* pi_guess = &p_guess_XYZs[i * 3];

        bx[i] = pi_guess[0];
        by[i] = pi_guess[1];
        bz[i] = pi_guess[2];
    }

    // for constraint "unit test"
    float constrained_test[3];
    if (nc > 0)
    {
        constrained_test[0] = p_guess_XYZs[cs[0] * 3 + 0];
        constrained_test[1] = p_guess_XYZs[cs[0] * 3 + 1];
        constrained_test[2] = p_guess_XYZs[cs[0] * 3 + 2];
    }

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

    // constraint "unit test"
    if (nc > 0)
    {
        assert(fabsf(constrained_test[0] - p_guess_XYZs[cs[0] * 3 + 0]) < 0.001f);
        assert(fabsf(constrained_test[1] - p_guess_XYZs[cs[0] * 3 + 1]) < 0.001f);
        assert(fabsf(constrained_test[2] - p_guess_XYZs[cs[0] * 3 + 2]) < 0.001f);
    }
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
