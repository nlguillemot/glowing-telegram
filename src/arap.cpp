#include "arap.h"

// only svd() is used from this header
// could probably inline it all to keep this file standalone,
// or better yet come up with a solution that doesn't require a full svd.
#include "svd3.h"

#include <assert.h>
#include <stdlib.h>

// CLAPACK imports
extern "C" int spptrf_(char* uplo, int* n, float* ap, int* info);
extern "C" int spptrs_(char *uplo, int* n, int* nrhs, float* ap, float* b, int* ldb, int* info);

void arap_factorize_system(
    int nv, int ne,
    const int* h_vfnpIDs,
    const float* e_ws,
    float* F)
{
    // number of packed weights
    int nw = nv * (nv + 1) / 2;

    // initialize packed weight matrix
    for (int w = 0; w < nw; w++)
    {
        F[w] = 0.0f;
    }

    // set the weights for each edge
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
        // two diagonals since each edge affects two rows (two vertices)
        int diag1 = vi + vi * (2 * nv - vi - 1) / 2;
        F[diag1] += wij;
        int diag2 = vj + vj * (2 * nv - vj - 1) / 2;
        F[diag2] += wij;
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
    for (int i = 0; i < nv; i++)
    {
        // build the covariance matrix of all outgoing edges in the 1-ring
        float Si[9] = {};
        int curr_hID = v_hIDs[i];
        do {
            int j = h_vfnpIDs[4 * curr_hID + 0];
            float wij = e_ws[curr_hID / 2];
            
            // wij * eij = wij * (pi - pj)
            const float* pi_bind_XYZ = &p_bind_XYZs[i * 3];
            const float* pj_bind_XYZ = &p_bind_XYZs[j * 3];
            float weij[3] = {
                wij * (pi_bind_XYZ[0] - pj_bind_XYZ[0]),
                wij * (pi_bind_XYZ[1] - pj_bind_XYZ[1]),
                wij * (pi_bind_XYZ[2] - pj_bind_XYZ[2]),
            };
            
            // eij' = (pi' - pj')
            const float* pi_guess_XYZ = &p_guess_XYZs[i * 3];
            const float* pj_guess_XYZ = &p_guess_XYZs[j * 3];
            float eij_guess[3] = {
                pi_guess_XYZ[0] - pj_guess_XYZ[0],
                pi_guess_XYZ[1] - pj_guess_XYZ[1],
                pi_guess_XYZ[2] - pj_guess_XYZ[2],
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
    }
}

static void update_positions(
    const float* p_bind_XYZs, float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    const float* F,
    int* ks, int nk,
    const float* Ris)
{
    // column major, three columns of size n
    float* b = (float*)malloc(3 * sizeof(float) * nv);
    float* bx = &b[nv * 0];
    float* by = &b[nv * 1];
    float* bz = &b[nv * 2];

    // initialize rhs
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
            const float* pi = &p_bind_XYZs[i * 3];
            const float* pj = &p_bind_XYZs[j * 3];

            // (pi - pj) / 2
            float pi_m_pj_2[3] = {
                (pi[0] - pj[0]) / 2.0f,
                (pi[1] - pj[1]) / 2.0f,
                (pi[2] - pj[2]) / 2.0f
            };

            // Ri + Rj
            float Ri_p_Rj[9] = {
                Ri[0] + Rj[0], Ri[1] + Rj[1], Ri[2] + Rj[2],
                Ri[3] + Rj[3], Ri[4] + Rj[4], Ri[5] + Rj[5],
                Ri[6] + Rj[6], Ri[7] + Rj[7], Ri[8] + Rj[8]
            };

            // 1/2 (Ri + Rj) (pi - pj)
            float rhs[3] = {
                Ri_p_Rj[0] * pi_m_pj_2[0] + Ri_p_Rj[3] * pi_m_pj_2[1] + Ri_p_Rj[6] * pi_m_pj_2[2],
                Ri_p_Rj[1] * pi_m_pj_2[0] + Ri_p_Rj[4] * pi_m_pj_2[1] + Ri_p_Rj[7] * pi_m_pj_2[2],
                Ri_p_Rj[2] * pi_m_pj_2[0] + Ri_p_Rj[5] * pi_m_pj_2[1] + Ri_p_Rj[8] * pi_m_pj_2[2]
            };

            // all neighbors are summed
            bx_tmp += wij * rhs[0];
            by_tmp += wij * rhs[1];
            bz_tmp += wij * rhs[2];

            curr_hID = h_vfnpIDs[4 * (curr_hID ^ 1) + 2];
        } while (curr_hID != v_hIDs[i]);

        // writeback to big ol' matrix
        bx[i] = bx_tmp;
        by[i] = by_tmp;
        bz[i] = bz_tmp;
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

    free(b);
}

void arap(
    const float* p_bind_XYZs, float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    const float* F,
    int* ks, int nk,
    int ni)
{
    // rotation matrices computed at each iteration
    float* Ris = (float*)malloc(sizeof(float) * 9 * nv);
    
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
            F,
            ks, nk,
            Ris);
    }

    free(Ris);
}