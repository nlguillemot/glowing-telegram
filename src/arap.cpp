#include "arap.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// For generating images of matrices for debugging
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// yanked in from the disgusting SVD headers
#include <algorithm>

// Intel Math Kernel Library
#include <mkl.h>

// --Incredibly-- over-engineered SVD implementation.
// Why is there no simple working svd3 library? This is ridiculous.
#define USE_SCALAR_IMPLEMENTATION
#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX
#include "Singular_Value_Decomposition_Preamble.hpp"

struct arap_system
{
    int nv;
    int nfree;

    // mapping between vertex IDs and rows of the system matrix ("free" aka unconstrained vertices)
    int* free_to_vid;
    int* vid_to_free;

    // storage for solver
    float* Ris;
    float* b;

    // factorized system matrix, compresessed dense lower triangular storage.
    float* F;
};

// outputs a tga image to visualize the (packed) matrix. useful for debugging.
static void matrix_to_image(const char* filename, int nv, const float* F)
{
    unsigned char* img = (unsigned char*)malloc(nv * nv * 4);
    for (int i = 0; i < nv; i++)
    {
        for (int j = 0; j < nv; j++)
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
            int c = int(roundf(F[lin])) + (1 << 23);
            img[(i * nv + j) * 4 + 0] = (c & 0xFF); // R
            img[(i * nv + j) * 4 + 1] = (c & 0xFF00) >> 8; // G
            img[(i * nv + j) * 4 + 2] = (c & 0xFF0000) >> 16; // B
            img[(i * nv + j) * 4 + 3] = 255; // A
        }
    }
    stbi_write_png(filename, nv, nv, 4, img, nv * 4);
    free(img);
}

arap_system* create_arap_system_matrix(
    int nv,
    const int* vcs,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws)
{
    arap_system* sys = (arap_system*)malloc(sizeof(arap_system));
    sys->nv = nv;

    // count the number of free vertices
    sys->nfree = 0;
    for (int i = 0; i < sys->nv; i++)
    {
        if (!vcs[i])
            sys->nfree++;
    }

    // list of free vertices (free vertex index -> vertexID)
    sys->free_to_vid = (int*)malloc(sizeof(int) * sys->nfree);
    
    // reverse mapping (vertexID -> free vertex index)
    sys->vid_to_free = (int*)malloc(sizeof(int) * sys->nv);

    for (int i = 0, curr_free = 0; i < sys->nv; i++)
    {
        if (vcs[i])
        {
            sys->vid_to_free[i] = -1;
            continue;
        }

        sys->free_to_vid[curr_free] = i;
        sys->vid_to_free[i] = curr_free;
        curr_free++;
    }

    // rotation matrices computed at each iteration
    sys->Ris = (float*)malloc(sizeof(float) * 9 * sys->nv);

    // storage for solver
    sys->b = (float*)malloc(sizeof(float) * 3 * sys->nfree);

    // number of weights in the system
    int nw = sys->nfree * (sys->nfree + 1) / 2;

    // matrix for weights that will be factorized
    sys->F = (float*)malloc(nw * sizeof(float));

    // initialize packed weight matrix
    memset(sys->F, 0, nw * sizeof(float));

    // initialize the rows for every free vertex in the laplacian
    for (int free_i = 0; free_i < sys->nfree; free_i++)
    {
        // the VertexID of this free vertex
        int vi = sys->free_to_vid[free_i];

        // go over every neighbor and add in their weights
        int curr_hID = v_hIDs[vi];
        do {
            // VertexID of neighbor
            int vj = h_vfnpIDs[4 * curr_hID + 0];

            // weight between these two vertices
            float wij = e_ws[curr_hID / 2];

            // constrained vertices don't take part in the laplacian matrix
            // this is accounted for later by fudging the rhs before solving
            if (!vcs[vj])
            {
                // compute linear index for packed matrix
                int row = free_i;
                int col = sys->vid_to_free[vj];
             
                // only set values in the row of this vertex, since this is building the matrix row-by-row.
                if (col < row)
                {
                    // set the weight between these two vertices.
                    // it's negative because rows of a laplacian matrix add to 0 (the diagonal is positive)
                    int lin = row + col * (2 * sys->nfree - col - 1) / 2;
                    sys->F[lin] = -wij;
                }
            }

            // add weight to diagonal
            int diag = free_i + free_i * (2 * sys->nfree - free_i - 1) / 2;
            sys->F[diag] += wij;

            curr_hID = h_vfnpIDs[4 * (curr_hID ^ 1) + 2];
        } while (curr_hID != v_hIDs[vi]);
    }


    // Output image of matrix (for debugging)
    matrix_to_image("arap_factorize_system_1.png", sys->nfree, sys->F);
    
    // cholesky factorize
    {
        char uplo = 'L';
        int n = sys->nfree;
        float* ap = sys->F;
        int info;
        spptrf_(&uplo, &n, ap, &info);
        assert(info == 0);
    }

    // Output image of matrix (for debugging)
    matrix_to_image("arap_factorize_system_2.png", sys->nfree, sys->F);

    return sys;
}

void destroy_arap_system_matrix(arap_system* sys)
{
    if (sys)
    {
        free(sys->free_to_vid);
        free(sys->vid_to_free);
        
        free(sys->Ris);
        free(sys->b);

        free(sys->F);
        
        free(sys);
    }
}

static void update_rotations(
    arap_system* sys,
    const float* p_bind_XYZs, const float* p_guess_XYZs,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws)
{
    // Compute new rigid rotation for every vertex
    for (int vi = 0; vi < sys->nv; vi++)
    {
        // the covariance matrix of all outgoing edges in the 1-ring
        float Si[9] = {};

        int curr_hID = v_hIDs[vi];
        do {
            int vj = h_vfnpIDs[4 * curr_hID + 0];
            float wij = e_ws[curr_hID / 2];
            
            // wij * eij = wij * (pi - pj)
            const float* pi_bind = &p_bind_XYZs[vi * 3];
            const float* pj_bind = &p_bind_XYZs[vj * 3];
            float weij[3] = {
                wij * (pi_bind[0] - pj_bind[0]),
                wij * (pi_bind[1] - pj_bind[1]),
                wij * (pi_bind[2] - pj_bind[2]),
            };
            
            // eij' = (pi' - pj')
            const float* pi_guess = &p_guess_XYZs[vi * 3];
            const float* pj_guess = &p_guess_XYZs[vj * 3];
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
        } while (curr_hID != v_hIDs[vi]);

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

        float Ui[9];
        Ui[0] = Su11.f;
        Ui[1] = Su21.f;
        Ui[2] = Su31.f;
        Ui[3] = Su12.f;
        Ui[4] = Su22.f;
        Ui[5] = Su32.f;
        Ui[6] = Su13.f;
        Ui[7] = Su23.f;
        Ui[8] = Su33.f;

        float Vi[9];
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
        float* Ri = &sys->Ris[vi * 9];
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
    arap_system* sys,
    const float* p_bind_XYZs, 
    float* p_guess_XYZs,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws)
{
    // column major, three columns of size n
    float* bx = &sys->b[sys->nfree * 0];
    float* by = &sys->b[sys->nfree * 1];
    float* bz = &sys->b[sys->nfree * 2];
    
    // initialize rhs
    for (int free_i = 0; free_i < sys->nfree; free_i++)
    {
        // get VertexID
        int vi = sys->free_to_vid[free_i];

        // summation of neighbors rhs
        float bx_tmp = 0.0f;
        float by_tmp = 0.0f;
        float bz_tmp = 0.0f;

        int curr_hID = v_hIDs[vi];
        do {
            int vj = h_vfnpIDs[4 * curr_hID + 0];
            float wij = e_ws[curr_hID / 2];

            const float* Ri = &sys->Ris[vi * 9];
            const float* Rj = &sys->Ris[vj * 9];
            const float* pi_bind = &p_bind_XYZs[vi * 3];
            const float* pj_bind = &p_bind_XYZs[vj * 3];

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

            if (sys->vid_to_free[vj] == -1)
            {
                // this neighbor is constrained, which means the rhs needs to be fudged
                // this is because the columns of the constrained vertices are set to zero in the system matrix,
                // which is necessary in order to keep the system matrix symmetric.
                // To counter this, the rhs is added to with the values that were originally erased.
                
                const float* pj_guess = &p_guess_XYZs[vj * 3];

                rhs[0] += wij * pj_guess[0];
                rhs[1] += wij * pj_guess[1];
                rhs[2] += wij * pj_guess[2];
            }

            // sum the neighbor's contribution to the rhs
            bx_tmp += rhs[0];
            by_tmp += rhs[1];
            bz_tmp += rhs[2];

            curr_hID = h_vfnpIDs[4 * (curr_hID ^ 1) + 2];
        } while (curr_hID != v_hIDs[vi]);

        // writeback to matrix
        bx[free_i] = bx_tmp;
        by[free_i] = by_tmp;
        bz[free_i] = bz_tmp;
    }

    // cholesky solve
    {
        char uplo = 'L';
        int n = sys->nfree;
        int nrhs = 3;
        float* ap = (float*)sys->F;
        int ldb = sys->nfree;
        int info;
        spptrs_(&uplo, &n, &nrhs, ap, sys->b, &ldb, &info);
        assert(info == 0);
    }

    // Update new guess
    for (int free_i = 0; free_i < sys->nfree; free_i++)
    {
        int vi = sys->free_to_vid[free_i];
        p_guess_XYZs[3 * vi + 0] = bx[free_i];
        p_guess_XYZs[3 * vi + 1] = by[free_i];
        p_guess_XYZs[3 * vi + 2] = bz[free_i];
    }
}

void arap(
    arap_system* sys, const float* p_bind_XYZs, float* p_guess_XYZs,
    const int* v_hIDs, const int* h_vfnpIDs, const float* e_ws, int ni)
{
    // iteratively refine guess by optimizing rotation and position
    for (int iter = 0; iter < ni; iter++)
    {
        update_rotations(
            sys, p_bind_XYZs, p_guess_XYZs, v_hIDs, h_vfnpIDs, e_ws);

        update_positions(
            sys, p_bind_XYZs, p_guess_XYZs, v_hIDs, h_vfnpIDs, e_ws);
    }
}
