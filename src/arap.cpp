#include "arap.h"

#include "svd3.h"

struct vec3
{
    union
    {
        struct
        {
            float x, y, z;
        };
        float e[3];
    };

    float& operator[](int i)
    {
        return e[i];
    }
};

vec3 operator-(vec3 a, vec3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

vec3 operator*(float s, vec3 v)
{
    v.x = s * v.x;
    v.y = s * v.y;
    v.z = s * v.z;
    return v;
}

struct mat3
{
    float m[9];

    mat3& operator+=(const mat3& other)
    {
        m[0] += other.m[0];
        m[1] += other.m[1];
        m[2] += other.m[2];
        m[3] += other.m[3];
        m[4] += other.m[4];
        m[5] += other.m[5];
        m[6] += other.m[6];
        m[7] += other.m[7];
        m[8] += other.m[8];
        return *this;
    }
};

mat3 transpose(const mat3& other)
{
    mat3 m;
    m.m[0] = other.m[0];
    m.m[1] = other.m[3];
    m.m[2] = other.m[6];
    m.m[3] = other.m[1];
    m.m[4] = other.m[4];
    m.m[5] = other.m[7];
    m.m[6] = other.m[2];
    m.m[7] = other.m[5];
    m.m[8] = other.m[8];
    return m;
}

mat3 operator*(const mat3& a, const mat3& b)
{
    mat3 m = {};
    int i = 0;
    for (int b_col = 0; b_col < 3; b_col++)
    {
        for (int a_row = 0; a_row < 3; a_row++)
        {
            m.m[i] += a.m[0 + a_row] * b.m[b_col * 3 + 0];
            m.m[i] += a.m[3 + a_row] * b.m[b_col * 3 + 1];
            m.m[i] += a.m[6 + a_row] * b.m[b_col * 3 + 2];
            i++;
        }
    }
    return m;
}

mat3 outerProduct(vec3 a, vec3 b)
{
    mat3 m;
    m.m[0] = a[0] * b[0];
    m.m[1] = a[1] * b[0];
    m.m[2] = a[2] * b[0];
    m.m[3] = a[0] * b[1];
    m.m[4] = a[1] * b[1];
    m.m[5] = a[2] * b[1];
    m.m[6] = a[0] * b[2];
    m.m[7] = a[1] * b[2];
    m.m[8] = a[2] * b[2];
    return m;
}

struct halfedge
{
    int VertexID;
    int FaceID;
    int NextHalfedgeID;
    int PrevHalfedgeID;
};

static void update_rotations(
    const vec3* ps_bind, const vec3* ps_guess, int nv,
    const int* v_hIDs,
    const halfedge* hs,
    const float* e_ws,
    mat3* Ris)
{
    for (int i = 0; i < nv; i++)
    {
        // build the covariance matrix of all outgoing edges in the 1-ring
        mat3 Si = {};
        int currHalfedgeID = v_hIDs[i];
        do {
            int j = hs[currHalfedgeID].VertexID;
            float wij = e_ws[currHalfedgeID / 2];
            vec3 weij = wij * (ps_bind[i] - ps_bind[j]);
            vec3 eij_guess = (ps_guess[i] - ps_guess[j]);

            Si += outerProduct(weij, eij_guess);
            
            currHalfedgeID = hs[currHalfedgeID ^ 1].NextHalfedgeID;
        } while (currHalfedgeID != v_hIDs[i]);

        // compute SVD of covariance matrix
        mat3 Ui, Sigma_i, Vi;
        svd(Si.m[0], Si.m[3], Si.m[6], Si.m[1], Si.m[4], Si.m[7], Si.m[2], Si.m[5], Si.m[8],
            Ui.m[0], Ui.m[3], Ui.m[6], Ui.m[1], Ui.m[4], Ui.m[7], Ui.m[2], Ui.m[5], Ui.m[8],
            Sigma_i.m[0], Sigma_i.m[3], Sigma_i.m[6], Sigma_i.m[1], Sigma_i.m[4], Sigma_i.m[7], Sigma_i.m[2], Sigma_i.m[5], Sigma_i.m[8],
            Vi.m[0], Vi.m[3], Vi.m[6], Vi.m[1], Vi.m[4], Vi.m[7], Vi.m[2], Vi.m[5], Vi.m[8]);
        
        Ris[i] = Vi * transpose(Ui);
    }
}

static void update_positions(
    const vec3* ps_bind, vec3* ps_guess, int nv,
    int* ks, int nk,
    const mat3* Ris)
{

}

void arap(
    const float* p_bind_XYZs, float* p_guess_XYZs, int nv,
    const int* v_hIDs,
    const int* h_vfnpIDs,
    const float* e_ws,
    int* ks, int nk,
    int ni)
{
    // TODO: Need pre-computed system matrix as input
    // bind pose and the current guess iteration
    const vec3* ps_bind = (const vec3*)p_bind_XYZs;
    vec3* ps_guess = (vec3*)p_guess_XYZs;
    const halfedge* hs = (const halfedge*)h_vfnpIDs;

    // rotation matrices computed at each iteration
    mat3* Ris = new mat3[nv];
    
    // iteratively refine guess by optimizing rotation and position
    for (int iter = 0; iter < ni; iter++)
    {
        update_rotations(
            ps_bind, ps_guess, nv,
            v_hIDs,
            hs, e_ws,
            Ris);

        update_positions(
            ps_bind, ps_guess, nv,
            ks, nk,
            Ris);
    }

    delete[] Ris;
}