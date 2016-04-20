#include "scene.h"

#include "apputil.h"
#include "renderer.h"
#include "app.h"
#include "flythrough_camera.h"
#include "arap.h"

#include "imgui.h"
#include "tiny_obj_loader.h"
#include "stb_image.h"
#include "halfedge.h"

#include <DirectXMath.h>

using namespace DirectX;

// to make HLSL compile as C++
using float4 = XMFLOAT4;
using float4x4 = XMFLOAT4X4;
using uint4 = XMUINT4;

#include "shaders/common.hlsl"

#include <unordered_map>
#include <random>
#include <algorithm>
#include <memory>

struct VertexPosition
{
    XMFLOAT3 Position;
};

struct VertexSelectorColor
{
    XMFLOAT4 Color;
};

static const XMFLOAT4 kUnselectedVertexColor = { 0.9f, 0.9f, 0.9f, 1.0f };
static const XMFLOAT4 kFixedVertexColor = { 0.8f, 0.1f, 0.1f, 1.0f };
static const XMFLOAT4 kHandleVertexColor = { 0.6f, 0.6f, 0.0f, 1.0f };

struct VertexTexCoord
{
    XMFLOAT2 TexCoord;
};

struct VertexNormal
{
    XMFLOAT3 Normal;
};

struct VertexTangent
{
    XMFLOAT4 Tangent;
};

struct VertexBitangent
{
    XMFLOAT3 Bitangent;
};

struct Texture
{
    std::string Name;
    ComPtr<ID3D11Resource> Resource;
    ComPtr<ID3D11ShaderResourceView> SRV;
};

struct Material
{
    std::string Name;
    XMFLOAT3 Ambient;
    XMFLOAT3 Diffuse;
    XMFLOAT3 Specular;
    float Shininess;
    float Opacity;
    int DiffuseTextureID;
    int SpecularTextureID;
    int BumpTextureID;
};

struct StaticMesh
{
    std::string Name;
    ComPtr<ID3D11Buffer> pPositionVertexBuffer;
    ComPtr<ID3D11Buffer> pSelectorColorVertexBuffer;
    ComPtr<ID3D11Buffer> pTexCoordVertexBuffer;
    ComPtr<ID3D11Buffer> pNormalVertexBuffer;
    ComPtr<ID3D11Buffer> pTangentVertexBuffer;
    ComPtr<ID3D11Buffer> pBitangentVertexBuffer;
    ComPtr<ID3D11Buffer> pIndexBuffer;
    int MaterialID; // the material this mesh was designed for
    UINT IndexCountPerInstance;
    UINT StartIndexLocation;

    std::shared_ptr<std::vector<XMFLOAT3>> BindPoseCPUPositions;
    std::shared_ptr<std::vector<XMFLOAT3>> CPUPositions;
    std::shared_ptr<HalfedgeMesh> Halfedge;
    std::shared_ptr<std::vector<float>> EdgeWeights;
    std::shared_ptr<std::vector<float>> ARAPSystemMatrix;
    std::shared_ptr<std::vector<int>> ControlVertexIDs;
    std::shared_ptr<std::vector<int>> VertexConstraintStatuses;
    bool SystemMatrixNeedsRebuild;
};

struct NodeTransform
{
    XMVECTOR Scaling;
    XMVECTOR RotationOrigin;
    XMVECTOR RotationQuaternion;
    XMVECTOR Translation;
};

struct StaticMeshNode
{
    int StaticMeshID;
};

enum SceneNodeType
{
    SCENENODETYPE_STATICMESH
};

struct SceneNode
{
    NodeTransform Transform;

    int MaterialID;
    
    SceneNodeType Type;

    union
    {
        StaticMeshNode AsStaticMesh;
    };
};

struct Scene
{
    std::vector<Texture> Textures;
    std::unordered_map<std::string, int> TextureNameToID;
    std::vector<Material> Materials;
    std::vector<StaticMesh> StaticMeshes;
    std::vector<SceneNode> SceneNodes;

    D3D11_VIEWPORT SceneViewport;
    XMMATRIX SceneProjection;
    XMMATRIX SceneWorldView;

    ComPtr<ID3D11Texture2D> pSceneDepthTex2D;
    ComPtr<ID3D11DepthStencilView> pSceneDepthDSV;
    ComPtr<ID3D11ShaderResourceView> pSceneDepthSRV;
    ComPtr<ID3D11SamplerState> pSceneDepthSampler;

    ComPtr<ID3D11Texture2D> pSceneNormalTex2D;
    ComPtr<ID3D11RenderTargetView> pSceneNormalRTV;
    ComPtr<ID3D11ShaderResourceView> pSceneNormalSRV;
    ComPtr<ID3D11SamplerState> pSceneNormalSampler;

    XMFLOAT3 CameraPos;
    XMFLOAT3 CameraLook;
    XMMATRIX CameraWorldViewProjection;
    ComPtr<ID3D11Buffer> pCameraBuffer;
    ComPtr<ID3D11Buffer> pViewportBuffer;
    ComPtr<ID3D11Buffer> pMaterialBuffer;
    ComPtr<ID3D11Buffer> pSceneNodeBuffer;

    ComPtr<ID3D11SamplerState> pDiffuseSampler;
    ComPtr<ID3D11SamplerState> pSpecularSampler;
    ComPtr<ID3D11SamplerState> pBumpSampler;

    ComPtr<ID3D11InputLayout> pSceneInputLayout;
    ComPtr<ID3D11RasterizerState> pSceneRasterizerState;
    ComPtr<ID3D11DepthStencilState> pSceneDepthStencilState;
    ComPtr<ID3D11BlendState> pSceneBlendState;

    Shader* SceneVS;
    Shader* SceneGS;
    Shader* ScenePS;

    ComPtr<ID3D11Texture2D> pSSAONoiseTexture;
    ComPtr<ID3D11ShaderResourceView> pSSAONoiseSRV;
    ComPtr<ID3D11SamplerState> pSSAONoiseSampler;

    ComPtr<ID3D11RasterizerState> pSSAORasterizerState;
    ComPtr<ID3D11DepthStencilState> pSSAODepthStencilState;
    ComPtr<ID3D11BlendState> pSSAOBlendState;

    Shader* SSAOVS;
    Shader* SSAOPS;

    ComPtr<ID3D11Texture2D> pSelectorTex2D;
    ComPtr<ID3D11RenderTargetView> pSelectorRTV;
    ComPtr<ID3D11UnorderedAccessView> pSelectorUAV;

    ComPtr<ID3D11InputLayout> pSelectorSphereInputLayout;
    ComPtr<ID3D11RasterizerState> pSelectorRasterizerState;
    ComPtr<ID3D11DepthStencilState> pSelectorDepthStencilState;
    ComPtr<ID3D11BlendState> pSelectorBlendState;

    ComPtr<ID3D11Buffer> pCurrSelectedVertexID;

    Shader* SelectorVS;
    Shader* SelectorGS;
    Shader* SelectorPS;

    ComPtr<ID3D11Buffer> pROISelectorBuffer;
    int ROISelectorBufferSizeInBytes;

    ComPtr<ID3D11InputLayout> pROIInputLayout;
    ComPtr<ID3D11RasterizerState> pROIRasterizerState;
    ComPtr<ID3D11DepthStencilState> pROIDepthStencilState;
    ComPtr<ID3D11BlendState> pROIBlendState;

    Shader* ROIVS;
    Shader* ROIGS;
    Shader* ROIPS;

    uint64_t LastTicks;
    int LastMouseX, LastMouseY;

    int ModelingSceneNodeID;
    int DragVertexID;

    bool ROISelectionActive;
    std::vector<XMFLOAT2> ROISelectionPoints;
};

Scene g_Scene;

static void SceneAddObjMesh(
    const char* filename, const char* mtlbasepath,
    const char* defaultName = NULL,
    std::vector<int>* newStaticMeshIDs = NULL,
    std::vector<int>* newMaterialIDs = NULL)
{
    ID3D11Device* dev = RendererGetDevice();
    ID3D11DeviceContext* dc = RendererGetDeviceContext();

    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    if (!tinyobj::LoadObj(shapes, materials, err, filename, mtlbasepath))
    {
        SimpleMessageBox_FatalError("Failed to load mesh: %s\nReason: %s", filename, err.c_str());
    }

    int firstMaterial = (int)g_Scene.Materials.size();

    for (tinyobj::material_t& material : materials)
    {
        Material m;
        m.Name = material.name;
        m.Ambient.x = material.ambient[0];
        m.Ambient.y = material.ambient[1];
        m.Ambient.z = material.ambient[2];
        m.Diffuse.x = material.diffuse[0];
        m.Diffuse.y = material.diffuse[1];
        m.Diffuse.z = material.diffuse[2];
        m.Specular.x = material.specular[0];
        m.Specular.y = material.specular[1];
        m.Specular.z = material.specular[2];
        m.Shininess = material.shininess;
        m.Opacity = material.dissolve;
        m.DiffuseTextureID = -1;
        m.SpecularTextureID = -1;
        m.BumpTextureID = -1;

        struct TextureToLoad
        {
            enum TextureToLoadType
            {
                TTLTYPE_DIFFUSE,
                TTLTYPE_SPECULAR,
                TTLTYPE_BUMP,
                TTLTYPE_Count
            };

            std::string Name;
            TextureToLoadType Type;
            int* pID;
        };

        TextureToLoad texturesToLoad[] = {
            TextureToLoad { material.diffuse_texname, TextureToLoad::TTLTYPE_DIFFUSE, &m.DiffuseTextureID },
            TextureToLoad { material.specular_texname, TextureToLoad::TTLTYPE_SPECULAR, &m.SpecularTextureID },
            TextureToLoad { material.bump_texname, TextureToLoad::TTLTYPE_BUMP, &m.BumpTextureID }
        };

        for (TextureToLoad& ttl : texturesToLoad)
        {
            if (ttl.Name.empty())
                continue;

            std::string texturePath = mtlbasepath + ttl.Name;
            auto foundTexture = g_Scene.TextureNameToID.find(texturePath);
            if (foundTexture == end(g_Scene.TextureNameToID))
            {
                static const DXGI_FORMAT kTextureTypeToFormat[TextureToLoad::TTLTYPE_Count] = {
                    DXGI_FORMAT_R8G8B8A8_TYPELESS,
                    DXGI_FORMAT_R8_TYPELESS,
                    DXGI_FORMAT_R8_TYPELESS
                };

                static const DXGI_FORMAT kTextureTypeToSRVFormat[TextureToLoad::TTLTYPE_Count] = {
                    DXGI_FORMAT_R8G8B8A8_UNORM_SRGB,
                    DXGI_FORMAT_R8_UNORM,
                    DXGI_FORMAT_R8_UNORM
                };

                static const int kTextureTypeToReqComp[TextureToLoad::TTLTYPE_Count] = {
                    4,
                    1,
                    1
                };

                int width, height, comp;
                int req_comp = kTextureTypeToReqComp[ttl.Type];
                stbi_uc* imgbytes = stbi_load(texturePath.c_str(), &width, &height, &comp, req_comp);
                if (imgbytes == NULL)
                {
                    SimpleMessageBox_FatalError("stbi_load(%s) failed.\nReason: %s", texturePath.c_str(), stbi_failure_reason());
                }

                ComPtr<ID3D11Texture2D> pTexture;

                D3D11_TEXTURE2D_DESC textureDesc = CD3D11_TEXTURE2D_DESC(kTextureTypeToFormat[ttl.Type], width, height);
                textureDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET;
                textureDesc.MiscFlags = D3D11_RESOURCE_MISC_GENERATE_MIPS;
                CHECKHR(dev->CreateTexture2D(&textureDesc, NULL, &pTexture));

                ComPtr<ID3D11ShaderResourceView> pSRV;

                CHECKHR(dev->CreateShaderResourceView(
                    pTexture.Get(),
                    &CD3D11_SHADER_RESOURCE_VIEW_DESC(D3D11_SRV_DIMENSION_TEXTURE2D, kTextureTypeToSRVFormat[ttl.Type]),
                    &pSRV));

                dc->UpdateSubresource(pTexture.Get(), 0, NULL, imgbytes, width * req_comp, width * height * req_comp);
                dc->GenerateMips(pSRV.Get());

                Texture texture;
                texture.Name = texturePath;
                texture.Resource = pTexture;
                texture.SRV = pSRV;

                int textureID = (int)g_Scene.Textures.size();
                g_Scene.Textures.push_back(std::move(texture));
                g_Scene.TextureNameToID[texturePath] = textureID;

                *ttl.pID = textureID;
            }
            else
            {
                *ttl.pID = foundTexture->second;
            }
        }

        if (newMaterialIDs)
            newMaterialIDs->push_back((int)g_Scene.Materials.size());

        g_Scene.Materials.push_back(std::move(m));
    }

    for (tinyobj::shape_t& shape : shapes)
    {
        tinyobj::mesh_t& mesh = shape.mesh;

        if (mesh.positions.size() % 3 != 0)
        {
            SimpleMessageBox_FatalError("Meshes must use 3D positions");
        }

        ComPtr<ID3D11Buffer> pPositionBuffer;
        ComPtr<ID3D11Buffer> pSelectorColorBuffer;
        ComPtr<ID3D11Buffer> pTexCoordBuffer;
        ComPtr<ID3D11Buffer> pNormalBuffer;
        ComPtr<ID3D11Buffer> pTangentBuffer;
        ComPtr<ID3D11Buffer> pBitangentBuffer;
        ComPtr<ID3D11Buffer> pIndexBuffer;

        std::shared_ptr<std::vector<XMFLOAT3>> cpuPositions = std::make_shared<std::vector<XMFLOAT3>>();
        std::shared_ptr<std::vector<XMFLOAT3>> bindPoseCPUPositions = std::make_shared<std::vector<XMFLOAT3>>();
        std::shared_ptr<HalfedgeMesh> halfedge = std::make_shared<HalfedgeMesh>();
        std::shared_ptr<std::vector<float>> edgeWeights = std::make_shared<std::vector<float>>();
        std::shared_ptr<std::vector<float>> arapSystemMatrix = std::make_shared<std::vector<float>>();
        std::shared_ptr<std::vector<int>> controlVertexIDs = std::make_shared<std::vector<int>>();
        std::shared_ptr<std::vector<int>> constraintedVertexStatuses = std::make_shared<std::vector<int>>();

        int numVertices = (int)mesh.positions.size() / 3;

        // initially nothing is constrained
        constraintedVertexStatuses->resize(numVertices, 0);

        if (!mesh.positions.empty())
        {
            D3D11_SUBRESOURCE_DATA positionVertexBufferData = {};
            positionVertexBufferData.pSysMem = mesh.positions.data();

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexPosition) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE), 
                &positionVertexBufferData, 
                &pPositionBuffer));

            cpuPositions->resize(numVertices);
            memcpy(cpuPositions->data(), mesh.positions.data(), numVertices * sizeof(XMFLOAT3));
            bindPoseCPUPositions->resize(numVertices);
            memcpy(bindPoseCPUPositions->data(), mesh.positions.data(), numVertices * sizeof(XMFLOAT3));
        }

        {
            std::vector<XMFLOAT4> selectorColorInitialData(numVertices, kUnselectedVertexColor);

            D3D11_SUBRESOURCE_DATA selectorColorVertexBufferData = {};
            selectorColorVertexBufferData.pSysMem = selectorColorInitialData.data();

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexSelectorColor) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
                &selectorColorVertexBufferData,
                &pSelectorColorBuffer));
        }

        if (!mesh.texcoords.empty())
        {
            if (mesh.texcoords.size() != numVertices * 2)
                SimpleMessageBox_FatalError("TexCoord conversion required (Expected 2D, got %dD)", mesh.texcoords.size() / numVertices);

            // flip all the texcoord.y (GL -> DX convention)
            for (int i = 1; i < mesh.texcoords.size(); i += 2)
            {
                mesh.texcoords[i] = 1.0f - mesh.texcoords[i];
            }

            D3D11_SUBRESOURCE_DATA texcoordVertexBufferData = {};
            texcoordVertexBufferData.pSysMem = mesh.texcoords.data();

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexTexCoord) * numVertices,  D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_IMMUTABLE), 
                &texcoordVertexBufferData, 
                &pTexCoordBuffer));
        }

        if (!mesh.normals.empty())
        {
            if (mesh.normals.size() != numVertices * 3)
                SimpleMessageBox_FatalError("Normal conversion required (Expected 3D, got %dD)", mesh.normals.size() / numVertices);

            D3D11_SUBRESOURCE_DATA normalVertexBufferData = {};
            normalVertexBufferData.pSysMem = mesh.normals.data();

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexNormal) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_IMMUTABLE),
                &normalVertexBufferData, 
                &pNormalBuffer));
        }

        UINT numIndices = (UINT)mesh.indices.size();

        static_assert(sizeof(mesh.indices[0]) == sizeof(UINT32), "Expecting UINT32 indices");

        D3D11_SUBRESOURCE_DATA indexBufferData = {};
        indexBufferData.pSysMem = mesh.indices.data();

        CHECKHR(dev->CreateBuffer(
            &CD3D11_BUFFER_DESC(sizeof(UINT32) * numIndices, D3D11_BIND_INDEX_BUFFER, D3D11_USAGE_IMMUTABLE), 
            &indexBufferData, 
            &pIndexBuffer));

        const int numFaces = (int)shape.mesh.indices.size() / 3;

        *halfedge = HalfedgeFromIndexedTriangles(
            mesh.positions.data(), numVertices,
            shape.mesh.indices.data(), numFaces);

        // Compute cotangent weights for every halfedge
        edgeWeights->resize(halfedge->Halfedges.size() / 2);
        for (int halfedgeID = 0; halfedgeID < (int)halfedge->Halfedges.size(); halfedgeID += 2)
        {
            HalfedgeMesh::Halfedge h1 = halfedge->Halfedges[halfedgeID];
            HalfedgeMesh::Halfedge h2 = halfedge->Halfedges[halfedgeID + 1];

            XMVECTOR v1 = XMLoadFloat3((XMFLOAT3*)&mesh.positions[3 * h2.VertexID]);
            XMVECTOR v2 = XMLoadFloat3((XMFLOAT3*)&mesh.positions[3 * h1.VertexID]);

            float weight = 0.0f;

            if (h1.FaceID != -1)
            {
                XMVECTOR v3 = XMLoadFloat3((XMFLOAT3*)&mesh.positions[3 * halfedge->Halfedges[h1.NextHalfedgeID].VertexID]);

                XMVECTOR a = v1 - v3;
                XMVECTOR b = v2 - v3;
                weight += XMVectorGetX(XMVector3Dot(a, b) / XMVector3Length(XMVector3Cross(a, b)) / XMVector3Length(a) / XMVector3Length(b));
            }

            if (h2.FaceID != -1)
            {
                XMVECTOR v3 = XMLoadFloat3((XMFLOAT3*)&mesh.positions[3 * halfedge->Halfedges[h2.NextHalfedgeID].VertexID]);

                XMVECTOR a = v1 - v3;
                XMVECTOR b = v2 - v3;
                weight += XMVectorGetX(XMVector3Dot(a, b) / XMVector3Length(XMVector3Cross(a, b)) / XMVector3Length(a) / XMVector3Length(b));
            }

            weight /= 2.0f;

            if (weight < 0.0f)
            {
                // hack to deal with negative weights, caused by cotangent weights with alpha + beta > pi
                weight = 0.0f;
            }

            (*edgeWeights)[halfedgeID / 2] = weight;
        }

        arapSystemMatrix->resize(ARAP_PACKED_SYSTEM_SIZE_IN_FLOATS(numVertices));

        // Generate tangents if possible.
        // Note: The handedness of the local coordinate system is stored as +/-1 in the w-coordinate
        // Based on:
        //  Lengyel, Eric. "Computing Tangent Space Basis Vectors for an Arbitrary Mesh".
        //  Terathon Software 3D Graphics Library, 2001. http://www.terathon.com/code/tangent.html
        if (!mesh.positions.empty() && !mesh.texcoords.empty() && !mesh.normals.empty())
        {
            XMFLOAT3* tan1 = new XMFLOAT3[numVertices * 2];
            XMFLOAT3* tan2 = tan1 + numVertices;
            ZeroMemory(tan1, numVertices * sizeof(XMFLOAT3) * 2);

            for (int face = 0; face < numFaces; face++)
            {
                int i1 = shape.mesh.indices[face * 3 + 0];
                int i2 = shape.mesh.indices[face * 3 + 1];
                int i3 = shape.mesh.indices[face * 3 + 2];

                const XMFLOAT3& v1 = (const XMFLOAT3&)shape.mesh.positions[i1 * 3];
                const XMFLOAT3& v2 = (const XMFLOAT3&)shape.mesh.positions[i2 * 3];
                const XMFLOAT3& v3 = (const XMFLOAT3&)shape.mesh.positions[i3 * 3];

                const XMFLOAT2& w1 = (const XMFLOAT2&)shape.mesh.texcoords[i1 * 2];
                const XMFLOAT2& w2 = (const XMFLOAT2&)shape.mesh.texcoords[i2 * 2];
                const XMFLOAT2& w3 = (const XMFLOAT2&)shape.mesh.texcoords[i3 * 2];

                float x1 = v2.x - v1.x;
                float x2 = v3.x - v1.x;
                float y1 = v2.y - v1.y;
                float y2 = v3.y - v1.y;
                float z1 = v2.z - v1.z;
                float z2 = v3.z - v1.z;

                float s1 = w2.x - w1.x;
                float s2 = w3.x - w1.x;
                float t1 = w2.y - w1.y;
                float t2 = w3.y - w1.y;

                float r = 1.0f / (s1 * t2 - s2 * t1);
                XMVECTOR sdir = XMVectorSet(
                    (t2 * x1 - t1 * x2) * r, 
                    (t2 * y1 - t1 * y2) * r,
                    (t2 * z1 - t1 * z2) * r,
                    0.0f);
                XMVECTOR tdir = XMVectorSet(
                    (s1 * x2 - s2 * x1) * r, 
                    (s1 * y2 - s2 * y1) * r,
                    (s1 * z2 - s2 * z1) * r,
                    0.0f);

                XMStoreFloat3(&tan1[i1], XMVectorAdd(XMLoadFloat3(&tan1[i1]), sdir));
                XMStoreFloat3(&tan1[i2], XMVectorAdd(XMLoadFloat3(&tan1[i2]), sdir));
                XMStoreFloat3(&tan1[i3], XMVectorAdd(XMLoadFloat3(&tan1[i3]), sdir));

                XMStoreFloat3(&tan2[i1], XMVectorAdd(XMLoadFloat3(&tan2[i1]), tdir));
                XMStoreFloat3(&tan2[i2], XMVectorAdd(XMLoadFloat3(&tan2[i2]), tdir));
                XMStoreFloat3(&tan2[i3], XMVectorAdd(XMLoadFloat3(&tan2[i3]), tdir));
            }

            float* tangents = new float[numVertices * 4];
            float* bitangents = new float[numVertices * 3];

            for (int vertex = 0; vertex < numVertices; vertex++)
            {
                XMVECTOR n = XMLoadFloat3((const XMFLOAT3*)&shape.mesh.normals[vertex * 3]);
                XMVECTOR t = XMLoadFloat3(&tan1[vertex]);

                // Gram-Schmidt orthogonalize
                XMVECTOR tangentxyz = XMVector3Normalize(t - n * XMVector3Dot(n, t));
                XMStoreFloat3((XMFLOAT3*)&tangents[vertex * 4], tangentxyz);

                // Calculate handedness
                tangents[vertex * 4 + 3] = (XMVectorGetX(XMVector3Dot(XMVector3Cross(n, t), XMLoadFloat3(&tan2[vertex]))) < 0.0f) ? -1.0f : 1.0f;
                
                // bitangent
                XMStoreFloat3((XMFLOAT3*)&bitangents[vertex * 3], XMVector3Cross(n, tangentxyz) * tangents[vertex * 4 + 3]);
            }

            delete[] tan1;

            // Upload to GPU
            D3D11_SUBRESOURCE_DATA tangentVertexBufferData = {};
            tangentVertexBufferData.pSysMem = tangents;

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexTangent) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_IMMUTABLE),
                &tangentVertexBufferData,
                &pTangentBuffer));

            D3D11_SUBRESOURCE_DATA bitangentVertexBufferData = {};
            bitangentVertexBufferData.pSysMem = bitangents;

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexBitangent) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_IMMUTABLE),
                &bitangentVertexBufferData,
                &pBitangentBuffer));

            delete[] tangents;
            delete[] bitangents;
        }

        int firstFace = 0;

        for (int face = 0; face < numFaces; face++)
        {
            int currMTL = shape.mesh.material_ids[face];
            
            int nextMTL = -2; // -1 already used to indicate "no material" by tinyobjloader
            if (face + 1 < numFaces)
                nextMTL = shape.mesh.material_ids[face + 1];
            
            if (currMTL == nextMTL)
            {
                // still same material, don't need to output mesh yet
                continue;
            }

            StaticMesh sm = {};
            sm.Name = shape.name;
            if (sm.Name.empty() && defaultName)
                sm.Name = defaultName;
            sm.pPositionVertexBuffer = pPositionBuffer;
            sm.pSelectorColorVertexBuffer = pSelectorColorBuffer;
            sm.pTexCoordVertexBuffer = pTexCoordBuffer;
            sm.pNormalVertexBuffer = pNormalBuffer;
            sm.pTangentVertexBuffer = pTangentBuffer;
            sm.pBitangentVertexBuffer = pBitangentBuffer;
            sm.pIndexBuffer = pIndexBuffer;
            sm.MaterialID = firstMaterial + currMTL;
            sm.IndexCountPerInstance = (face + 1 - firstFace) * 3;
            sm.StartIndexLocation = firstFace * 3;
            sm.CPUPositions = cpuPositions;
            sm.BindPoseCPUPositions = bindPoseCPUPositions;
            sm.Halfedge = halfedge;
            sm.EdgeWeights = edgeWeights;
            sm.ARAPSystemMatrix = arapSystemMatrix;
            sm.ControlVertexIDs = controlVertexIDs;
            sm.VertexConstraintStatuses = constraintedVertexStatuses;
            sm.SystemMatrixNeedsRebuild = true;

            if (newStaticMeshIDs)
                newStaticMeshIDs->push_back((int)g_Scene.StaticMeshes.size());

            g_Scene.StaticMeshes.push_back(std::move(sm));

            // first face for next mesh
            firstFace = face + 1;
        }
    }
}

static int SceneAddStaticMeshSceneNode(int staticMeshID)
{
    const StaticMesh& staticMesh = g_Scene.StaticMeshes[staticMeshID];

    SceneNode sceneNode;
    sceneNode.Transform.Scaling = XMVectorSet(1.0f, 1.0f, 1.0f, 1.0f);
    sceneNode.Transform.RotationOrigin = XMVectorSet(0.0f, 0.0f, 0.0f, 1.0f);
    sceneNode.Transform.RotationQuaternion = XMQuaternionIdentity();
    sceneNode.Transform.Translation = XMVectorSet(0.0f, 0.0f, 0.0f, 1.0f);
    sceneNode.MaterialID = staticMesh.MaterialID;
    sceneNode.Type = SCENENODETYPE_STATICMESH;
    sceneNode.AsStaticMesh.StaticMeshID = staticMeshID;

    g_Scene.SceneNodes.push_back(std::move(sceneNode));
    return (int)g_Scene.SceneNodes.size() - 1;
}

void SceneInit()
{
    ID3D11Device* dev = RendererGetDevice();

    std::vector<std::string> meshesToLoad = {
        "armadillo_1k",
        "bar1",
        "bar2",
        "bar3",
        "cactus_highres",
        "cactus_small",
        "cylinder_small",
        // doesn't have perfect indexing "dino",
        "indorelax",
        "square_21",
        "square_21_spikes"
    };

    std::vector<int> newStaticMeshIDs;
    for (const std::string& meshToLoad : meshesToLoad)
    {
        std::string meshFolder = "assets/" + meshToLoad + "/";
        std::string meshFile = meshFolder + meshToLoad + ".obj";
        SceneAddObjMesh(meshFile.c_str(), meshFolder.c_str(), meshToLoad.c_str(), &newStaticMeshIDs);
    }

    for (int newStaticMeshID : newStaticMeshIDs)
    {
        int sceneNodeID = SceneAddStaticMeshSceneNode(newStaticMeshID);

        SceneNode& sceneNode = g_Scene.SceneNodes[sceneNodeID];

        std::string name = g_Scene.StaticMeshes[newStaticMeshID].Name;

        if (name == "armadillo_1k")
        {
            sceneNode.Transform.RotationQuaternion = XMQuaternionRotationRollPitchYaw(0.0f, XMConvertToRadians(180.0f), 0.0f);
        }
        else if (name == "bar1" || name == "bar2")
        {
            sceneNode.Transform.Scaling = XMVectorSet(1.0f / 150.0f, 1.0f / 150.0f, 1.0f / 150.0f, 1.0f);
            sceneNode.Transform.Translation = XMVectorSet(-0.1f, -0.4f, 0.3f, 0.0f);
        }
        else if (name == "cactus_highres" || name == "cactus_small")
        {
            sceneNode.Transform.Translation = XMVectorSet(0.5f, -0.5f, 0.7f, 0.0f);
            sceneNode.Transform.RotationQuaternion = XMQuaternionRotationRollPitchYaw(XMConvertToRadians(-90.0f), XMConvertToRadians(90.0f), 0.0f);
        }
        else if (name == "square_21")
        {
            sceneNode.Transform.Translation = XMVectorSet(-0.5f, -0.5f, 0.0f, 0.0f);
        }
        else if (name == "square_21_spikes")
        {
            sceneNode.Transform.Translation = XMVectorSet(-40.5f, -13.4f, -0.5f, 0.0f);
            sceneNode.Transform.RotationQuaternion = XMQuaternionRotationRollPitchYaw(XMConvertToRadians(-90.0f), 0.0f, 0.0f);
            sceneNode.Transform.RotationOrigin = -sceneNode.Transform.Translation;
        }
    }

    g_Scene.ModelingSceneNodeID = 0;

    g_Scene.SceneVS = RendererAddShader("scene.hlsl", "VSmain", "vs_5_0");
    g_Scene.SceneGS = RendererAddShader("scene.hlsl", "GSmain", "gs_5_0");
    g_Scene.ScenePS = RendererAddShader("scene.hlsl", "PSmain", "ps_5_0");
    g_Scene.SSAOVS = RendererAddShader("ssao.hlsl", "VSmain", "vs_5_0");
    g_Scene.SSAOPS = RendererAddShader("ssao.hlsl", "PSmain", "ps_5_0");
    g_Scene.SelectorVS = RendererAddShader("selector.hlsl", "VSmain", "vs_5_0");
    g_Scene.SelectorGS = RendererAddShader("selector.hlsl", "GSmain", "gs_5_0");
    g_Scene.SelectorPS = RendererAddShader("selector.hlsl", "PSmain", "ps_5_0");
    g_Scene.ROIVS = RendererAddShader("roi.hlsl", "VSmain", "vs_5_0");
    g_Scene.ROIGS = RendererAddShader("roi.hlsl", "GSmain", "gs_5_0");
    g_Scene.ROIPS = RendererAddShader("roi.hlsl", "PSmain", "ps_5_0");

    // XMStoreFloat3(&g_Scene.CameraPos, XMVectorSet(-100.0f, 100.0f, -100.0f, 1.0f));
    XMStoreFloat3(&g_Scene.CameraPos, XMVectorSet(0.0f, 0.4f, 0.7f, 1.0f));
    XMStoreFloat3(&g_Scene.CameraLook, XMVector3Normalize(XMVectorSet(0.0f, 0.0f, 0.0f, 1.0f) - XMLoadFloat3(&g_Scene.CameraPos)));

    CD3D11_SAMPLER_DESC sceneNormalSamplerDesc(D3D11_DEFAULT);
    sceneNormalSamplerDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sceneNormalSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
    sceneNormalSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
    CHECKHR(dev->CreateSamplerState(&sceneNormalSamplerDesc, &g_Scene.pSceneNormalSampler));

    CD3D11_SAMPLER_DESC sceneDepthSamplerDesc(D3D11_DEFAULT);
    sceneDepthSamplerDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sceneDepthSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
    sceneDepthSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
    CHECKHR(dev->CreateSamplerState(&sceneDepthSamplerDesc, &g_Scene.pSceneDepthSampler));

    CHECKHR(dev->CreateBuffer(
        &CD3D11_BUFFER_DESC(sizeof(PerCameraData), D3D11_BIND_CONSTANT_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
        NULL,
        &g_Scene.pCameraBuffer));

    CHECKHR(dev->CreateBuffer(
        &CD3D11_BUFFER_DESC(sizeof(PerViewportData), D3D11_BIND_CONSTANT_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
        NULL,
        &g_Scene.pViewportBuffer));

    CHECKHR(dev->CreateBuffer(
        &CD3D11_BUFFER_DESC(sizeof(PerMaterialData), D3D11_BIND_CONSTANT_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
        NULL, 
        &g_Scene.pMaterialBuffer));

    CHECKHR(dev->CreateBuffer(
        &CD3D11_BUFFER_DESC(sizeof(PerSceneNodeData), D3D11_BIND_CONSTANT_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE), 
        NULL, 
        &g_Scene.pSceneNodeBuffer));

    CD3D11_SAMPLER_DESC diffuseSamplerDesc(D3D11_DEFAULT);
    diffuseSamplerDesc.Filter = D3D11_FILTER_ANISOTROPIC;
    diffuseSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
    diffuseSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
    diffuseSamplerDesc.MaxAnisotropy = 8;
    CHECKHR(dev->CreateSamplerState(&diffuseSamplerDesc, &g_Scene.pDiffuseSampler));

    CD3D11_SAMPLER_DESC specularSamplerDesc(D3D11_DEFAULT);
    specularSamplerDesc.Filter = D3D11_FILTER_ANISOTROPIC;
    specularSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
    specularSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
    specularSamplerDesc.MaxAnisotropy = 8;
    CHECKHR(dev->CreateSamplerState(&specularSamplerDesc, &g_Scene.pSpecularSampler));

    CD3D11_SAMPLER_DESC bumpSamplerDesc(D3D11_DEFAULT);
    bumpSamplerDesc.Filter = D3D11_FILTER_ANISOTROPIC;
    bumpSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
    bumpSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
    bumpSamplerDesc.MaxAnisotropy = 8;
    CHECKHR(dev->CreateSamplerState(&bumpSamplerDesc, &g_Scene.pBumpSampler));

    D3D11_INPUT_ELEMENT_DESC sceneInputElementDescs[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 1, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 2, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "TANGENT", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 3, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "BITANGENT", 0, DXGI_FORMAT_R32G32B32_FLOAT, 4, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
    };
    CHECKHR(dev->CreateInputLayout(
        sceneInputElementDescs, _countof(sceneInputElementDescs),
        g_Scene.SceneVS->Blob->GetBufferPointer(), g_Scene.SceneVS->Blob->GetBufferSize(),
        &g_Scene.pSceneInputLayout));

    D3D11_RASTERIZER_DESC sceneRasterizerDesc = CD3D11_RASTERIZER_DESC(D3D11_DEFAULT);
    sceneRasterizerDesc.FrontCounterClockwise = TRUE;
    sceneRasterizerDesc.CullMode = D3D11_CULL_NONE;
    CHECKHR(dev->CreateRasterizerState(&sceneRasterizerDesc, &g_Scene.pSceneRasterizerState));

    D3D11_DEPTH_STENCIL_DESC sceneDepthStencilDesc = CD3D11_DEPTH_STENCIL_DESC(D3D11_DEFAULT);
    CHECKHR(dev->CreateDepthStencilState(&sceneDepthStencilDesc, &g_Scene.pSceneDepthStencilState));

    D3D11_BLEND_DESC sceneBlendDesc = CD3D11_BLEND_DESC(D3D11_DEFAULT);
    CHECKHR(dev->CreateBlendState(&sceneBlendDesc, &g_Scene.pSceneBlendState));

    D3D11_RASTERIZER_DESC ssaoRasterizerDesc = CD3D11_RASTERIZER_DESC(D3D11_DEFAULT);
    CHECKHR(dev->CreateRasterizerState(&ssaoRasterizerDesc, &g_Scene.pSSAORasterizerState));

    D3D11_DEPTH_STENCIL_DESC ssaoDepthStencilDesc = CD3D11_DEPTH_STENCIL_DESC(D3D11_DEFAULT);
    ssaoDepthStencilDesc.DepthEnable = FALSE;
    CHECKHR(dev->CreateDepthStencilState(&ssaoDepthStencilDesc, &g_Scene.pSSAODepthStencilState));

    D3D11_BLEND_DESC ssaoBlendDesc = CD3D11_BLEND_DESC(D3D11_DEFAULT);
    ssaoBlendDesc.RenderTarget[0].BlendEnable = TRUE;
    ssaoBlendDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
    ssaoBlendDesc.RenderTarget[0].DestBlend = D3D11_BLEND_ONE;
    ssaoBlendDesc.RenderTarget[0].BlendOp = D3D11_BLEND_OP_REV_SUBTRACT;
    CHECKHR(dev->CreateBlendState(&ssaoBlendDesc, &g_Scene.pSSAOBlendState));

    // Noise texture (for SSAO)
    {
        int noiseSize = 64;
        std::vector<XMFLOAT4> noiseData(noiseSize * noiseSize);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<float> d(0.0f, 1.0f);
        for (XMFLOAT4& noise : noiseData)
        {
            noise.x = d(gen);
            noise.y = d(gen);
            noise.z = d(gen);
            noise.w = d(gen);
        }
        
        D3D11_SUBRESOURCE_DATA noiseInitialData = {};
        noiseInitialData.pSysMem = noiseData.data();
        noiseInitialData.SysMemPitch = noiseSize * sizeof(XMFLOAT4);
        noiseInitialData.SysMemSlicePitch = noiseSize * noiseSize * sizeof(XMFLOAT4);
        
        CHECKHR(dev->CreateTexture2D(
            &CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32G32B32A32_FLOAT, noiseSize, noiseSize, 1, 1, D3D11_BIND_SHADER_RESOURCE),
            &noiseInitialData,
            &g_Scene.pSSAONoiseTexture));

        CHECKHR(dev->CreateShaderResourceView(
            g_Scene.pSSAONoiseTexture.Get(),
            &CD3D11_SHADER_RESOURCE_VIEW_DESC(D3D11_SRV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 1),
            &g_Scene.pSSAONoiseSRV));

        CD3D11_SAMPLER_DESC noiseSamplerDesc(D3D11_DEFAULT);
        noiseSamplerDesc.Filter = D3D11_FILTER_MIN_MAG_LINEAR_MIP_POINT;
        noiseSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
        noiseSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
        CHECKHR(dev->CreateSamplerState(&noiseSamplerDesc, &g_Scene.pSSAONoiseSampler));
    }

    // Selector
    {
        D3D11_INPUT_ELEMENT_DESC selectorInputElementDescs[] = {
            { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
            { "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
        };

        CHECKHR(dev->CreateInputLayout(
            selectorInputElementDescs, _countof(selectorInputElementDescs),
            g_Scene.SelectorVS->Blob->GetBufferPointer(), g_Scene.SelectorVS->Blob->GetBufferSize(),
            &g_Scene.pSelectorSphereInputLayout));

        D3D11_RASTERIZER_DESC selectorRasterizerDesc = CD3D11_RASTERIZER_DESC(D3D11_DEFAULT);
        CHECKHR(dev->CreateRasterizerState(&selectorRasterizerDesc, &g_Scene.pSelectorRasterizerState));

        D3D11_DEPTH_STENCIL_DESC selectorDepthStencilDesc = CD3D11_DEPTH_STENCIL_DESC(D3D11_DEFAULT);
        CHECKHR(dev->CreateDepthStencilState(&selectorDepthStencilDesc, &g_Scene.pSelectorDepthStencilState));

        D3D11_BLEND_DESC selectorBlendDesc = CD3D11_BLEND_DESC(D3D11_DEFAULT);
        selectorBlendDesc.IndependentBlendEnable = TRUE;
        selectorBlendDesc.RenderTarget[0].BlendEnable = TRUE;
        selectorBlendDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
        selectorBlendDesc.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
        CHECKHR(dev->CreateBlendState(&selectorBlendDesc, &g_Scene.pSelectorBlendState));

        CurrSelectionData selection = {};
        selection.VertexID.x = selection.VertexID.y = selection.VertexID.z = selection.VertexID.w = UINT_MAX;

        D3D11_SUBRESOURCE_DATA initialSelect = {};
        initialSelect.pSysMem = &selection;
        initialSelect.SysMemPitch = initialSelect.SysMemSlicePitch = sizeof(CurrSelectionData);
        CHECKHR(dev->CreateBuffer(
            &CD3D11_BUFFER_DESC(sizeof(CurrSelectionData), D3D11_BIND_CONSTANT_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE), 
            &initialSelect, 
            &g_Scene.pCurrSelectedVertexID));
    }

    // ROI
    {
        D3D11_INPUT_ELEMENT_DESC roiInputElementDescs[] = {
            { "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
        };

        CHECKHR(dev->CreateInputLayout(
            roiInputElementDescs, _countof(roiInputElementDescs),
            g_Scene.ROIVS->Blob->GetBufferPointer(), g_Scene.ROIVS->Blob->GetBufferSize(),
            &g_Scene.pROIInputLayout));

        D3D11_RASTERIZER_DESC roiRasterizerDesc = CD3D11_RASTERIZER_DESC(D3D11_DEFAULT);
        CHECKHR(dev->CreateRasterizerState(&roiRasterizerDesc, &g_Scene.pROIRasterizerState));

        D3D11_DEPTH_STENCIL_DESC roiDepthStencilDesc = CD3D11_DEPTH_STENCIL_DESC(D3D11_DEFAULT);
        roiDepthStencilDesc.DepthEnable = FALSE;
        CHECKHR(dev->CreateDepthStencilState(&roiDepthStencilDesc, &g_Scene.pROIDepthStencilState));

        D3D11_BLEND_DESC roiBlendDesc = CD3D11_BLEND_DESC(D3D11_DEFAULT);
        roiBlendDesc.RenderTarget[0].BlendEnable = TRUE;
        roiBlendDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
        roiBlendDesc.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
        CHECKHR(dev->CreateBlendState(&roiBlendDesc, &g_Scene.pROIBlendState));

        // decent initial size
        g_Scene.ROISelectorBufferSizeInBytes = sizeof(XMFLOAT2) * 10;
        CHECKHR(dev->CreateBuffer(
            &CD3D11_BUFFER_DESC(g_Scene.ROISelectorBufferSizeInBytes, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
            NULL,
            &g_Scene.pROISelectorBuffer));
    }

    g_Scene.LastMouseX = INT_MIN;
    g_Scene.LastMouseY = INT_MIN;

    g_Scene.DragVertexID = -1;
}

void SceneResize(
    int windowWidth, int windowHeight,
    int renderWidth, int renderHeight)
{
    ID3D11Device* dev = RendererGetDevice();

    CHECKHR(dev->CreateTexture2D(
        &CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32_TYPELESS, renderWidth, renderHeight, 1, 1, D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE), 
        NULL, 
        &g_Scene.pSceneDepthTex2D));

    CHECKHR(dev->CreateDepthStencilView(
        g_Scene.pSceneDepthTex2D.Get(), 
        &CD3D11_DEPTH_STENCIL_VIEW_DESC(D3D11_DSV_DIMENSION_TEXTURE2D, DXGI_FORMAT_D32_FLOAT), 
        &g_Scene.pSceneDepthDSV));
    
    CHECKHR(dev->CreateShaderResourceView(
        g_Scene.pSceneDepthTex2D.Get(),
        &CD3D11_SHADER_RESOURCE_VIEW_DESC(D3D11_SRV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32_FLOAT, 0, 1),
        &g_Scene.pSceneDepthSRV));

    CHECKHR(dev->CreateTexture2D(
        &CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32G32B32A32_FLOAT, renderWidth, renderHeight, 1, 1, D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE),
        NULL,
        &g_Scene.pSceneNormalTex2D));

    CHECKHR(dev->CreateRenderTargetView(
        g_Scene.pSceneNormalTex2D.Get(),
        &CD3D11_RENDER_TARGET_VIEW_DESC(D3D11_RTV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 0, 1),
        &g_Scene.pSceneNormalRTV));

    CHECKHR(dev->CreateShaderResourceView(
        g_Scene.pSceneNormalTex2D.Get(),
        &CD3D11_SHADER_RESOURCE_VIEW_DESC(D3D11_SRV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 1),
        &g_Scene.pSceneNormalSRV));

    CHECKHR(dev->CreateTexture2D(
        &CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32_UINT, renderWidth, renderHeight, 1, 1, D3D11_BIND_RENDER_TARGET | D3D11_BIND_UNORDERED_ACCESS, D3D11_USAGE_DEFAULT, 0),
        NULL,
        &g_Scene.pSelectorTex2D));

    CHECKHR(dev->CreateRenderTargetView(
        g_Scene.pSelectorTex2D.Get(),
        &CD3D11_RENDER_TARGET_VIEW_DESC(D3D11_RTV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32_UINT, 0, 0, 1),
        &g_Scene.pSelectorRTV));

    CHECKHR(dev->CreateUnorderedAccessView(
        g_Scene.pSelectorTex2D.Get(),
        &CD3D11_UNORDERED_ACCESS_VIEW_DESC(D3D11_UAV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32_UINT, 0, 0),
        &g_Scene.pSelectorUAV));

    g_Scene.SceneViewport = CD3D11_VIEWPORT(0.0f, 0.0f, (FLOAT)renderWidth, (FLOAT)renderHeight);

    float aspectWbyH = g_Scene.SceneViewport.Width / g_Scene.SceneViewport.Height;
    // g_Scene.SceneProjection = XMMatrixPerspectiveFovLH(XMConvertToRadians(90.0f), aspectWbyH, 1.0f, 5000.0f);
    g_Scene.SceneProjection = XMMatrixPerspectiveFovLH(XMConvertToRadians(90.0f), aspectWbyH, 0.001f, 10.0f);
}

static UINT32 PickSelector(int x, int y)
{
    if (x < g_Scene.SceneViewport.TopLeftX || x >= g_Scene.SceneViewport.TopLeftX + g_Scene.SceneViewport.Width ||
        y < g_Scene.SceneViewport.TopLeftY || y >= g_Scene.SceneViewport.TopLeftY + g_Scene.SceneViewport.Height)
    {
        return UINT_MAX;
    }

    ID3D11Device* dev = RendererGetDevice();
    ID3D11DeviceContext* dc = RendererGetDeviceContext();

    ComPtr<ID3D11Texture2D> pickTex2D;
    CHECKHR(dev->CreateTexture2D(&CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32_UINT, 1, 1, 1, 1, 0, D3D11_USAGE_STAGING, D3D11_CPU_ACCESS_READ), NULL, &pickTex2D));

    D3D11_BOX srcBox = {};
    dc->CopySubresourceRegion(pickTex2D.Get(), 0, 0, 0, 0, g_Scene.pSelectorTex2D.Get(), 0, &CD3D11_BOX(x, y, 0, x + 1, y + 1, 1));

    D3D11_MAPPED_SUBRESOURCE mapped;
    CHECKHR(dc->Map(pickTex2D.Get(), 0, D3D11_MAP_READ, 0, &mapped));

    UINT32 pickVertexID = *(UINT32*)mapped.pData;

    dc->Unmap(pickTex2D.Get(), 0);

    return pickVertexID;
}

bool SceneHandleEvent(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if (msg == WM_KEYDOWN && wParam == VK_SPACE)
    {
        if (g_Scene.ROISelectionActive)
        {
            // mark all the points inside the selection as ROI, outside as non-ROI
            if (g_Scene.ROISelectionPoints.size() >= 3 && g_Scene.ModelingSceneNodeID != -1)
            {
                SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
                StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

                // tie the loop
                {
                    XMFLOAT2 firstPoint = g_Scene.ROISelectionPoints.front();
                    g_Scene.ROISelectionPoints.push_back(firstPoint);
                }

                // Convert points to clip space
                std::vector<XMVECTOR> roiClipSpacePoints(g_Scene.ROISelectionPoints.size());
                for (int pointID = 0; pointID < (int)g_Scene.ROISelectionPoints.size(); pointID++)
                {
                    // grab normalized coordinates (texcoord coordinates)
                    XMVECTOR p = XMLoadFloat2(&g_Scene.ROISelectionPoints[pointID]);

                    // convert texcoord coordinates to clip space
                    p = XMVectorSetY(p, 1.0f - XMVectorGetY(p)) * XMVectorSet(2.0f, 2.0f, 1.0f, 1.0f) - XMVectorSet(1.0f, 1.0f, 0.0f, 0.0f);
                    
                    roiClipSpacePoints[pointID] = p;
                }

                for (int vertexID = 0; vertexID < (int)staticMesh.CPUPositions->size(); vertexID++)
                {
                    // convert vertex to clip space
                    XMMATRIX transform = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, sceneNode.Transform.Translation);
                    XMVECTOR modelPos = XMVectorSetW(XMLoadFloat3(&(*staticMesh.CPUPositions)[vertexID]), 1.0f);
                    XMVECTOR worldPos = XMVector4Transform(modelPos, transform);
                    XMVECTOR clipPos = XMVector4Transform(worldPos, g_Scene.CameraWorldViewProjection);

                    float cx = XMVectorGetX(clipPos);
                    float cy = XMVectorGetY(clipPos);
                    float cz = XMVectorGetZ(clipPos);
                    float cw = XMVectorGetW(clipPos);

                    // clipping for stuff behind the camera
                    if (cz < 0.0f)
                    {
                        // clipped, so constrained
                        (*staticMesh.VertexConstraintStatuses)[vertexID] = 1;
                        continue;
                    }

                    // convert to screen coordinates
                    float v_sx = g_Scene.SceneViewport.TopLeftX + (cx / cw + 1.0f) / 2.0f * g_Scene.SceneViewport.Width;
                    float v_sy = g_Scene.SceneViewport.TopLeftY + g_Scene.SceneViewport.Height - (cy / cw + 1.0f) / 2.0f * g_Scene.SceneViewport.Height;

                    int numIntersections = 0;

                    for (int lineID = 0; lineID < (int)g_Scene.ROISelectionPoints.size() - 1; lineID++)
                    {
                        // grab clip space line coordinates
                        XMVECTOR p0 = roiClipSpacePoints[lineID];
                        XMVECTOR p1 = roiClipSpacePoints[lineID + 1];

                        // convert to screen coordinates
                        float v_p0x = g_Scene.SceneViewport.TopLeftX + (XMVectorGetX(p0) + 1.0f) / 2.0f * g_Scene.SceneViewport.Width;
                        float v_p0y = g_Scene.SceneViewport.TopLeftY + g_Scene.SceneViewport.Height - (XMVectorGetY(p0) + 1.0f) / 2.0f * g_Scene.SceneViewport.Height;
                        float v_p1x = g_Scene.SceneViewport.TopLeftX + (XMVectorGetX(p1) + 1.0f) / 2.0f * g_Scene.SceneViewport.Width;
                        float v_p1y = g_Scene.SceneViewport.TopLeftY + g_Scene.SceneViewport.Height - (XMVectorGetY(p1) + 1.0f) / 2.0f * g_Scene.SceneViewport.Height;

                        // make p0 always the higher line (for convenience)
                        if (v_p0y > v_p1y)
                        {
                            std::swap(v_p0x, v_p1x);
                            std::swap(v_p0y, v_p1y);
                        }

                        // check if the y-range of the line overlaps the vertex
                        if (v_p0y <= v_sy && v_sy <= v_p1y)
                        {
                            // find intersection of line and vertex through similar triangles
                            float y_vdist = v_sy - v_p0y;
                            float y_len = v_p1y - v_p0y;
                            float x_len = v_p1x - v_p0x;
                            float x_vdist = y_vdist / y_len * x_len;

                            // here's the intersection
                            float isect_x = v_p0x + x_vdist;
                            float isect_y = v_sy;

                            if (isect_x < v_sx)
                            {
                                numIntersections++;
                            }
                        }
                    }

                    // odd number of lines to the left of the vertex means the vertex is inside
                    if (numIntersections % 2 == 1)
                    {
                        // inside: not constrained
                        (*staticMesh.VertexConstraintStatuses)[vertexID] = 0;
                    }
                    else
                    {
                        // outside: constrained
                        (*staticMesh.VertexConstraintStatuses)[vertexID] = 1;
                    }
                }

                // set colors based on constraint status
                ID3D11Device* dev = RendererGetDevice();
                ID3D11DeviceContext* dc = RendererGetDeviceContext();

                D3D11_MAPPED_SUBRESOURCE mappedColors;
                CHECKHR(dc->Map(staticMesh.pSelectorColorVertexBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedColors));

                VertexSelectorColor* colors = (VertexSelectorColor*)mappedColors.pData;

                for (int vertexID = 0; vertexID < (int)staticMesh.CPUPositions->size(); vertexID++)
                {
                    if ((*staticMesh.VertexConstraintStatuses)[vertexID])
                        colors[vertexID] = VertexSelectorColor{ kFixedVertexColor };
                    else
                        colors[vertexID] = VertexSelectorColor{ kUnselectedVertexColor };
                }

                dc->Unmap(staticMesh.pSelectorColorVertexBuffer.Get(), 0);

                staticMesh.ControlVertexIDs->clear();
                staticMesh.SystemMatrixNeedsRebuild = true;
            }

            g_Scene.ROISelectionActive = false;
            g_Scene.ROISelectionPoints.clear();
        }
    }

    if (msg == WM_LBUTTONDOWN)
    {
        if (!AppIsKeyPressed(VK_RBUTTON)) // don't pick when in camera mode
        {
            if (g_Scene.ROISelectionActive)
            {
                // Adding a point to the ROI selection
                RECT clientRect;
                CHECKWIN32(GetClientRect(hWnd, &clientRect));

                XMFLOAT2 normalizedPoint = {
                    (float)GET_X_LPARAM(lParam) / (clientRect.right - clientRect.left),
                    (float)GET_Y_LPARAM(lParam) / (clientRect.bottom - clientRect.top)
                };

                g_Scene.ROISelectionPoints.push_back(normalizedPoint);
            }
            else
            {
                // Vertex picking mode
                UINT32 pickVertexID = PickSelector(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
                printf("pickVertexID: %u\n", pickVertexID);

                if (pickVertexID == UINT_MAX)
                {
                    g_Scene.DragVertexID = -1;
                }
                else
                {
                    if (g_Scene.ModelingSceneNodeID != -1)
                    {
                        SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
                        StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

                        auto found = std::find(begin(*staticMesh.ControlVertexIDs), end(*staticMesh.ControlVertexIDs), (int)pickVertexID);

                        if (AppIsKeyPressed('C')) // toggle control vertices
                        {
                            ID3D11Device* dev = RendererGetDevice();
                            ID3D11DeviceContext* dc = RendererGetDeviceContext();

                            XMFLOAT4 newColor;

                            // toggle control status
                            if (found == end(*staticMesh.ControlVertexIDs))
                            {
                                printf("Adding %d to control vertices\n", (int)pickVertexID);
                                staticMesh.ControlVertexIDs->push_back((int)pickVertexID);
                                (*staticMesh.VertexConstraintStatuses)[pickVertexID] = 1;
                                newColor = kHandleVertexColor;
                            }
                            else
                            {
                                printf("Removing %d from control vertices\n", (int)pickVertexID);
                                staticMesh.ControlVertexIDs->erase(found);
                                (*staticMesh.VertexConstraintStatuses)[pickVertexID] = 0;
                                newColor = kUnselectedVertexColor;
                            }

                            D3D11_SUBRESOURCE_DATA newColorData = {};
                            newColorData.pSysMem = &newColor;

                            ComPtr<ID3D11Buffer> pNewColorBuffer;
                            CHECKHR(dev->CreateBuffer(
                                &CD3D11_BUFFER_DESC(sizeof(VertexSelectorColor), 0, D3D11_USAGE_STAGING, D3D11_CPU_ACCESS_WRITE),
                                &newColorData,
                                &pNewColorBuffer));

                            dc->CopySubresourceRegion(staticMesh.pSelectorColorVertexBuffer.Get(), 0, sizeof(VertexSelectorColor) * pickVertexID, 0, 0, pNewColorBuffer.Get(), 0, NULL);

                            staticMesh.SystemMatrixNeedsRebuild = true;
                        }
                        else
                        {
                            // new vertex to drag for modeling
                            if (found != end(*staticMesh.ControlVertexIDs))
                                g_Scene.DragVertexID = (int)pickVertexID;
                        }
                    }
                }
            }
        }
    }
    
    return false;
}

static void SceneShowToolboxGUI()
{
    ImGuiIO& io = ImGui::GetIO();
    int w = int(io.DisplaySize.x / io.DisplayFramebufferScale.x);
    int h = int(io.DisplaySize.y / io.DisplayFramebufferScale.y);

    int toolboxW = 300, toolboxH = 300;

    ImGui::SetNextWindowSize(ImVec2((float)toolboxW, (float)toolboxH), ImGuiSetCond_Always);
    ImGui::SetNextWindowPos(ImVec2((float)w - toolboxW, 0), ImGuiSetCond_Always);
    if (ImGui::Begin("Toolbox", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize))
    {
        if (g_Scene.ModelingSceneNodeID != -1)
        {
            SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
            StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

            if (ImGui::Button("Select Region of Interest"))
            {
                g_Scene.ROISelectionActive = true;
                g_Scene.ROISelectionPoints.clear();
            }

            if (ImGui::Button("Refactorize ARAP"))
            {
                arap_factorize_system(
                    (int)staticMesh.CPUPositions->size(),
                    staticMesh.VertexConstraintStatuses->data(),
                    (const int*)staticMesh.Halfedge->Vertices.data(),
                    (const int*)staticMesh.Halfedge->Halfedges.data(),
                    staticMesh.EdgeWeights->data(),
                    staticMesh.ARAPSystemMatrix->data());

                staticMesh.SystemMatrixNeedsRebuild = false;
            }

            if (ImGui::Button("Reset Constraints"))
            {
                staticMesh.ControlVertexIDs->clear();
                staticMesh.VertexConstraintStatuses->clear();
                staticMesh.VertexConstraintStatuses->resize(staticMesh.CPUPositions->size(), 0);
                // Reset colors
                {
                    ID3D11DeviceContext* dc = RendererGetDeviceContext();

                    D3D11_MAPPED_SUBRESOURCE mappedColors;
                    CHECKHR(dc->Map(staticMesh.pSelectorColorVertexBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedColors));

                    VertexSelectorColor* colors = (VertexSelectorColor*)mappedColors.pData;

                    for (int vertexID = 0; vertexID < (int)staticMesh.CPUPositions->size(); vertexID++)
                    {
                        colors[vertexID] = VertexSelectorColor{ kUnselectedVertexColor };
                    }

                    dc->Unmap(staticMesh.pSelectorColorVertexBuffer.Get(), 0);
                }
                staticMesh.SystemMatrixNeedsRebuild = true;
            }

            if (ImGui::Button("Reset Bind Pose"))
            {
                ID3D11DeviceContext* dc = RendererGetDeviceContext();

                D3D11_MAPPED_SUBRESOURCE mapped;
                CHECKHR(dc->Map(staticMesh.pPositionVertexBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped));

                memcpy(staticMesh.CPUPositions->data(), staticMesh.BindPoseCPUPositions->data(), staticMesh.CPUPositions->size() * sizeof(XMFLOAT3));
                memcpy(mapped.pData, staticMesh.CPUPositions->data(), staticMesh.CPUPositions->size() * sizeof(XMFLOAT3));

                dc->Unmap(staticMesh.pPositionVertexBuffer.Get(), 0);
            }

            std::vector<const char*> modelNames;
            for (int i = 0; i < (int)g_Scene.SceneNodes.size(); i++)
            {
                assert(g_Scene.SceneNodes[i].Type == SCENENODETYPE_STATICMESH);
                const StaticMesh& sm = g_Scene.StaticMeshes[g_Scene.SceneNodes[i].AsStaticMesh.StaticMeshID];
                modelNames.push_back(sm.Name.c_str());
            }
            ImGui::Text("Model selection");
            if (ImGui::ListBox("##modelselect", &g_Scene.ModelingSceneNodeID, modelNames.data(), (int)modelNames.size()))
            {
                g_Scene.DragVertexID = -1;
            }
        }
    }
    ImGui::End();

    int tutorialW = 600, tutorialH = 190;

    ImGui::SetNextWindowSize(ImVec2((float)tutorialW, (float)tutorialH), ImGuiSetCond_Always);
    ImGui::SetNextWindowPos(ImVec2(0, (float)h - tutorialH), ImGuiSetCond_Always);
    if (ImGui::Begin("Tutorial", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize))
    {
        if (g_Scene.ModelingSceneNodeID == -1)
        {
            ImGui::Text("A model needs to be loaded");
        }
        else
        {
            SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
            StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

            bool anyConstrained = std::any_of(begin(*staticMesh.VertexConstraintStatuses), end(*staticMesh.VertexConstraintStatuses), [](int s) { return s; });

            if (g_Scene.ROISelectionActive)
            {
                ImGui::Text("Left click to create polygon selection vertices.");
                ImGui::Text("Press spacebar when done.");
            }
            else if (anyConstrained && !staticMesh.ControlVertexIDs->empty() && staticMesh.SystemMatrixNeedsRebuild)
            {
                ImGui::Text("Constraints have changed, the ARAP system matrix needs to be rebuilt.");
                ImGui::Text("Click \"Refactorize ARAP\" in the Toolbox before dragging handles.");
            }
            else if (staticMesh.ControlVertexIDs->empty())
            {
                if (anyConstrained)
                {
                    ImGui::Text("Select some handle points by left-clicking while holding the \"C\" key.");
                }
                else
                {
                    ImGui::Text("Click \"Select Region of Interest\" to define which points to deform");
                }
            }
            else
            {
                ImGui::Text("Left click drag handle points to deform the mesh");
            }

            ImGui::Text("-----");
            ImGui::Text("Camera controls:");
            ImGui::Text("Hold right click - Activate camera");
            if (AppIsKeyPressed(VK_RBUTTON))
            {
                ImGui::Text("Mouse movement - Look around");
                ImGui::Text("WASD - Forward/Left/Backward/Right");
                ImGui::Text("Spacebar - Go up");
                ImGui::Text("Left Ctrl - Go down");
            }
        }
    }
    ImGui::End();
}

void ScenePaint(ID3D11RenderTargetView* pBackBufferRTV)
{
    SceneShowToolboxGUI();

    uint64_t currTicks;
    QueryPerformanceCounter((LARGE_INTEGER*)&currTicks);
    
    if (g_Scene.LastTicks == 0) {
        g_Scene.LastTicks = currTicks;
    }
    
    uint64_t deltaTicks = currTicks - g_Scene.LastTicks;

    uint64_t ticksPerSecond;
    QueryPerformanceFrequency((LARGE_INTEGER*)&ticksPerSecond);

    int currMouseX, currMouseY;
    AppGetClientCursorPos(&currMouseX, &currMouseY);
    
    // Initialize the last mouse position on the first update
    if (g_Scene.LastMouseX == INT_MIN)
        g_Scene.LastMouseX = currMouseX;
    if (g_Scene.LastMouseY == INT_MIN)
        g_Scene.LastMouseY = currMouseY;

    int deltaMouseX = currMouseX - g_Scene.LastMouseX;
    int deltaMouseY = currMouseY - g_Scene.LastMouseY;

    if (!AppIsKeyPressed(VK_LBUTTON))
    {
        g_Scene.DragVertexID = -1;
    }

    ID3D11Device* dev = RendererGetDevice();
    ID3D11DeviceContext* dc = RendererGetDeviceContext();

    // Update camera
    XMFLOAT3 cameraUp;
    XMStoreFloat3(&cameraUp, XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f));
    static const float kCameraUp[3] = { 0.0f, 1.0f, 0.0f };
    {
        float activated = AppIsKeyPressed(VK_RBUTTON) ? 1.0f : 0.0f;
        XMFLOAT4X4 worldView;
        if (activated)
            flythrough_camera_update(
                &g_Scene.CameraPos.x,
                &g_Scene.CameraLook.x,
                &cameraUp.x,
                &worldView.m[0][0],
                deltaTicks / (float)ticksPerSecond,
                1.0f * (AppIsKeyPressed(VK_LSHIFT) ? 3.0f : 1.0f) * activated,
                0.5f * activated,
                80.0f,
                deltaMouseX, deltaMouseY,
                AppIsKeyPressed('W'), AppIsKeyPressed('A'), AppIsKeyPressed('S'), AppIsKeyPressed('D'),
                AppIsKeyPressed(VK_SPACE), AppIsKeyPressed(VK_LCONTROL),
                FLYTHROUGH_CAMERA_LEFT_HANDED_BIT);
        else
            flythrough_camera_look_to(&g_Scene.CameraPos.x, &g_Scene.CameraLook.x, &cameraUp.x, &worldView.m[0][0], FLYTHROUGH_CAMERA_LEFT_HANDED_BIT);

        g_Scene.SceneWorldView = XMLoadFloat4x4(&worldView);

        D3D11_MAPPED_SUBRESOURCE mappedCamera;
        CHECKHR(dc->Map(g_Scene.pCameraBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedCamera));

        XMMATRIX worldViewProjection = XMMatrixMultiply(XMLoadFloat4x4(&worldView), g_Scene.SceneProjection);
        g_Scene.CameraWorldViewProjection = worldViewProjection;

        PerCameraData* camera = (PerCameraData*)mappedCamera.pData;
        XMStoreFloat4x4(&camera->WorldViewProjection, XMMatrixTranspose(worldViewProjection));
        XMStoreFloat4(&camera->WorldPosition, XMVectorSetW(XMLoadFloat3(&g_Scene.CameraPos), 1.0f));
        XMStoreFloat4(&camera->LookDirection, XMLoadFloat3(&g_Scene.CameraLook));
        XMStoreFloat4(&camera->Up, XMLoadFloat3(&cameraUp));

        dc->Unmap(g_Scene.pCameraBuffer.Get(), 0);
    }

    // Update viewport
    {
        D3D11_MAPPED_SUBRESOURCE mappedViewport;
        CHECKHR(dc->Map(g_Scene.pViewportBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedViewport));

        PerViewportData* viewport = (PerViewportData*)mappedViewport.pData;
        XMStoreFloat4(&viewport->Size, XMVectorSet(g_Scene.SceneViewport.Width, g_Scene.SceneViewport.Height, 0.0, 0.0));

        dc->Unmap(g_Scene.pViewportBuffer.Get(), 0);
    }

    // Update dragged vertex and mesh with it
    if (g_Scene.ModelingSceneNodeID != -1 && g_Scene.DragVertexID != -1)
    {
        SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
        StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

        if (staticMesh.SystemMatrixNeedsRebuild)
        {
            arap_factorize_system(
                (int)staticMesh.CPUPositions->size(),
                staticMesh.VertexConstraintStatuses->data(),
                (const int*)staticMesh.Halfedge->Vertices.data(),
                (const int*)staticMesh.Halfedge->Halfedges.data(),
                staticMesh.EdgeWeights->data(),
                staticMesh.ARAPSystemMatrix->data());

            staticMesh.SystemMatrixNeedsRebuild = false;
        }

        XMFLOAT3* cpuPosBuf = staticMesh.CPUPositions->data();
        XMFLOAT3 currPos = cpuPosBuf[g_Scene.DragVertexID];

        XMVECTOR up = XMLoadFloat3(&cameraUp);
        XMVECTOR forward = XMLoadFloat3(&g_Scene.CameraLook);
        XMVECTOR across = XMVector3Normalize(XMVector3Cross(forward, up));
        XMVECTOR upward = XMVector3Normalize(XMVector3Cross(across, forward));

        XMMATRIX transform = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, sceneNode.Transform.Translation);

        XMVECTOR clipPos = XMVector3Transform(XMLoadFloat3(&currPos), XMMatrixMultiply(XMMatrixMultiply(transform, g_Scene.SceneWorldView), g_Scene.SceneProjection));

        XMVECTOR oldUnprojected = XMVector3Unproject(
            XMVectorSet((float)g_Scene.LastMouseX, (float)g_Scene.LastMouseY, XMVectorGetZ(clipPos) / XMVectorGetW(clipPos), 1.0f),
            g_Scene.SceneViewport.TopLeftX, g_Scene.SceneViewport.TopLeftY, g_Scene.SceneViewport.Width, g_Scene.SceneViewport.Height, g_Scene.SceneViewport.MinDepth, g_Scene.SceneViewport.MaxDepth,
            g_Scene.SceneProjection, g_Scene.SceneWorldView, transform);

        XMVECTOR newUnprojected = XMVector3Unproject(
            XMVectorSet((float)currMouseX, (float)currMouseY, XMVectorGetZ(clipPos) / XMVectorGetW(clipPos), 1.0f),
            g_Scene.SceneViewport.TopLeftX, g_Scene.SceneViewport.TopLeftY, g_Scene.SceneViewport.Width, g_Scene.SceneViewport.Height, g_Scene.SceneViewport.MinDepth, g_Scene.SceneViewport.MaxDepth,
            g_Scene.SceneProjection, g_Scene.SceneWorldView, transform);

        XMVECTOR movement = newUnprojected - oldUnprojected;

        // set positions of control vertices
        for (int c = 0; c < (int)staticMesh.ControlVertexIDs->size(); c++)
        {
            int i = (*staticMesh.ControlVertexIDs)[c];
            XMStoreFloat3(&cpuPosBuf[i], XMLoadFloat3(&cpuPosBuf[i]) + movement);
        }

        D3D11_MAPPED_SUBRESOURCE mapped;
        CHECKHR(dc->Map(staticMesh.pPositionVertexBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped));

        // update the mesh as rigidly as possible based on the constrained vertices
        arap(
            (int)staticMesh.CPUPositions->size(), &staticMesh.BindPoseCPUPositions->data()->x, &staticMesh.CPUPositions->data()->x,
            staticMesh.VertexConstraintStatuses->data(),
            (const int*)staticMesh.Halfedge->Vertices.data(),
            (const int*)staticMesh.Halfedge->Halfedges.data(),
            staticMesh.EdgeWeights->data(),
            staticMesh.ARAPSystemMatrix->data(),
            DEFAULT_NUM_ARAP_ITERATIONS);

        memcpy(mapped.pData, cpuPosBuf, staticMesh.CPUPositions->size() * sizeof(XMFLOAT3));

        dc->Unmap(staticMesh.pPositionVertexBuffer.Get(), 0);
    }

    const float kClearColor[] = {
        std::pow(100.0f / 255.0f, 2.2f),
        std::pow(149.0f / 255.0f, 2.2f),
        std::pow(237.0f / 255.0f, 2.2f),
        1.0f
    };
    dc->ClearRenderTargetView(pBackBufferRTV, kClearColor);
    dc->ClearDepthStencilView(g_Scene.pSceneDepthDSV.Get(), D3D11_CLEAR_DEPTH, 1.0f, 0);

    // Scene rendering pass
    {
        Material defaultMtl = {};
        XMStoreFloat3(&defaultMtl.Ambient, XMVectorSet(0.1f, 0.1f, 0.1f, 0.0f));
        XMStoreFloat3(&defaultMtl.Diffuse, XMVectorSet(0.9f, 0.9f, 0.9f, 0.0f));
        XMStoreFloat3(&defaultMtl.Specular, XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f));
        defaultMtl.Shininess = 1.0f;
        defaultMtl.Opacity = 1.0f;
        defaultMtl.DiffuseTextureID = -1;
        defaultMtl.SpecularTextureID = -1;
        defaultMtl.BumpTextureID = -1;

        ID3D11RenderTargetView* sceneRTVs[] = { pBackBufferRTV, g_Scene.pSceneNormalRTV.Get() };
        ID3D11DepthStencilView* sceneDSV = g_Scene.pSceneDepthDSV.Get();
        dc->OMSetRenderTargets(_countof(sceneRTVs), sceneRTVs, sceneDSV);

        dc->VSSetShader(g_Scene.SceneVS->VS, NULL, 0);
        dc->GSSetShader(g_Scene.SceneGS->GS, NULL, 0);
        dc->PSSetShader(g_Scene.ScenePS->PS, NULL, 0);
        dc->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        dc->IASetInputLayout(g_Scene.pSceneInputLayout.Get());
        dc->RSSetState(g_Scene.pSceneRasterizerState.Get());
        dc->OMSetDepthStencilState(g_Scene.pSceneDepthStencilState.Get(), 0);
        dc->OMSetBlendState(g_Scene.pSceneBlendState.Get(), NULL, UINT_MAX);
        dc->RSSetViewports(1, &g_Scene.SceneViewport);

        ID3D11Buffer* cameraCBV = g_Scene.pCameraBuffer.Get();
        dc->VSSetConstantBuffers(SCENE_CAMERA_BUFFER_SLOT, 1, &cameraCBV);
        dc->PSSetConstantBuffers(SCENE_CAMERA_BUFFER_SLOT, 1, &cameraCBV);
        
        ID3D11SamplerState* diffuseSMP = g_Scene.pDiffuseSampler.Get();
        dc->PSSetSamplers(SCENE_DIFFUSE_SAMPLER_SLOT, 1, &diffuseSMP);

        ID3D11SamplerState* specularSMP = g_Scene.pSpecularSampler.Get();
        dc->PSSetSamplers(SCENE_SPECULAR_SAMPLER_SLOT, 1, &specularSMP);

        ID3D11SamplerState* bumpSMP = g_Scene.pBumpSampler.Get();
        dc->PSSetSamplers(SCENE_BUMP_SAMPLER_SLOT, 1, &bumpSMP);

        int currMaterialID = -2;
        for (int sceneNodeID = 0; sceneNodeID < (int)g_Scene.SceneNodes.size(); sceneNodeID++)
        {
            if (sceneNodeID != g_Scene.ModelingSceneNodeID)
                continue;

            SceneNode& sceneNode = g_Scene.SceneNodes[sceneNodeID];

            // Update Material CBV
            if (currMaterialID != sceneNode.MaterialID)
            {
                const Material* material;

                if (sceneNode.MaterialID != -1)
                {
                    material = &g_Scene.Materials[sceneNode.MaterialID];
                    currMaterialID = sceneNode.MaterialID;
                }
                else
                {
                     material = &defaultMtl;
                     currMaterialID = -1;
                }

                D3D11_MAPPED_SUBRESOURCE mapped;
                dc->Map(g_Scene.pMaterialBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);

                PerMaterialData* materialData = (PerMaterialData*)mapped.pData;
                XMStoreFloat4(&materialData->Ambient, XMLoadFloat3(&material->Ambient));
                XMStoreFloat4(&materialData->Diffuse, XMLoadFloat3(&material->Diffuse));
                XMStoreFloat4(&materialData->Specular, XMLoadFloat3(&material->Specular));
                XMStoreFloat4(&materialData->Shininess, XMVectorReplicate(material->Shininess));
                XMStoreFloat4(&materialData->Opacity, XMVectorReplicate(material->Opacity));
                XMStoreFloat4(&materialData->HasDiffuse, XMVectorReplicate(material->DiffuseTextureID != -1 ? 1.0f : 0.0f));
                XMStoreFloat4(&materialData->HasSpecular, XMVectorReplicate(material->SpecularTextureID != -1 ? 1.0f : 0.0f));
                XMStoreFloat4(&materialData->HasBump, XMVectorReplicate(material->BumpTextureID != -1 ? 1.0f : 0.0f));

                dc->Unmap(g_Scene.pMaterialBuffer.Get(), 0);

                ID3D11Buffer* materialCBV = g_Scene.pMaterialBuffer.Get();
                dc->PSSetConstantBuffers(SCENE_MATERIAL_BUFFER_SLOT, 1, &materialCBV);

                ID3D11ShaderResourceView* diffuseSRV = NULL;
                if (material->DiffuseTextureID != -1)
                diffuseSRV = g_Scene.Textures[material->DiffuseTextureID].SRV.Get();
                dc->PSSetShaderResources(SCENE_DIFFUSE_TEXTURE_SLOT, 1, &diffuseSRV);

                ID3D11ShaderResourceView* specularSRV = NULL;
                if (material->SpecularTextureID != -1)
                specularSRV = g_Scene.Textures[material->SpecularTextureID].SRV.Get();
                dc->PSSetShaderResources(SCENE_SPECULAR_TEXTURE_SLOT, 1, &specularSRV);

                ID3D11ShaderResourceView* bumpSRV = NULL;
                if (material->BumpTextureID != -1)
                bumpSRV = g_Scene.Textures[material->BumpTextureID].SRV.Get();
                dc->PSSetShaderResources(SCENE_BUMP_TEXTURE_SLOT, 1, &bumpSRV);
            }

            // Update SceneNode CBV
            {
                D3D11_MAPPED_SUBRESOURCE mapped;
                dc->Map(g_Scene.pSceneNodeBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);

                PerSceneNodeData* sceneNodeData = (PerSceneNodeData*)mapped.pData;

                XMMATRIX worldMatrix = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, sceneNode.Transform.Translation);
                XMStoreFloat4x4(&sceneNodeData->WorldTransform, XMMatrixTranspose(worldMatrix));

                XMMATRIX normalMatrix = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, XMVectorZero());
                XMStoreFloat4x4(&sceneNodeData->NormalTransform, XMMatrixTranspose(normalMatrix));

                if (sceneNode.Type == SCENENODETYPE_STATICMESH &&
                    !g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID].pNormalVertexBuffer)
                {
                    XMStoreFloat4(&sceneNodeData->AutoGenNormals, XMVectorReplicate(1.0f));
                }
                else
                {
                    XMStoreFloat4(&sceneNodeData->AutoGenNormals, XMVectorReplicate(0.0f));
                }

                dc->Unmap(g_Scene.pSceneNodeBuffer.Get(), 0);

                ID3D11Buffer* sceneNodeCBV = g_Scene.pSceneNodeBuffer.Get();
                dc->VSSetConstantBuffers(SCENE_SCENENODE_BUFFER_SLOT, 1, &sceneNodeCBV);
                dc->GSSetConstantBuffers(SCENE_SCENENODE_BUFFER_SLOT, 1, &sceneNodeCBV);
            }

            if (sceneNode.Type == SCENENODETYPE_STATICMESH)
            {
                StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

                ID3D11Buffer* staticMeshVertexBuffers[] = {
                    staticMesh.pPositionVertexBuffer.Get(),
                    staticMesh.pTexCoordVertexBuffer.Get(),
                    staticMesh.pNormalVertexBuffer.Get(),
                    staticMesh.pTangentVertexBuffer.Get(),
                    staticMesh.pBitangentVertexBuffer.Get()
                };
                UINT staticMeshStrides[] = {
                    sizeof(VertexPosition),
                    sizeof(VertexTexCoord),
                    sizeof(VertexNormal),
                    sizeof(VertexTangent),
                    sizeof(VertexBitangent)
                };
                UINT staticMeshOffsets[] = {
                    0, 0, 0, 0, 0
                };
                dc->IASetVertexBuffers(0, _countof(staticMeshVertexBuffers), staticMeshVertexBuffers, staticMeshStrides, staticMeshOffsets);
                dc->IASetIndexBuffer(staticMesh.pIndexBuffer.Get(), DXGI_FORMAT_R32_UINT, 0);

                dc->DrawIndexed(staticMesh.IndexCountPerInstance, staticMesh.StartIndexLocation, 0);
            }
        }

        dc->OMSetRenderTargets(0, NULL, NULL);
        dc->VSSetShader(NULL, NULL, 0);
        dc->GSSetShader(NULL, NULL, 0);
        dc->PSSetShader(NULL, NULL, 0);
    }

    // update current selection
    {
        D3D11_MAPPED_SUBRESOURCE mapped;
        CHECKHR(dc->Map(g_Scene.pCurrSelectedVertexID.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped));

        CurrSelectionData selection = {};
        selection.VertexID.x = selection.VertexID.y = selection.VertexID.z = selection.VertexID.w = g_Scene.DragVertexID == UINT_MAX ? PickSelector(currMouseX, currMouseY) : g_Scene.DragVertexID;

        CurrSelectionData* pSelection = (CurrSelectionData*)mapped.pData;
        *pSelection = selection;

        dc->Unmap(g_Scene.pCurrSelectedVertexID.Get(), 0);
    }

    UINT kSelectorClear[4] = { UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX };
    dc->ClearUnorderedAccessViewUint(g_Scene.pSelectorUAV.Get(), kSelectorClear);

    // Draw selection markers
    if (g_Scene.ModelingSceneNodeID != -1)
    {
        ID3D11RenderTargetView* selectorRTVs[] = { pBackBufferRTV, g_Scene.pSceneNormalRTV.Get(), g_Scene.pSelectorRTV.Get() };
        ID3D11DepthStencilView* selectorDSV = g_Scene.pSceneDepthDSV.Get();
        dc->OMSetRenderTargets(_countof(selectorRTVs), selectorRTVs, selectorDSV);

        const SceneNode& sceneNode = g_Scene.SceneNodes[g_Scene.ModelingSceneNodeID];
        const StaticMesh& staticMesh = g_Scene.StaticMeshes[sceneNode.AsStaticMesh.StaticMeshID];

        dc->IASetInputLayout(g_Scene.pSelectorSphereInputLayout.Get());
        dc->RSSetState(g_Scene.pSelectorRasterizerState.Get());
        dc->OMSetDepthStencilState(g_Scene.pSelectorDepthStencilState.Get(), 0);
        dc->OMSetBlendState(g_Scene.pSelectorBlendState.Get(), NULL, UINT_MAX);
        dc->RSSetViewports(1, &g_Scene.SceneViewport);

        dc->VSSetShader(g_Scene.SelectorVS->VS, NULL, 0);
        dc->GSSetShader(g_Scene.SelectorGS->GS, NULL, 0);
        dc->PSSetShader(g_Scene.SelectorPS->PS, NULL, 0);

        ID3D11Buffer* cameraCBV = g_Scene.pCameraBuffer.Get();
        dc->GSSetConstantBuffers(SELECTOR_CAMERA_BUFFER_SLOT, 1, &cameraCBV);
        dc->PSSetConstantBuffers(SELECTOR_CAMERA_BUFFER_SLOT, 1, &cameraCBV);

        ID3D11Buffer* currSelectedCBV = g_Scene.pCurrSelectedVertexID.Get();
        dc->GSSetConstantBuffers(SELECTOR_CURRSELECTION_BUFFER_SLOT, 1, &currSelectedCBV);

        // SceneNode CBV
        {
            D3D11_MAPPED_SUBRESOURCE mapped;
            dc->Map(g_Scene.pSceneNodeBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);

            PerSceneNodeData* sceneNodeData = (PerSceneNodeData*)mapped.pData;

            XMMATRIX worldMatrix = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, sceneNode.Transform.Translation);
            XMStoreFloat4x4(&sceneNodeData->WorldTransform, XMMatrixTranspose(worldMatrix));

            XMMATRIX normalMatrix = XMMatrixAffineTransformation(sceneNode.Transform.Scaling, sceneNode.Transform.RotationOrigin, sceneNode.Transform.RotationQuaternion, XMVectorZero());
            XMStoreFloat4x4(&sceneNodeData->NormalTransform, XMMatrixTranspose(normalMatrix));

            dc->Unmap(g_Scene.pSceneNodeBuffer.Get(), 0);

            ID3D11Buffer* sceneNodeCBV = g_Scene.pSceneNodeBuffer.Get();
            dc->GSSetConstantBuffers(SELECTOR_SCENENODE_BUFFER_SLOT, 1, &sceneNodeCBV);
        }

        ID3D11Buffer* selectorBuffers[] = { staticMesh.pPositionVertexBuffer.Get(), staticMesh.pSelectorColorVertexBuffer.Get() };
        UINT selectorStrides[] = { sizeof(VertexPosition), sizeof(VertexSelectorColor) };
        UINT selectorOffsets[] = { 0, 0 };
        dc->IASetVertexBuffers(0, _countof(selectorBuffers), selectorBuffers, selectorStrides, selectorOffsets);
        dc->IASetIndexBuffer(staticMesh.pIndexBuffer.Get(), DXGI_FORMAT_R32_UINT, 0);
        dc->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
        dc->DrawIndexed(staticMesh.IndexCountPerInstance, staticMesh.StartIndexLocation, 0);

        dc->OMSetRenderTargets(0, NULL, NULL);

        dc->VSSetShader(NULL, NULL, 0);
        dc->GSSetShader(NULL, NULL, 0);
        dc->PSSetShader(NULL, NULL, 0);
    }

    // SSAO pass
    {
        ID3D11RenderTargetView* ssaoRTVs[] = { pBackBufferRTV };
        dc->OMSetRenderTargets(_countof(ssaoRTVs), ssaoRTVs, NULL);

        dc->VSSetShader(g_Scene.SSAOVS->VS, NULL, 0);
        dc->PSSetShader(g_Scene.SSAOPS->PS, NULL, 0);
        dc->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        dc->IASetInputLayout(NULL);
        dc->RSSetState(g_Scene.pSSAORasterizerState.Get());
        dc->OMSetDepthStencilState(g_Scene.pSSAODepthStencilState.Get(), 0);
        dc->OMSetBlendState(g_Scene.pSSAOBlendState.Get(), NULL, UINT_MAX);
        dc->RSSetViewports(1, &g_Scene.SceneViewport);

        ID3D11ShaderResourceView* normalsSRV = g_Scene.pSceneNormalSRV.Get();
        dc->PSSetShaderResources(SSAO_NORMAL_TEXTURE_SLOT, 1, &normalsSRV);

        ID3D11SamplerState* normalsSMP = g_Scene.pSceneNormalSampler.Get();
        dc->PSSetSamplers(SSAO_NORMAL_SAMPLER_SLOT, 1, &normalsSMP);

        ID3D11ShaderResourceView* depthSRV = g_Scene.pSceneDepthSRV.Get();
        dc->PSSetShaderResources(SSAO_DEPTH_TEXTURE_SLOT, 1, &depthSRV);

        ID3D11SamplerState* depthSMP = g_Scene.pSceneDepthSampler.Get();
        dc->PSSetSamplers(SSAO_DEPTH_SAMPLER_SLOT, 1, &depthSMP);

        ID3D11ShaderResourceView* noiseSRV = g_Scene.pSSAONoiseSRV.Get();
        dc->PSSetShaderResources(SSAO_NOISE_TEXTURE_SLOT, 1, &noiseSRV);

        ID3D11SamplerState* noiseSMP = g_Scene.pSSAONoiseSampler.Get();
        dc->PSSetSamplers(SSAO_NOISE_SAMPLER_SLOT, 1, &noiseSMP);

        dc->IASetVertexBuffers(0, 0, NULL, NULL, NULL);
        dc->IASetIndexBuffer(NULL, DXGI_FORMAT_UNKNOWN, 0);

        dc->Draw(3, 0);

        ID3D11ShaderResourceView* nullSRV = NULL;
        dc->PSSetShaderResources(SSAO_NORMAL_TEXTURE_SLOT, 1, &nullSRV);
        dc->PSSetShaderResources(SSAO_DEPTH_TEXTURE_SLOT, 1, &nullSRV);
        dc->PSSetShaderResources(SSAO_NOISE_TEXTURE_SLOT, 1, &nullSRV);

        dc->OMSetRenderTargets(0, NULL, NULL);
        dc->VSSetShader(NULL, NULL, 0);
        dc->PSSetShader(NULL, NULL, 0);
    }

    // Draw ROI selection
    if (!g_Scene.ROISelectionPoints.empty())
    {
        ID3D11RenderTargetView* roiRTVs[] = { pBackBufferRTV };
        dc->OMSetRenderTargets(_countof(roiRTVs), roiRTVs, NULL);

        dc->VSSetShader(g_Scene.ROIVS->VS, NULL, 0);
        dc->GSSetShader(g_Scene.ROIGS->GS, NULL, 0);
        dc->PSSetShader(g_Scene.ROIPS->PS, NULL, 0);
        dc->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINESTRIP);
        dc->IASetInputLayout(g_Scene.pROIInputLayout.Get());
        dc->RSSetState(g_Scene.pROIRasterizerState.Get());
        dc->OMSetDepthStencilState(g_Scene.pROIDepthStencilState.Get(), 0);
        dc->OMSetBlendState(g_Scene.pROIBlendState.Get(), NULL, UINT_MAX);
        dc->RSSetViewports(1, &g_Scene.SceneViewport);

        ID3D11Buffer* viewportBuffer = g_Scene.pViewportBuffer.Get();
        dc->GSSetConstantBuffers(ROI_VIEWPORT_BUFFER_SLOT, 1, &viewportBuffer);

        int sizeRequired = (int)g_Scene.ROISelectionPoints.size() * sizeof(XMFLOAT2);
        if (g_Scene.ROISelectorBufferSizeInBytes < sizeRequired)
        {
            g_Scene.ROISelectorBufferSizeInBytes = sizeRequired * 2;

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC((UINT)g_Scene.ROISelectorBufferSizeInBytes, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_DYNAMIC, D3D11_CPU_ACCESS_WRITE),
                NULL,
                &g_Scene.pROISelectorBuffer));
        }

        D3D11_MAPPED_SUBRESOURCE mappedVertices;
        CHECKHR(dc->Map(g_Scene.pROISelectorBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedVertices));

        memcpy(mappedVertices.pData, g_Scene.ROISelectionPoints.data(), sizeRequired);

        dc->Unmap(g_Scene.pROISelectorBuffer.Get(), 0);
        
        ID3D11Buffer* roiVertexBuffers[] = { g_Scene.pROISelectorBuffer.Get() };
        UINT roiVertexStrides[] = { sizeof(XMFLOAT2) };
        UINT roiVertexOffsets[] = { 0 };
        dc->IASetVertexBuffers(0, _countof(roiVertexBuffers), roiVertexBuffers, roiVertexStrides, roiVertexOffsets);
        dc->IASetIndexBuffer(NULL, DXGI_FORMAT_UNKNOWN, 0);

        dc->Draw((UINT)g_Scene.ROISelectionPoints.size(), 0);

        dc->OMSetRenderTargets(0, NULL, NULL);
        dc->VSSetShader(NULL, NULL, 0);
        dc->GSSetShader(NULL, NULL, 0);
        dc->PSSetShader(NULL, NULL, 0);
    }

    g_Scene.LastTicks = currTicks;
    g_Scene.LastMouseX = currMouseX;
    g_Scene.LastMouseY = currMouseY;
}