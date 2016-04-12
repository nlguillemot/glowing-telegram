#include "scene.h"

#include "apputil.h"
#include "renderer.h"
#include "app.h"
#include "flythrough_camera.h"

#include "imgui.h"
#include "tiny_obj_loader.h"
#include "stb_image.h"

#include <DirectXMath.h>

using namespace DirectX;

// to make HLSL compile as C++
using float4 = XMFLOAT4;
using float4x4 = XMFLOAT4X4;

#include "shaders/common.hlsl"

#include <unordered_map>
#include <random>

struct VertexPosition
{
    XMFLOAT3 Position;
};

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
    ComPtr<ID3D11Buffer> pTexCoordVertexBuffer;
    ComPtr<ID3D11Buffer> pNormalVertexBuffer;
    ComPtr<ID3D11Buffer> pTangentVertexBuffer;
    ComPtr<ID3D11Buffer> pBitangentVertexBuffer;
    ComPtr<ID3D11Buffer> pIndexBuffer;
    int MaterialID; // the material this mesh was designed for
    UINT IndexCountPerInstance;
    UINT StartIndexLocation;
};

struct NodeTransform
{
    XMVECTOR Scale;
    XMVECTOR Quaternion;
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
    ComPtr<ID3D11Buffer> pCameraBuffer;
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
    Shader* ScenePS;

    ComPtr<ID3D11Texture2D> pSSAONoiseTexture;
    ComPtr<ID3D11ShaderResourceView> pSSAONoiseSRV;
    ComPtr<ID3D11SamplerState> pSSAONoiseSampler;

    ComPtr<ID3D11RasterizerState> pSSAORasterizerState;
    ComPtr<ID3D11DepthStencilState> pSSAODepthStencilState;
    ComPtr<ID3D11BlendState> pSSAOBlendState;

    Shader* SSAOVS;
    Shader* SSAOPS;

    uint64_t LastTicks;
    int LastMouseX, LastMouseY;
};

Scene g_Scene;

static void SceneAddObjMesh(
    const char* filename, const char* mtlbasepath,
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
        ComPtr<ID3D11Buffer> pTexCoordBuffer;
        ComPtr<ID3D11Buffer> pNormalBuffer;
        ComPtr<ID3D11Buffer> pTangentBuffer;
        ComPtr<ID3D11Buffer> pBitangentBuffer;
        ComPtr<ID3D11Buffer> pIndexBuffer;

        int numVertices = (int)mesh.positions.size() / 3;

        if (!mesh.positions.empty())
        {
            D3D11_SUBRESOURCE_DATA positionVertexBufferData = {};
            positionVertexBufferData.pSysMem = mesh.positions.data();

            CHECKHR(dev->CreateBuffer(
                &CD3D11_BUFFER_DESC(sizeof(VertexPosition) * numVertices, D3D11_BIND_VERTEX_BUFFER, D3D11_USAGE_IMMUTABLE), 
                &positionVertexBufferData, 
                &pPositionBuffer));
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

        if (mesh.indices.empty())
        {
            SimpleMessageBox_FatalError("Expected indices");
        }

        static_assert(sizeof(mesh.indices[0]) == sizeof(UINT32), "Expecting UINT32 indices");

        D3D11_SUBRESOURCE_DATA indexBufferData = {};
        indexBufferData.pSysMem = mesh.indices.data();

        CHECKHR(dev->CreateBuffer(
            &CD3D11_BUFFER_DESC(sizeof(UINT32) * numIndices, D3D11_BIND_INDEX_BUFFER, D3D11_USAGE_IMMUTABLE), 
            &indexBufferData, 
            &pIndexBuffer));

        const int numFaces = (int)shape.mesh.indices.size() / 3;

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

            StaticMesh sm;
            sm.Name = shape.name;
            sm.pPositionVertexBuffer = pPositionBuffer;
            sm.pTexCoordVertexBuffer = pTexCoordBuffer;
            sm.pNormalVertexBuffer = pNormalBuffer;
            sm.pTangentVertexBuffer = pTangentBuffer;
            sm.pBitangentVertexBuffer = pBitangentBuffer;
            sm.pIndexBuffer = pIndexBuffer;
            sm.MaterialID = firstMaterial + currMTL;
            sm.IndexCountPerInstance = (face + 1 - firstFace) * 3;
            sm.StartIndexLocation = firstFace * 3;

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
    sceneNode.Transform.Scale = XMVectorSet(1.0f, 1.0f, 1.0f, 0.0f);
    sceneNode.Transform.Quaternion = XMQuaternionIdentity();
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
        "dragon"
    };

    std::vector<int> newStaticMeshIDs;
    for (const std::string& meshToLoad : meshesToLoad)
    {
        std::string meshFolder = "assets/" + meshToLoad + "/";
        std::string meshFile = meshFolder + meshToLoad + ".obj";
        SceneAddObjMesh(meshFile.c_str(), meshFolder.c_str(), &newStaticMeshIDs);
    }

    for (int newStaticMeshID : newStaticMeshIDs)
    {
        SceneAddStaticMeshSceneNode(newStaticMeshID);
    }

    g_Scene.SceneVS = RendererAddShader("scene.hlsl", "VSmain", "vs_5_0");
    g_Scene.ScenePS = RendererAddShader("scene.hlsl", "PSmain", "ps_5_0");
    g_Scene.SSAOVS = RendererAddShader("ssao.hlsl", "VSmain", "vs_5_0");
    g_Scene.SSAOPS = RendererAddShader("ssao.hlsl", "PSmain", "ps_5_0");

    XMStoreFloat3(&g_Scene.CameraPos, XMVectorSet(-0.5f, 0.5f, -0.5f, 1.0f));
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

    D3D11_INPUT_ELEMENT_DESC sceneInputElements[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 1, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 2, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "TANGENT", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 3, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "BITANGENT", 0, DXGI_FORMAT_R32G32B32_FLOAT, 4, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
    };
    CHECKHR(dev->CreateInputLayout(
        sceneInputElements, _countof(sceneInputElements),
        g_Scene.SceneVS->Blob->GetBufferPointer(), g_Scene.SceneVS->Blob->GetBufferSize(),
        &g_Scene.pSceneInputLayout));

    D3D11_RASTERIZER_DESC sceneRasterizerDesc = CD3D11_RASTERIZER_DESC(D3D11_DEFAULT);
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
        std::vector<XMFLOAT3> noiseData(noiseSize * noiseSize);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<float> d(0.0f, 1.0f);
        for (XMFLOAT3& noise : noiseData)
        {
            noise.x = d(gen);
            noise.y = d(gen);
            noise.z = d(gen);
        }
        
        D3D11_SUBRESOURCE_DATA noiseInitialData = {};
        noiseInitialData.pSysMem = noiseData.data();
        noiseInitialData.SysMemPitch = noiseSize * sizeof(XMFLOAT3);
        noiseInitialData.SysMemSlicePitch = noiseSize * noiseSize * sizeof(XMFLOAT3);
        
        CHECKHR(dev->CreateTexture2D(
            &CD3D11_TEXTURE2D_DESC(DXGI_FORMAT_R32G32B32_FLOAT, noiseSize, noiseSize, 1, 1, D3D11_BIND_SHADER_RESOURCE),
            &noiseInitialData,
            &g_Scene.pSSAONoiseTexture));

        CHECKHR(dev->CreateShaderResourceView(
            g_Scene.pSSAONoiseTexture.Get(),
            &CD3D11_SHADER_RESOURCE_VIEW_DESC(D3D11_SRV_DIMENSION_TEXTURE2D, DXGI_FORMAT_R32G32B32_FLOAT, 0, 1),
            &g_Scene.pSSAONoiseSRV));

        CD3D11_SAMPLER_DESC noiseSamplerDesc(D3D11_DEFAULT);
        noiseSamplerDesc.Filter = D3D11_FILTER_MIN_MAG_LINEAR_MIP_POINT;
        noiseSamplerDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
        noiseSamplerDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
        CHECKHR(dev->CreateSamplerState(&noiseSamplerDesc, &g_Scene.pSSAONoiseSampler));
    }

    g_Scene.LastMouseX = INT_MIN;
    g_Scene.LastMouseY = INT_MIN;
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

    g_Scene.SceneViewport = CD3D11_VIEWPORT(0.0f, 0.0f, (FLOAT)renderWidth, (FLOAT)renderHeight);
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

    ID3D11Device* dev = RendererGetDevice();
    ID3D11DeviceContext* dc = RendererGetDeviceContext();

    // Update camera
    {
        float activated = GetAsyncKeyState(VK_RBUTTON) ? 1.0f : 0.0f;
        float up[3] = { 0.0f, 1.0f, 0.0f };
        XMFLOAT4X4 worldView;
        flythrough_camera_update(
            &g_Scene.CameraPos.x,
            &g_Scene.CameraLook.x,
            up,
            &worldView.m[0][0],
            deltaTicks / (float)ticksPerSecond,
            1.0f * (GetAsyncKeyState(VK_LSHIFT) ? 3.0f : 1.0f) * activated,
            0.5f * activated,
            80.0f,
            currMouseX - g_Scene.LastMouseX, currMouseY - g_Scene.LastMouseY,
            GetAsyncKeyState('W'), GetAsyncKeyState('A'), GetAsyncKeyState('S'), GetAsyncKeyState('D'),
            GetAsyncKeyState(VK_SPACE), GetAsyncKeyState(VK_LCONTROL),
            FLYTHROUGH_CAMERA_LEFT_HANDED_BIT);

        D3D11_MAPPED_SUBRESOURCE mappedCamera;
        CHECKHR(dc->Map(g_Scene.pCameraBuffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedCamera));

        float aspectWbyH = g_Scene.SceneViewport.Width / g_Scene.SceneViewport.Height;
        XMMATRIX viewProjection = XMMatrixPerspectiveFovLH(XMConvertToRadians(90.0f), aspectWbyH, 0.0001f, 20.0f);
        XMMATRIX worldViewProjection = XMMatrixMultiply(XMLoadFloat4x4(&worldView), viewProjection);

        PerCameraData* camera = (PerCameraData*)mappedCamera.pData;
        XMStoreFloat4x4(&camera->WorldViewProjection, XMMatrixTranspose(worldViewProjection));
        XMStoreFloat4(&camera->WorldPosition, XMVectorSetW(XMLoadFloat3(&g_Scene.CameraPos),1.0f));

        dc->Unmap(g_Scene.pCameraBuffer.Get(), 0);
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
        XMStoreFloat3(&defaultMtl.Ambient, XMVectorSet(0.5f, 0.5f, 0.5f, 0.0f));
        XMStoreFloat3(&defaultMtl.Diffuse, XMVectorSet(0.5f, 0.5f, 0.5f, 0.0f));
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

                XMMATRIX worldMatrix = XMMatrixIdentity();
                worldMatrix = XMMatrixMultiply(worldMatrix, XMMatrixScalingFromVector(sceneNode.Transform.Scale));
                worldMatrix = XMMatrixMultiply(worldMatrix, XMMatrixRotationQuaternion(sceneNode.Transform.Quaternion));
                worldMatrix = XMMatrixMultiply(worldMatrix, XMMatrixTranslationFromVector(sceneNode.Transform.Translation));
                XMStoreFloat4x4(&sceneNodeData->WorldTransform, XMMatrixTranspose(worldMatrix));

                XMMATRIX normalMatrix = XMMatrixIdentity();
                normalMatrix = XMMatrixMultiply(normalMatrix, XMMatrixScalingFromVector(XMVectorReciprocal(sceneNode.Transform.Scale)));
                normalMatrix = XMMatrixMultiply(normalMatrix, XMMatrixRotationQuaternion(sceneNode.Transform.Quaternion));
                XMStoreFloat4x4(&sceneNodeData->NormalTransform, XMMatrixTranspose(normalMatrix));

                dc->Unmap(g_Scene.pSceneNodeBuffer.Get(), 0);

                ID3D11Buffer* sceneNodeCBV = g_Scene.pSceneNodeBuffer.Get();
                dc->VSSetConstantBuffers(SCENE_SCENENODE_BUFFER_SLOT, 1, &sceneNodeCBV);
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
        dc->IASetIndexBuffer(NULL, (DXGI_FORMAT)0, 0);

        dc->Draw(3, 0);

        ID3D11ShaderResourceView* nullSRV = NULL;
        dc->PSSetShaderResources(SSAO_NORMAL_TEXTURE_SLOT, 1, &nullSRV);
        dc->PSSetShaderResources(SSAO_DEPTH_TEXTURE_SLOT, 1, &nullSRV);
        dc->PSSetShaderResources(SSAO_NOISE_TEXTURE_SLOT, 1, &nullSRV);

        dc->OMSetRenderTargets(0, NULL, NULL);
    }

    g_Scene.LastTicks = currTicks;
    g_Scene.LastMouseX = currMouseX;
    g_Scene.LastMouseY = currMouseY;
}