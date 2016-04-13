#ifndef COMMON_HLSL
#define COMMON_HLSL

struct PerCameraData
{
    float4x4 WorldViewProjection;
    float4 WorldPosition;
    float4 LookDirection;
    float4 Up;
};

struct PerViewportData
{
    float4 Size;
};

struct PerMaterialData
{
    float4 Ambient;
    float4 Diffuse;
    float4 Specular;
    float4 Shininess; // all same
    float4 Opacity; // all same

    float4 HasDiffuse; // all same
    float4 HasSpecular; // all same
    float4 HasBump; // all same
};

struct PerSceneNodeData
{
    float4x4 WorldTransform;
    float4x4 NormalTransform;
    float4 AutoGenNormals;
};

struct CurrSelectionData
{
    uint4 VertexID;
    uint4 Captured;
};

#define SCENE_CAMERA_BUFFER_SLOT 0
#define SCENE_MATERIAL_BUFFER_SLOT 1
#define SCENE_SCENENODE_BUFFER_SLOT 2

#define SCENE_DIFFUSE_TEXTURE_SLOT 0
#define SCENE_SPECULAR_TEXTURE_SLOT 1
#define SCENE_BUMP_TEXTURE_SLOT 2

#define SCENE_DIFFUSE_SAMPLER_SLOT 0
#define SCENE_SPECULAR_SAMPLER_SLOT 1
#define SCENE_BUMP_SAMPLER_SLOT 2

#define SSAO_NORMAL_TEXTURE_SLOT 0
#define SSAO_DEPTH_TEXTURE_SLOT 1
#define SSAO_NOISE_TEXTURE_SLOT 2

#define SSAO_NORMAL_SAMPLER_SLOT 0
#define SSAO_DEPTH_SAMPLER_SLOT 1
#define SSAO_NOISE_SAMPLER_SLOT 2

#define SELECTOR_CAMERA_BUFFER_SLOT 0
#define SELECTOR_SCENENODE_BUFFER_SLOT 1
#define SELECTOR_CURRSELECTION_BUFFER_SLOT 2

#define BUFFER_REGISTER(slot) b##slot
#define TEXTURE_REGISTER(slot) t##slot
#define SAMPLER_REGISTER(slot) s##slot

#endif // COMMON_HLSL