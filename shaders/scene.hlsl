#include "common.hlsl"

struct VSIn
{
    float4 Position : POSITION;
    float2 TexCoord : TEXCOORD;
    float3 Normal : NORMAL;
    float4 Tangent : TANGENT;
    float3 Bitangent : BITANGENT;
};

struct VSOut
{
    float4 Position : SV_Position;
    float3 WorldPosition : WORLDPOSITION;
    float2 TexCoord : TEXCOORD;
    float3 WorldNormal : WORLDNORMAL;
    float4 WorldTangent : WORLDTANGENT;
    float3 WorldBitangent : WORLDBITANGENT;
};

struct GSOut
{
    VSOut vs;

    nointerpolation float4 EdgePlane0 : EDGEPLANE0;
    nointerpolation float4 EdgePlane1 : EDGEPLANE1;
    nointerpolation float4 EdgePlane2 : EDGEPLANE2;
};

struct PSOut
{
    float4 Color : SV_Target0;
    float3 Normal : SV_Target1;
};

cbuffer CameraBuffer : register(BUFFER_REGISTER(SCENE_CAMERA_BUFFER_SLOT))
{
    PerCameraData Camera;
};

cbuffer MaterialBuffer : register(BUFFER_REGISTER(SCENE_MATERIAL_BUFFER_SLOT))
{
    PerMaterialData Material;
};

cbuffer SceneNodeBuffer : register(BUFFER_REGISTER(SCENE_SCENENODE_BUFFER_SLOT))
{
    PerSceneNodeData SceneNode;
};

Texture2D DiffuseTexture : register(TEXTURE_REGISTER(SCENE_DIFFUSE_TEXTURE_SLOT));
SamplerState DiffuseSampler : register(SAMPLER_REGISTER(SCENE_DIFFUSE_SAMPLER_SLOT));

Texture2D SpecularTexture : register(TEXTURE_REGISTER(SCENE_SPECULAR_TEXTURE_SLOT));
SamplerState SpecularSampler : register(SAMPLER_REGISTER(SCENE_SPECULAR_SAMPLER_SLOT));

Texture2D BumpTexture : register(TEXTURE_REGISTER(SCENE_BUMP_TEXTURE_SLOT));
SamplerState BumpSampler : register(SAMPLER_REGISTER(SCENE_BUMP_SAMPLER_SLOT));

VSOut VSmain(VSIn input)
{
    VSOut output;
    output.WorldPosition = mul(input.Position, SceneNode.WorldTransform).xyz;
    output.Position = mul(float4(output.WorldPosition, 1), Camera.WorldViewProjection);
    output.TexCoord = input.TexCoord;
    output.WorldNormal = normalize(mul(float4(input.Normal, 0), SceneNode.NormalTransform).xyz);
    output.WorldTangent = float4(normalize(mul(float4(input.Tangent.xyz, 0), SceneNode.NormalTransform).xyz), input.Tangent.w);
    output.WorldBitangent = normalize(mul(float4(input.Bitangent, 0), SceneNode.NormalTransform).xyz);
    return output;
}

[maxvertexcount(3)]
void GSmain(triangle VSOut input[3], inout TriangleStream<GSOut> output)
{
    float3 genN = normalize(cross(input[1].WorldPosition - input[0].WorldPosition, input[2].WorldPosition - input[0].WorldPosition));

    GSOut gs;
    float3 n0 = normalize(cross(genN, input[1].WorldPosition - input[0].WorldPosition));
    float3 n1 = normalize(cross(genN, input[2].WorldPosition - input[1].WorldPosition));
    float3 n2 = normalize(cross(genN, input[0].WorldPosition - input[2].WorldPosition));
    gs.EdgePlane0 = float4(n0, -dot(input[0].WorldPosition, n0));
    gs.EdgePlane1 = float4(n1, -dot(input[1].WorldPosition, n1));
    gs.EdgePlane2 = float4(n2, -dot(input[2].WorldPosition, n2));

    [unroll]
    for (int i = 0; i < 3; i++)
    {
        VSOut vs = input[i];
        if (SceneNode.AutoGenNormals.x != 0.0)
        {
            vs.WorldNormal = genN;
        }

        gs.vs = vs;
        output.Append(gs);
    }
    output.RestartStrip();
}

static const float kWireframeThickness = 0.001;

PSOut PSmain(GSOut gs, bool isFrontFace : SV_IsFrontFace)
{
    PSOut output;

    VSOut input = gs.vs;
    
    float4 diffuseMap;
    if (Material.HasDiffuse.x)
        diffuseMap = DiffuseTexture.Sample(DiffuseSampler, input.TexCoord);
    else
        diffuseMap = float4(1, 1, 1, 1);

    float specularMap;
    if (Material.HasSpecular.x)
        specularMap = SpecularTexture.Sample(SpecularSampler, input.TexCoord).r;
    else
        specularMap = 0.0;

    float3 N = normalize(input.WorldNormal);
    if (!isFrontFace)
        N = -N;
    float3 P = input.WorldPosition;
    float3 C = Camera.WorldPosition.xyz;

    if (Material.HasBump.x)
    {
        int3 boff = int3(-1, 0, 1);
        float b11 = BumpTexture.Sample(BumpSampler, input.TexCoord).r;
        float b01 = BumpTexture.Sample(BumpSampler, input.TexCoord, boff.xy).x;
        float b21 = BumpTexture.Sample(BumpSampler, input.TexCoord, boff.zy).x;
        float b10 = BumpTexture.Sample(BumpSampler, input.TexCoord, boff.yx).x;
        float b12 = BumpTexture.Sample(BumpSampler, input.TexCoord, boff.yz).x;
        float3 va = normalize(float3(2, 0, (b21 - b01) / input.Position.w * 3000));
        float3 vb = normalize(float3(0, 2, (b12 - b10) / input.Position.w * 3000));
        float4 bump = float4(cross(va, vb), b11);

        float3 worldTangent = normalize(input.WorldTangent.xyz) * input.WorldTangent.w;
        float3 worldBitangent = normalize(input.WorldBitangent);
        float3 worldNormal = normalize(input.WorldNormal);
        float3x3 tangentFrame = float3x3(worldTangent, worldBitangent, worldNormal);
        N = mul(transpose(tangentFrame), bump.xyz);
    }

    float3 V = normalize(C - P);
    float3 L = V;
    float G = max(0, dot(N, L));
    float3 R = reflect(-L, N);
    float S = pow(max(0, dot(R, V)), Material.Shininess.r);

    float4 ambient = diffuseMap * Material.Ambient;
    float4 diffuse = diffuseMap * Material.Diffuse * G;
    float4 specular = specularMap * Material.Specular * S;
    
    output.Color = float4(ambient.xyz + diffuse.xyz + specular.xyz, 1.0);
    output.Normal = N;

    float dist = min(dot(gs.EdgePlane0, float4(P, 1)), min(dot(gs.EdgePlane1, float4(P, 1)), dot(gs.EdgePlane2, float4(P, 1))));
    if (dist < kWireframeThickness)
    {
        output.Color *= float4((dist / kWireframeThickness).xxx, 0);
    }

    return output;
}