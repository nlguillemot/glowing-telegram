#include "common.hlsl"

struct VSIn
{
    float4 Position : POSITION;
    uint VertexID : SV_VertexID;
};

struct VSOut
{
    float4 Position : POSITION;
    uint VertexID : VERTEXID;
};

struct GSOut
{
    float4 Position : SV_Position;
    float2 TexCoord : TEXCOORD;
    nointerpolation float3 WorldCenter : WORLDCENTER;
    nointerpolation uint VertexID : VERTEXID;
};

struct PSOut
{
    float4 Color : SV_Target0;
    float3 Normal : SV_Target1;
    uint VertexID : SV_Target2;
    float Depth : SV_Depth;
};

cbuffer CameraBuffer : register(BUFFER_REGISTER(SELECTOR_CAMERA_BUFFER_SLOT))
{
    PerCameraData Camera;
};

cbuffer SceneNodeBuffer : register(BUFFER_REGISTER(SELECTOR_SCENENODE_BUFFER_SLOT))
{
    PerSceneNodeData SceneNode;
};

VSOut VSmain(VSIn input)
{
    VSOut output;
    output.Position = input.Position;
    output.VertexID = input.VertexID;
    return output;
}

static const float WorldRadius = 0.004;

[maxvertexcount(4)]
void GSmain(point VSOut input[1], inout TriangleStream<GSOut> output)
{
    float3 worldPosition = mul(input[0].Position, SceneNode.WorldTransform).xyz;

    float3 across = normalize(cross(Camera.LookDirection.xyz, Camera.Up.xyz));
    float3 upward = normalize(cross(across, Camera.LookDirection.xyz));

    GSOut o;

    o.WorldCenter = worldPosition;
    o.VertexID = input[0].VertexID;

    o.Position = mul(float4(worldPosition + WorldRadius * (-across + -upward), 1.0), Camera.WorldViewProjection);
    o.TexCoord = float2(0, 0);
    output.Append(o);

    o.Position = mul(float4(worldPosition + WorldRadius * (-across + +upward), 1.0), Camera.WorldViewProjection);
    o.TexCoord = float2(0, 1);
    output.Append(o);
    
    o.Position = mul(float4(worldPosition + WorldRadius * (+across + -upward), 1.0), Camera.WorldViewProjection);
    o.TexCoord = float2(1, 0);
    output.Append(o);

    o.Position = mul(float4(worldPosition + WorldRadius * (+across + +upward), 1.0), Camera.WorldViewProjection);
    o.TexCoord = float2(1, 1);
    output.Append(o);

    output.RestartStrip();
}

PSOut PSmain(GSOut input)
{
    float c = length(input.TexCoord - float2(0.5, 0.5)) / 0.5;
    if (c > 1.0)
    {
        discard;
    }

    c = 1.0 - c * c;

    PSOut output;
    
    float3 across = normalize(cross(Camera.LookDirection.xyz, Camera.Up.xyz));
    float3 upward = normalize(cross(across, Camera.LookDirection.xyz));

    float3 pX = across * (input.TexCoord.x - 0.5) / 0.5;
    float3 pY = upward * (input.TexCoord.y - 0.5) / 0.5;
    float3 pZ = normalize(-Camera.LookDirection.xyz) * c;
    float3 toPos = WorldRadius * (pX + pY + pZ);
    float3 worldPos = input.WorldCenter + toPos;
    float4 clipPos = mul(float4(worldPos, 1.0), Camera.WorldViewProjection);
    output.Normal = normalize(toPos);
    output.Depth = clipPos.z / clipPos.w;
    output.Color = float4(0, dot(output.Normal, normalize(-Camera.LookDirection.xyz)), 0, 0.5);
    output.Color.rgb *= output.Color.a;
    output.VertexID = input.VertexID;
    return output;
}