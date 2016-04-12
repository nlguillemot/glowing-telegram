#include "common.hlsl"

struct VSIn
{
    float4 Position : POSITION;
};

struct VSOut
{
    float4 Position : SV_Position;
};

struct GSOut
{
    float4 Position : SV_Position;
    float2 TexCoord : TEXCOORD;
};

struct PSOut
{
    float4 Color : SV_Target0;
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
    float3 worldPosition = mul(input.Position, SceneNode.WorldTransform).xyz;
    output.Position = mul(float4(worldPosition, 1), Camera.WorldViewProjection);
    return output;
}

[maxvertexcount(4)]
void GSmain(point VSOut input[1], inout TriangleStream<GSOut> output)
{
    float2 off2 = float2(0.01, 0.01);

    GSOut o;
    
    o.Position.zw = float2(input[0].Position.z / input[0].Position.w, 1.0);
    
    o.Position.xy = input[0].Position.xy / input[0].Position.w + float2(-off2.x, -off2.y);
    o.TexCoord = float2(0, 0);
    output.Append(o);

    o.Position.xy = input[0].Position.xy / input[0].Position.w + float2(-off2.x, off2.y);
    o.TexCoord = float2(0, 1);
    output.Append(o);
    
    o.Position.xy = input[0].Position.xy / input[0].Position.w + float2(off2.x, -off2.y);
    o.TexCoord = float2(1, 0);
    output.Append(o);

    o.Position.xy = input[0].Position.xy / input[0].Position.w + float2(off2.x, off2.y);
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

    PSOut output;
    output.Color = float4(0, c * c, 0, 1);
    return output;
}