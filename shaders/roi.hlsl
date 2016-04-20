#include "common.hlsl"

struct VSIn
{
    float2 Position : POSITION;
};

struct VSOut
{
    float2 Position : POSITION;
};

struct GSOut
{
    float4 Position : SV_Position;
    float2 TexCoord : TEXCOORD;
    float LineLength : LINELENGTH;
};

struct PSOut
{
    float4 Color : SV_Target0;
};

cbuffer ViewportBuffer : register(BUFFER_REGISTER(ROI_VIEWPORT_BUFFER_SLOT))
{
    PerViewportData Viewport;
};

static const float LineRadius = 4.0;

VSOut VSmain(VSIn input)
{
    VSOut output;
    output.Position = input.Position;
    return output;
}

[maxvertexcount(4)]
void GSmain(line VSOut input[2], inout TriangleStream<GSOut> output)
{
    GSOut o;

    // screen coordinates with origin at bottom left
    float2 screen0 = input[0].Position * Viewport.Size.xy;
    screen0.y = Viewport.Size.y - screen0.y;
    float2 screen1 = input[1].Position * Viewport.Size.xy;
    screen1.y = Viewport.Size.y - screen1.y;

    float2 delta = screen1 - screen0;
    
    float2 n = normalize(float2(-delta.y, delta.x));
    float r = LineRadius;

    o.LineLength = length(screen1 - screen0);

    screen0 = screen0 + -normalize(delta) * LineRadius;
    screen1 = screen1 + normalize(delta) * LineRadius;

    o.Position = float4((screen0 + n * r) / Viewport.Size.xy * float2(2.0, 2.0) - float2(1.0, 1.0), 0.0, 1.0);
    o.TexCoord = float2(0, 1);
    output.Append(o);

    o.Position = float4((screen1 + n * r) / Viewport.Size.xy * float2(2.0, 2.0) - float2(1.0, 1.0), 0.0, 1.0);
    o.TexCoord = float2(0, 0);
    output.Append(o);

    o.Position = float4((screen0 - n * r) / Viewport.Size.xy * float2(2.0, 2.0) - float2(1.0, 1.0), 0.0, 1.0);
    o.TexCoord = float2(1, 1);
    output.Append(o);

    o.Position = float4((screen1 - n * r) / Viewport.Size.xy * float2(2.0, 2.0) - float2(1.0, 1.0), 0.0, 1.0);
    o.TexCoord = float2(1, 0);
    output.Append(o);

    output.RestartStrip();
}

PSOut PSmain(GSOut input)
{
    PSOut output;
    
    float horizDist = (1.0 - abs(input.TexCoord.x - 0.5) / 0.5);

    output.Color = float4(0.0, 0.5, 0.3, horizDist);

    float distFromMiddle = abs(input.TexCoord.y - 0.5) * input.LineLength;
    float capDistanceFromMiddle = input.LineLength / 2.0 - LineRadius;
    float capRadiusDist = length(float2((input.TexCoord.x - 0.5) / 0.5 * LineRadius, distFromMiddle - capDistanceFromMiddle)) / LineRadius;

    float distFromLine = horizDist;
    if (distFromMiddle > capDistanceFromMiddle)
    {
        distFromLine = 1.0 - capRadiusDist;
    }
    output.Color.a = min(output.Color.a, distFromLine);

    output.Color.r = min(1.0, max(0.0, 1.0 - capRadiusDist));
    output.Color.rgb *= distFromLine * distFromLine * distFromLine * distFromLine; // add a dark outline

    output.Color.rgb *= output.Color.a; // premultiplied alpha
    return output;
}