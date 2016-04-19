#include "common.hlsl"

struct VSIn
{
    uint VertexID : SV_VertexID;
};

struct VSOut
{
    float4 Position : SV_Position;
    float2 TexCoord : TEXCOORD;
};

struct PSOut
{
    float3 Occlusion : SV_Target0;
};

Texture2D NormalTexture : register(TEXTURE_REGISTER(SSAO_NORMAL_TEXTURE_SLOT));
Texture2D DepthTexture : register(TEXTURE_REGISTER(SSAO_DEPTH_TEXTURE_SLOT));
Texture2D NoiseTexture : register(TEXTURE_REGISTER(SSAO_NOISE_TEXTURE_SLOT));

SamplerState NormalSampler : register(SAMPLER_REGISTER(SSAO_NORMAL_SAMPLER_SLOT));
SamplerState DepthSampler : register(SAMPLER_REGISTER(SSAO_DEPTH_SAMPLER_SLOT));
SamplerState NoiseSampler : register(SAMPLER_REGISTER(SSAO_NOISE_SAMPLER_SLOT));

// SSAO parameters
static const float SSAO_Strength = 0.25;
static const float2 SSAO_Offset = float2(1280.0, 720.0) / 4.0;
static const float SSAO_Falloff = 0.0000000002;
static const float SSAO_Radius = 0.4;

#define NUM_SAMPLES 10
static const float SSAO_InvSamples = 1.0 / (float)NUM_SAMPLES;

// AO sampling directions 
static const float3 AO_SAMPLES[NUM_SAMPLES] =
{
    float3(-0.010735935,  0.0164701800,  0.0062425877),
    float3(-0.065333690,  0.3647007000, -0.1374632100),
    float3(-0.653923500, -0.0167263880, -0.5300095700),
    float3(0.409582850,  0.0052428036, -0.5591124000),
    float3(-0.146536600,  0.0989926700,  0.1557167900),
    float3(-0.441221120, -0.5458797000,  0.0491253200),
    float3(0.037555660, -0.1096134500, -0.3304027300),
    float3(0.019100213,  0.2965278300,  0.0662376660),
    float3(0.876532300,  0.0112360040,  0.2826596200),
    float3(0.292644350, -0.4079423800,  0.1596416700)
};

// http://www.slideshare.net/DevCentralAMD/vertex-shader-tricks-bill-bilodeau
VSOut VSmain(VSIn input)
{
    VSOut output;
    
    output.Position.x = (float)(input.VertexID / 2) * 4.0 - 1.0;
    output.Position.y = (float)(input.VertexID % 2) * 4.0 - 1.0;
    output.Position.z = 0.0;
    output.Position.w = 1.0;

    output.TexCoord.x = (float)(input.VertexID / 2) * 2.0;
    output.TexCoord.y = 1.0 - (float)(input.VertexID % 2) * 2.0;

    return output;
}

// http://www.justinmclaine.com/shader---ssao.html
PSOut PSmain(VSOut input)
{
    PSOut output;

    float depth = DepthTexture.Sample(DepthSampler, input.TexCoord).x;
    if (depth == 1.0)
    {
        discard;
    }

    float3 normal = NormalTexture.Sample(NormalSampler, input.TexCoord).xyz;

    float3 fres = normalize(NoiseTexture.Sample(NoiseSampler, input.TexCoord * SSAO_Offset).xyz);

    float bl = 0.0;
    float radD = SSAO_Radius / depth;

    float3 ray;
    float3 sampleNormal;
    float  depthDiff;

    [unroll]
    for (int i = 0; i < NUM_SAMPLES; ++i)
    {
        ray = radD * reflect(AO_SAMPLES[i], fres);
        
        float2 sampleLoc = input.TexCoord.xy + sign(dot(ray, normal)) * ray.xy;

        sampleNormal = NormalTexture.Sample(NormalSampler, sampleLoc).xyz;
        depthDiff = depth - DepthTexture.Sample(DepthSampler, sampleLoc).x;
        
        bl += step(SSAO_Falloff, depthDiff) * (1.0 - dot(sampleNormal, normal)) *
            (1.0 - smoothstep(SSAO_Falloff, SSAO_Strength, depthDiff));
    }

    float ao = bl * SSAO_InvSamples;
    output.Occlusion = float3(ao, ao, ao);
    return output;
}