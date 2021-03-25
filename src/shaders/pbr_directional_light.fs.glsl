#version 330

in vec3 vViewSpacePosition;
in vec3 vViewSpaceNormal;
in vec2 vTexCoords;

uniform vec3 uLightDirection;
uniform vec3 uLightIntensity;

uniform sampler2D uBaseColorTexture;
uniform vec4 uBaseColorFactor;

uniform sampler2D uEmissiveTexture;
uniform vec3 uEmissiveFactor;

uniform sampler2D uOcclusionTexture;
uniform float uOcclusionFactor;
uniform int uOcclusionEnabled;

uniform float uMetallicFactor;
uniform float uRoughnessFactor;
uniform sampler2D uMetallicRoughnessTexture;

out vec3 fColor;

// Constants
const float GAMMA = 2.2;
const float INV_GAMMA = 1. / GAMMA;
const float M_PI = 3.141592653589793;
const float M_1_PI = 1.0 / M_PI;

// We need some simple tone mapping functions
// Basic gamma = 2.2 implementation
// Stolen here:
// https://github.com/KhronosGroup/glTF-Sample-Viewer/blob/master/src/shaders/tonemapping.glsl

// linear to sRGB approximation
// see http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
vec3 LINEARtoSRGB(vec3 color) { return pow(color, vec3(INV_GAMMA)); }

// sRGB to linear approximation
// see http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
vec4 SRGBtoLINEAR(vec4 srgbIn)
{
  return vec4(pow(srgbIn.xyz, vec3(GAMMA)), srgbIn.w);
}

int heaviside(float x) {
    return (x > 0) ? 1 : 0;
}

// https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#metal-brdf-and-dielectric-brdf
void main()
{
  vec3 N = normalize(vViewSpaceNormal);
  vec3 L = uLightDirection;
  vec3 V = normalize(-vViewSpacePosition);
  vec3 H = normalize(L + V);

  // Read color texture.
  vec4 baseColorFromTexture =
      SRGBtoLINEAR(texture(uBaseColorTexture, vTexCoords));
  vec4 baseColor = baseColorFromTexture * uBaseColorFactor;
  float NdotL = clamp(dot(N, L), 0., 1.);

  // Read metallic/roughness texture.
  vec4 metallicRoughnessFromTexture = texture(uMetallicRoughnessTexture, vTexCoords);
  float metallic = metallicRoughnessFromTexture.b * uMetallicFactor;
  float roughness = metallicRoughnessFromTexture.g * uRoughnessFactor;

  // Some 'constants'.
  float dielectricSpecular = 0.04;
  vec3 black = vec3(0.);
  float alpha = roughness * roughness;


  // Compute specular BRDF.
  float D_numerator = alpha * alpha * heaviside(dot(N, H));
  float D_denom = M_PI * pow((dot(N, H) * dot(N, H) * (alpha * alpha - 1) + 1), 2);
  float D = D_numerator / D_denom;

  float V_spec = 0;
  float V_numerator_1 = heaviside(dot(H, L));
  float V_numerator_2 = heaviside(dot(H, V));
  float V_denom_1 = abs(dot(N, L)) + sqrt(alpha * alpha + (1 - alpha * alpha) * dot(N, L) * dot(N, L));
  float V_denom_2 = abs(dot(N, V)) + sqrt(alpha * alpha + (1 - alpha * alpha) * dot(N, V) * dot(N, V));
  if(V_denom_1 > 0 && V_denom_2 > 0) {
      float V_1 = V_numerator_1 / V_denom_1;
      float V_2 = V_numerator_2 / V_denom_2;
      V_spec = V_1 * V_2;
  }

  float specular_brdf = V_spec * D;


  // Mix BRDFs computation (see the link above).
  vec3 c_diff = mix(baseColor.rgb * (1 - dielectricSpecular), black, metallic);
  vec3 f0 = mix(vec3(dielectricSpecular), baseColor.rgb, metallic);

  float baseFresnelFactor = 1 - abs(dot(V, H));
  float fresnelFactor = baseFresnelFactor * baseFresnelFactor; // power 2
  fresnelFactor *= fresnelFactor; // power 4
  fresnelFactor *= baseFresnelFactor; // power 5
  vec3 F = f0 + ((1 - f0) * fresnelFactor);

  vec3 f_diffuse = (1 - F) * (M_1_PI) * c_diff;
  vec3 f_specular = F * specular_brdf;

  vec3 material = f_diffuse + f_specular;

  // Emissive.
  vec4 emissiveFromTexture =
      SRGBtoLINEAR(texture(uEmissiveTexture, vTexCoords));
  vec3 emissive = emissiveFromTexture.rgb * uEmissiveFactor;


  vec3 color = material * uLightIntensity * NdotL + emissive;

  // Occlusion (based on the provided solution).
  if(uOcclusionEnabled != 0) {
    float ao = texture2D(uOcclusionTexture, vTexCoords).r;
    color = mix(color, color * ao, uOcclusionFactor);
  }

  // Output.
  fColor = LINEARtoSRGB(color);
}