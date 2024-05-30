#ifndef MODULES_OPENDT_LIB
#define MODULES_OPENDT_LIB


const int EOTF_lin     = 0; // Linear
const int EOTF_srgb    = 1; // 2.2 Power - sRGB Display
const int EOTF_rec1886 = 2; // 2.4 Power - Rec .1886
const int EOTF_dci     = 3; // 2.6 Power - DCI
const int EOTF_pq      = 4; // ST 2084 PQ
const int EOTF_hlg     = 5; // HLG

const mat3 matrix_rec2020_to_xyz = mat3(
    vec3(0.636958122253f, 0.144616916776f, 0.168880969286f),
    vec3(0.262700229883f, 0.677998125553f, 0.059301715344f),
    vec3(0.000000000000f, 0.028072696179, 1.060985088348f)
);

const mat3 in_to_xyz = matrix_rec2020_to_xyz;
mat3 xyz_to_display = inverse(matrix_rec2020_to_xyz);
mat3 xyz_to_in = xyz_to_display;
const mat3 display_to_xyz = matrix_rec2020_to_xyz;

/* Math helper functions ----------------------------*/

// Multiply 3x3 matrix m and vec3 vector v
vec3 vdot(mat3 m, vec3 v) {
    return v * m;
}

// Safe division of float a by float b
float sdivf(float a, float b) {
    if (b == 0.0) return 0.0;
    else return a / b;
}

// Safe division of vec3 a by float b
vec3 sdivf3f(vec3 a, float b) {
    return vec3(sdivf(a.x, b), sdivf(a.y, b), sdivf(a.z, b));
}

// Safe element-wise division of vec3 a by vec3 b
vec3 sdivf3f3(vec3 a, vec3 b) {
    return vec3(sdivf(a.x, b.x), sdivf(a.y, b.y), sdivf(a.z, b.z));
}

// Safe power function raising float a to power float b
#define spowf(a, b) (((a) <= 0.0) ? (a) : pow(a, b))

// Safe power function raising vec3 a to power float b
vec3 spowf3(vec3 a, float b) { return vec3(spowf(a.x, b), spowf(a.y, b), spowf(a.z, b)); }

#define maxf3 max_of
#define minf3 min_of

// Clamp vec3 a to max value mx
vec3 clampmaxf3(vec3 a, float mx) { return vec3(min(a.x, mx), min(a.y, mx), min(a.z, mx)); }

// Clamp vec3 a to min value mn
vec3 clampminf3(vec3 a, float mn) { return vec3(max(a.x, mn), max(a.y, mn), max(a.z, mn)); }

// Clamp each component of vec3 a to be between float mn and float mx
vec3 clampf3(vec3 a, float mn, float mx) { 
    return vec3(min(max(a.x, mn), mx), min(max(a.y, mn), mx), min(max(a.z, mn), mx));
}



vec3 eotf_hlg(vec3 rgb, int inverse) {
  /* Apply the HLG Forward or Inverse EOTF. Implements the full ambient surround illumination model
      ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
      ITU-R Rep BT.2390-8: https://www.itu.int/pub/R-REP-BT.2390
      Perceptual Quantiser (PQ) to Hybrid Log-Gamma (HLG) Transcoding: 
        https://www.bbc.co.uk/rd/sites/50335ff370b5c262af000004/assets/592eea8006d63e5e5200f90d/BBC_HDRTV_PQ_HLG_Transcode_v2.pdf
  */

    const float HLG_Lw = 1000.0;
    // const float HLG_Lb = 0.0;
    const float HLG_Ls = 5.0;
    const float h_a = 0.17883277;
    const float h_b = 1.0 - 4.0 * 0.17883277;
    const float h_c = 0.5 - h_a * log(4.0 * h_a);
    const float h_g = 1.2 * pow(1.111, log2(HLG_Lw / 1000.0)) * pow(0.98, log2(max(1e-6, HLG_Ls) / 5.0));
    if (inverse == 1) {
        float Yd = dot(rgb, vec3(0.2627, 0.6780, 0.0593));
        // HLG Inverse OOTF
        rgb = rgb * pow(Yd, (1.0 - h_g) / h_g);
        // HLG OETF
        rgb.x = rgb.x <= 1.0 / 12.0 ? sqrt(3.0 * rgb.x) : h_a * log(12.0 * rgb.x - h_b) + h_c;
        rgb.y = rgb.y <= 1.0 / 12.0 ? sqrt(3.0 * rgb.y) : h_a * log(12.0 * rgb.y - h_b) + h_c;
        rgb.z = rgb.z <= 1.0 / 12.0 ? sqrt(3.0 * rgb.z) : h_a * log(12.0 * rgb.z - h_b) + h_c;
    } else {
        // HLG Inverse OETF
        rgb.x = rgb.x <= 0.5 ? rgb.x * rgb.x / 3.0 : (exp((rgb.x - h_c) / h_a) + h_b) / 12.0;
        rgb.y = rgb.y <= 0.5 ? rgb.y * rgb.y / 3.0 : (exp((rgb.y - h_c) / h_a) + h_b) / 12.0;
        rgb.z = rgb.z <= 0.5 ? rgb.z * rgb.z / 3.0 : (exp((rgb.z - h_c) / h_a) + h_b) / 12.0;
        // HLG OOTF
        float Ys = dot(rgb, vec3(0.2627, 0.6780, 0.0593));
        rgb = rgb * pow(Ys, h_g - 1.0);
    }
    return rgb;
}

vec3 eotf_pq(vec3 rgb, int inverse) {
  /* Apply the ST-2084 PQ Forward or Inverse EOTF
      ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
      ITU-R Rep BT.2390-9 https://www.itu.int/pub/R-REP-BT.2390
      Note: in the spec there is a normalization for peak display luminance. 
      For this function we assume the input is already normalized such that 1.0 = 10,000 nits
  */

    // const float Lp = 1.0;
    const float m1 = 2610.0 / 16384.0;
    const float m2 = 2523.0 / 32.0;
    const float c1 = 107.0 / 128.0;
    const float c2 = 2413.0 / 128.0;
    const float c3 = 2392.0 / 128.0;

    if (inverse == 1) {
        // rgb /= Lp;
        rgb = spowf3(rgb, m1);
        rgb = spowf3((c1 + c2 * rgb) / (1.0 + c3 * rgb), m2);
    } else {
        rgb = spowf3(rgb, 1.0 / m2);
        rgb = spowf3((rgb - c1) / (c2 - c3 * rgb), 1.0 / m1);
        // rgb *= Lp;
    }
    return rgb;
}

#endif // MODULES_OPENDT_LIB