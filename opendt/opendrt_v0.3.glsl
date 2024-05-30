#ifndef MODULES_OPENDT_OPENDRT_v03
#define MODULES_OPENDT_OPENDRT_v03



/*  OpenDRT -------------------------------------------------/
      v0.3.2
      Written by Jed Smith
      https://github.com/jedypod/open-display-transform

      License: GPL v3
-------------------------------------------------*/

//
// Converted to GLSL by sw-52
//

#include "/modules/opendt/lib.glsl"



/*#define OPENDRT_03_LP     100.0
#define OPENDRT_03_LG      10.0
#define OPENDRT_03_LG_BOOST 0.12
#define OPENDRT_03_CONTRAST 1.40
#define OPENDRT_03_TOE      0.001

#define OPENDRT_03_PC_P  0.30
#define OPENDRT_03_PB    0.50
#define OPENDRT_03_HS_R  0.30
#define OPENDRT_03_HS_G -0.10
#define OPENDRT_03_HS_B -0.30

#define OPENDRT_03_SAT_F   0.40
#define OPENDRT_03_SAT_W_R 0.15
#define OPENDRT_03_SAT_W_G 0.50
#define OPENDRT_03_SAT_W_B 0.35
#define OPENDRT_03_DN_W_C  0.70
#define OPENDRT_03_DN_W_M  0.60
#define OPENDRT_03_DN_W_Y  0.80*/

/*
// Tonescale Parameters
DEFINE_UI_PARAMS(Lp, Lp, DCTLUI_SLIDER_FLOAT, 100.0, 100.0, 1000.0, 0.0)
DEFINE_UI_PARAMS(Lg, Lg, DCTLUI_SLIDER_FLOAT, 10.0, 3.0, 30.0, 0.0)
DEFINE_UI_PARAMS(Lgb, Lg boost, DCTLUI_SLIDER_FLOAT, 0.12, 0.0, 0.5, 0.0)
DEFINE_UI_PARAMS(p, contrast, DCTLUI_SLIDER_FLOAT, 1.4, 1.0, 2.0, 0.0)
DEFINE_UI_PARAMS(toe, toe, DCTLUI_SLIDER_FLOAT, 0.001, 0.0, 0.02, 0.0)

// Color Parameters
DEFINE_UI_PARAMS(pc_p, purity compress, DCTLUI_SLIDER_FLOAT, 0.3, 0.0, 1.0, 0.0)
DEFINE_UI_PARAMS(pb, purity boost, DCTLUI_SLIDER_FLOAT, 0.3, 0.0, 1.0, 0.0)
DEFINE_UI_PARAMS(hs_r, hueshift r, DCTLUI_SLIDER_FLOAT, 0.3, -1.0, 1.0, 0.0)
DEFINE_UI_PARAMS(hs_g, hueshift g, DCTLUI_SLIDER_FLOAT, 0.0, -1.0, 1.0, 0.0)
DEFINE_UI_PARAMS(hs_b, hueshift b, DCTLUI_SLIDER_FLOAT, -0.3, -1.0, 1.0, 0.0)

// Encoding / IO
DEFINE_UI_PARAMS(in_gamut, in gamut, DCTLUI_COMBO_BOX, 14, {i_xyz, i_ap0, i_ap1, i_p3d65, i_rec2020, i_rec709, i_awg3, i_awg4, i_rwg, i_sgamut3, i_sgamut3cine, i_vgamut, i_bmdwg, i_egamut, i_davinciwg}, {XYZ, ACES 2065-1, ACEScg, P3D65, Rec.2020, Rec.709, Arri Wide Gamut 3, Arri Wide Gamut 4, Red Wide Gamut RGB, Sony SGamut3, Sony SGamut3Cine, Panasonic V-Gamut, Blackmagic Wide Gamut, Filmlight E-Gamut, DaVinci Wide Gamut})
DEFINE_UI_PARAMS(in_oetf, in transfer function, DCTLUI_COMBO_BOX, 0, {ioetf_linear, ioetf_davinci_intermediate, ioetf_filmlight_tlog, ioetf_arri_logc3, ioetf_arri_logc4, ioetf_panasonic_vlog, ioetf_sony_slog3, ioetf_fuji_flog}, {Linear, Davinci Intermediate, Filmlight T-Log, Arri LogC3, Arri LogC4, Panasonic V-Log, Sony S-Log3, Fuji F-Log})
DEFINE_UI_PARAMS(display_gamut, display gamut, DCTLUI_COMBO_BOX, 0, {Rec709, P3D65, Rec2020}, {Rec.709, P3 D65, Rec.2020})
DEFINE_UI_PARAMS(EOTF, display eotf, DCTLUI_COMBO_BOX, 2, {lin, srgb, rec1886, dci, pq, hlg}, {Linear, 2.2 Power sRGB Display, 2.4 Power Rec .1886, 2.6 Power DCI, ST 2084 PQ, HLG})
*/

// hypotf3 = length(vec3 x)

// Input gamut conversion matrix (CAT02 chromatic adaptation to D65)
//#define in_to_xyz matrix_rec2020_to_xyz
//#define xyz_to_display matrix_xyz_to_rec2020
#define eotf EOTF_lin

/*
// Gamut Conversion Matrices
const mat3 matrix_xyz_to_xyz = mat3(1.0);
const mat3 matrix_ap0_to_xyz = mat3(vec3(0.93863094875f, -0.00574192055f, 0.017566898852f), vec3(0.338093594922f, 0.727213902811f, -0.065307497733f), vec3(0.000723121511f, 0.000818441849f, 1.0875161874f));
const mat3 matrix_ap1_to_xyz = mat3(vec3(0.652418717672f, 0.127179925538f, 0.170857283842f), vec3(0.268064059194f, 0.672464478993f, 0.059471461813f), vec3(-0.00546992851f, 0.005182799977f, 1.08934487929f));
const mat3 matrix_rec709_to_xyz = mat3(vec3(0.412390917540f, 0.357584357262f, 0.180480793118f), vec3(0.212639078498f, 0.715168714523f, 0.072192311287f), vec3(0.019330825657f, 0.119194783270f, 0.950532138348f));
const mat3 matrix_p3d65_to_xyz = mat3(vec3(0.486571133137f, 0.265667706728f, 0.198217317462f), vec3(0.228974640369f, 0.691738605499f, 0.079286918044f), vec3(-0.000000000000f, 0.045113388449, 1.043944478035f));
const mat3 matrix_rec2020_to_xyz = mat3(vec3(0.636958122253f, 0.144616916776f, 0.168880969286f), vec3(0.262700229883f, 0.677998125553f, 0.059301715344f), vec3(0.000000000000f, 0.028072696179, 1.060985088348f));
const mat3 matrix_arriwg3_to_xyz = mat3(vec3(0.638007619284f, 0.214703856337f, 0.097744451431f), vec3(0.291953779f, 0.823841041511f, -0.11579482051f), vec3(0.002798279032f, -0.067034235689f, 1.15329370742f));
const mat3 matrix_arriwg4_to_xyz = mat3(vec3(0.704858320407f, 0.12976029517f, 0.115837311474f), vec3(0.254524176404f, 0.781477732712f, -0.036001909116f), vec3(0.0f, 0.0f, 1.08905775076f));
const mat3 matrix_redwg_to_xyz = mat3(vec3(0.735275208950f, 0.068609409034f, 0.146571278572f), vec3(0.286694079638f, 0.842979073524f, -0.129673242569f), vec3(-0.079680845141f, -0.347343206406, 1.516081929207f));
const mat3 matrix_sonysgamut3_to_xyz = mat3(vec3(0.706482713192f, 0.128801049791f, 0.115172164069f), vec3(0.270979670813f, 0.786606411221f, -0.057586082034f), vec3(-0.009677845386f, 0.004600037493f, 1.09413555865f));
const mat3 matrix_sonysgamut3cine_to_xyz = mat3(vec3(0.599083920758f, 0.248925516115f, 0.102446490178f), vec3(0.215075820116f, 0.885068501744f, -0.100144321859f), vec3(-0.032065849545f, -0.027658390679f, 1.14878199098f));
const mat3 matrix_bmdwg_to_xyz = mat3(vec3(0.606538414955f, 0.220412746072f, 0.123504832387f), vec3(0.267992943525f, 0.832748472691f, -0.100741356611f), vec3(-0.029442556202f, -0.086612440646, 1.205112814903f));
const mat3 matrix_egamut_to_xyz = mat3(vec3(0.705396831036f, 0.164041340351f, 0.081017754972f), vec3(0.280130714178f, 0.820206701756f, -0.100337378681f), vec3(-0.103781513870f, -0.072907261550, 1.265746593475f));
const mat3 matrix_davinciwg_to_xyz = mat3(vec3(0.700622320175f, 0.148774802685f, 0.101058728993f), vec3(0.274118483067f, 0.873631775379f, -0.147750422359f), vec3(-0.098962903023f, -0.137895315886, 1.325916051865f));

const mat3 matrix_xyz_to_rec709 = mat3(vec3(3.2409699419f, -1.53738317757f, -0.498610760293f), vec3(-0.969243636281f, 1.87596750151f, 0.041555057407f), vec3(0.055630079697f, -0.203976958889f, 1.05697151424f));
const mat3 matrix_xyz_to_p3d65 = mat3(vec3(2.49349691194f, -0.931383617919f, -0.402710784451f), vec3(-0.829488969562f, 1.76266406032f, 0.023624685842f), vec3(0.035845830244f, -0.076172389268f, 0.956884524008f));
const mat3 matrix_xyz_to_rec2020 = mat3(vec3(1.71665118797f, -0.355670783776f, -0.253366281374f), vec3(-0.666684351832f, 1.61648123664f, 0.015768545814f), vec3(0.017639857445f, -0.042770613258f, 0.942103121235f));
*/



/* Functions for the OpenDRT Transform ---------------------------------------- */

float compress_powerptoe(float x, const float p, const float x0, const float t0, int inv) {
    /* Variable slope compression function.
        p: Slope of the compression curve. Controls how compressed values are distributed. 
            p=0.0 is a clip. p=1.0 is a hyperbolic curve.
        x0: Compression amount. How far to reach outside of the gamut boundary to pull values in.
        t0: Threshold point within gamut to start compression. t0=0.0 is a clip.
        https://www.desmos.com/calculator/igy3az7maq
    */
    // Precalculations for Purity Compress intersection constraint at (-x0, 0)
    const float m0 = spowf((t0 + max(1e-6, x0)) / t0, 1.0 / p) - 1.0;
    const float m = spowf(m0, -p) * (t0 * spowf(m0, p) - t0 - max(1e-6, x0));

  float i = inv == 1 ? -1.0 : 1.0;
  return x > t0 ? x : (x - t0) * spowf(1.0 + i * spowf((t0 - x) / (t0 - m), 1.0 / p), -p) + t0;
}

float hyperbolic_compress(float x, float m, float s, float p, int inv) {
    if (inv == 0) {
        return spowf(m * x / (x + s), p);
    } else {
        float ip = 1.0 / p;
        return spowf(s * x, ip) / (m - spowf(x, ip));
    }
}

#define quadratic_toe_compress(x, toe, inv) ((toe) == 0.0 ? (x) : ((inv) == 0 ? spowf((x), 2.0) / ((x) + (toe)) : ((x) + sqrt((x) * (4.0 * (toe) + (x)))) / 2.0))

/*float quadratic_toe_compress(float x, float toe, int inv) {
    if (toe == 0.0) return x;
    if (inv == 0) {
        return spowf(x, 2.0) / (x + toe);
    } else {
        return (x + sqrt(x * (4.0 * toe + x))) / 2.0;
    }
}*/

float tonescale(float x, const float Lp, const float Lg, const float Lgb, const float p, const float toe, int inv) {
    // input scene-linear peak x intercept
    const float px = 256.0 * log(Lp) / log(100.0) - 128.0;
    // output display-linear peak y intercept
    const float py = Lp / 100.0;
    // input scene-linear middle grey x intercept
    const float gx = 0.18;
    // output display-linear middle grey y intercept
    const float gy = Lg / 100.0 * (1.0 + Lgb * log(py) / log(2.0));
    // s0 and s are input x scale for middle grey intersection constraint
    // m0 and m are output y scale for peak white intersection constraint
    const float s0 = quadratic_toe_compress(gy, toe, 1);
    const float m0 = quadratic_toe_compress(py, toe, 1);
    const float ip = 1.0 / p;
    const float s1 = pow(s0, ip);
    const float m1 = pow(m0, ip);
    const float s = (px * gx * (m1 - s1)) / (px * s1 - gx * m1);
    const float m = m1 * (s + px) / px;

    if (inv == 0) {
        x = hyperbolic_compress(x, m, s, p, 0);
        return quadratic_toe_compress(x, toe, 0) / py;
    } else {
        x = quadratic_toe_compress(x * py, toe, 1);
        return hyperbolic_compress(x, m, s, p, 1);
    }
}

/*// https://www.desmos.com/calculator/gfubm2kvlu
float powerp(float x, float p, float m) {
    float y = x <= 0.0f ? x : x*spowf(spowf(x/m, 1.0/p) + 1.0, -p);
    return y;
}

// https://www.desmos.com/calculator/jrff9lrztn
float powerptoe(float x, float p, float m, float t0) {
    float y = x > t0 ? x : (x-t0)*spowf(spowf((t0-x)/(t0-m), 1.0/p) + 1.0, -p) + t0;
    return y;
}*/

vec3 opendrtransform_v03(vec3 rgb) {
    // **************************************************
    // Parameter Setup
    // --------------------------------------------------

    const float Lp  = OPENDRT_03_LP;       // Max Output        100.0 [100.0 - 1000.0]
    const float Lg  = OPENDRT_03_LG;       // Middle Gray       10.0  [3.0 - 30.0]
    const float Lgb = OPENDRT_03_LG_BOOST; // Middle Gray Boost 0.12  [0.0 - 0.5]
    const float p   = OPENDRT_03_CONTRAST; // Contrast          1.4   [1.0 - 2.0]
    const float toe = OPENDRT_03_TOE;      // Toe               0.001 [0.0 - 0.02]

    const float pc_p = OPENDRT_03_PC_P; // Purity Compress  0.3 [0.0 - 1.0]
    const float pb   = OPENDRT_03_PB;   // Purity Boost     0.5 [0.0 - 1.0]
    const float hs_r = OPENDRT_03_HS_R; // Hue Shift Red     0.3 [-1.0 - 1.0]
    const float hs_g = OPENDRT_03_HS_G; // Hue Shift Green  -0.1 [-1.0 - 1.0]
    const float hs_b = OPENDRT_03_HS_B; // Hue Shift Blue   -0.3 [-1.0 - 1.0]

    // Hue Shift RGB controls
    const vec3 hs = vec3(hs_r, hs_g, hs_b);

    /* Display Scale ---------------*
          Remap peak white in display linear depending on the selected inverse EOTF.
          In our tonescale model, 1.0 is 100 nits, and as we scale up peak display luminance (Lp),
          we multiply up by the same amount. So if Lp=1,000, peak output of the tonescale model
          will be 10.0.

          So in ST2084 PQ, 1.0 is 10,000 nits, so we need to divide by 100 to fit out output into the
          container.

          Similarly in HLG, 1.0 is 1,000 nits, so we need to divide by 10.

          If we are in an SDR mode, instead we just scale the peak so it hits display 1.0.
      */
    // const float ds = eotf == 4 ? 0.01f : eotf == 5 ? 0.1f : 100.0f/Lp;
    const float ds = eotf == 4 ? Lp / 10000.0 : (eotf == 5 ? Lp / 1000.0 : 1.0f);



    /* Parameters which _could_ be tweaked but are not exposed 
        ------------------------------------------ */
    // "Saturation" amount
    const float sat_f = OPENDRT_03_SAT_F; // 0.4
    // "Saturation" weights
    const vec3 sat_w = vec3(OPENDRT_03_SAT_W_R, OPENDRT_03_SAT_W_G, OPENDRT_03_SAT_W_B); // 0.15 0.5 0.35
    // Density weights CMY
    const vec3 dn_w = vec3(OPENDRT_03_DN_W_C, OPENDRT_03_DN_W_M, OPENDRT_03_DN_W_Y); // 0.7 0.6 0.8


    /* Rendering Code ------------------------------------------ */

    // Convert into display gamut
    /*
    rgb = rgb * in_to_xyz;
    rgb = rgb * xyz_to_display;
    */

    //  "Desaturate" to control shape of color volume in the norm ratios (Desaturate in scare quotes because the weights are creative)
    float sat_L = dot(rgb, sat_w);
    rgb = sat_L * (1.0 - sat_f) + rgb * sat_f;


    // Norm and RGB Ratios
    float norm = length(clampminf3(rgb, 0.0)) / sqrt(3.0);
    rgb = sdivf3f(rgb, norm);
    rgb = clampminf3(rgb, -2.0); // Prevent bright pixels from crazy values in shadow grain


    /* Purity Compression --------------------------------------- */
    // rgb purity compression strength, tuned for common camera gamuts
    /*vec3 pc_rats = rgb; TODO: Remove

    // minrgb with out of gamut values compressed
    pc_rats.x = compress_powerptoe(pc_rats.x, 0.25, 0.06, 1.0, 0);
    pc_rats.y = compress_powerptoe(pc_rats.y, 0.25, 0.2, 1.0, 0);
    pc_rats.z = compress_powerptoe(pc_rats.z, 0.25, 0.06, 1.0, 0);

    // minrgb with out of gamut values compressed
    float pc_rats_mn = minf3(pc_rats);

    // remove achromatic
    pc_rats = pc_rats - pc_rats_mn;
    // restore full range for tonescale curve
    pc_rats_mn = max(pc_rats_mn, 0.0) * norm;

    // the factor to mix back to maxrgb
    vec3 pc_s = vec3(0.04, 0.2, 0.08);
    float pc_f = dot(pc_s, pc_rats);

    // minrgb mix to maxrgb by pc_f
    float ccf = max(0.0, pc_rats_mn * (1.0 - pc_f) + norm * pc_f);*/


    /* Tonescale Parameters
          ----------------------
        For the tonescale compression function, we use one inspired by the wisdom shared by Daniele Siragusano
        on the tonescale thread on acescentral: https://community.acescentral.com/t/output-transform-tone-scale/3498/224

        This is a variation which puts the power function _after_ the display-linear scale, which allows a simpler and exact
        solution for the intersection constraints. The resulting function is pretty much identical to Daniele's but simpler.
        Here is a desmos graph with the math. https://www.desmos.com/calculator/hglnae2ame

        And for more info on the derivation, see the "Michaelis-Menten Constrained" Tonescale Function here:
        https://colab.research.google.com/drive/1aEjQDPlPveWPvhNoEfK4vGH5Tet8y1EB#scrollTo=Fb_8dwycyhlQ

        For the user parameter space, we include the following creative controls:
        - Lp: display peak luminance. This sets the display device peak luminance and allows rendering for HDR.
        - contrast: This is a pivoted power function applied after the hyperbolic compress function,
            which keeps middle grey and peak white the same but increases contrast in between.
        - flare: Applies a parabolic toe compression function after the hyperbolic compression function.
            This compresses values near zero without clipping. Used for flare or glare compensation.
        - gb: Grey Boost. This parameter controls how many stops to boost middle grey per stop of peak luminance increase.   // stops to boost Lg per stop of Lp increase

        Notes on the other non user-facing parameters:
        - (px, py): This is the peak luminance intersection constraint for the compression function.
            px is the input scene-linear x-intersection constraint. That is, the scene-linear input value
            which is mapped to py through the compression function. By default this is set to 128 at Lp=100, and 256 at Lp=1000.
            Here is the regression calculation using a logarithmic function to match: https://www.desmos.com/calculator/chdqwettsj
        - (gx, gy): This is the middle grey intersection constraint for the compression function.
            Scene-linear input value gx is mapped to display-linear output gy through the function.
            Why is gy set to 0.11696 at Lp=100? This matches the position of middle grey through the Rec709 system.
            We use this value for consistency with the Arri and TCAM Rec.1886 display rendering transforms.
      */

    // input scene-linear peak x intercept
    const float px = 256.0 * log(Lp) / log(100.0) - 128.0;
    // output display-linear peak y intercept
    const float py = Lp / 100.0;
    // input scene-linear middle grey x intercept
    const float gx = 0.18;
    // output display-linear middle grey y intercept
    const float gy = Lg / 100.0 * (1.0 + Lgb * log2(py));
    // s0 and s are input x scale for middle grey intersection constraint
    // m0 and m are output y scale for peak white intersection constraint
    const float s0 = quadratic_toe_compress(gy, toe, 1);
    const float m0 = quadratic_toe_compress(py, toe, 1);
    const float ip = 1.0 / p;
    const float s1 = pow(s0, ip);
    const float m1 = pow(m0, ip);
    const float s = (px * gx * (m1 - s1)) / (px * s1 - gx * m1);
    const float m = m1 * (s + px) / px;
    
    norm = max(0.0, norm);
    norm = hyperbolic_compress(norm, m, s, p, 0);
    norm = quadratic_toe_compress(norm, toe, 0) / py;

    // Apply purity boost
    float pb_m0 = 1.0 + pb;
    float pb_m1 = 2.0 - pb_m0;
    float pb_f = norm*(pb_m1 - pb_m0) + pb_m0;
    // Lerp from weights on bottom end to 1.0 at top end of tonescale
    float pb_L = dot(rgb, vec3(0.25, 0.7, 0.05)) * (1.0 - norm) + norm;
    float rats_mn = max(0.0, minf3(rgb));
    rgb = (rgb * pb_f + pb_L * (1.0 - pb_f)) * rats_mn + rgb * (1.0 - rats_mn);

    /* Purity Compression --------------------------------------- */
    // Apply purity compress using ccf by lerping to 1.0 in rgb ratios (peak achromatic)
    float ccf = norm / (spowf(m, p) / py); // normalize to enforce 0-1
    ccf = spowf(1.0 - ccf, pc_p);
    rgb = rgb * ccf + (1.0 - ccf);

    // "Density" - scale down intensity of colors to better fit in display-referred gamut volume 
    // and reduce discontinuities in high intensity high purity tristimulus.
    vec3 dn_r = clampminf3(1.0 - rgb, 0.0);
    // rgb = rgb * (dn_w.x * dn_r.x + 1.0 - dn_r.x) * (dn_w.y * dn_r.y + 1.0 - dn_r.y) * (dn_w.z * dn_r.z + 1.0 - dn_r.z);
    vec3 dn_wr = (dn_w * dn_r + 1.0 - dn_r);
    rgb *= dn_wr.x * dn_wr.y * dn_wr.z;

    /* Chroma Compression Hue Shift ------------------------------------------ *
        Since we compress chroma by lerping in a straight line towards 1.0 in rgb ratios, this can result in perceptual hue shifts
        due to the Abney effect. For example, pure blue compressed in a straight line towards achromatic appears to shift in hue towards purple.
        To combat this, and to add another important user control for image appearance, we add controls to curve the hue paths 
        as they move towards achromatic. We include only controls for primary colors: RGB. In my testing, it was of limited use to
        control hue paths for CMY.

        To accomplish this, we use the inverse of the chroma compression factor multiplied by the RGB hue angles as a factor
        for a lerp between the various rgb components.
    */

    float hs_mx = maxf3(rgb);
    vec3 hs_rgb = sdivf3f(rgb, hs_mx);
    float hs_mn = minf3(hs_rgb);
    hs_rgb = hs_rgb - hs_mn;
    // Narrow hue angles
    hs_rgb = vec3(
        min(1.0, max(0.0, hs_rgb.x - (hs_rgb.y + hs_rgb.z))),
        min(1.0, max(0.0, hs_rgb.y - (hs_rgb.x + hs_rgb.z))),
        min(1.0, max(0.0, hs_rgb.z - (hs_rgb.x + hs_rgb.y)))
    );
    hs_rgb = hs_rgb * (1.0 - ccf);

    // Apply hue shift to RGB Ratios
    vec3 rats_hs = vec3(
        rgb.x + hs_rgb.z * hs.z - hs_rgb.y * hs.y,
        rgb.y + hs_rgb.x * hs.x - hs_rgb.z * hs.z,
        rgb.z + hs_rgb.y * hs.y - hs_rgb.x * hs.x
    );

    // Mix hue shifted RGB ratios by ts, so that we shift where highlights were chroma compressed plus a bit.
    rgb = rgb * (1.0 - ccf) + rats_hs * ccf;

    // "Re-Saturate" using an inverse lerp
    sat_L = dot(rgb, sat_w);
    rgb = (sat_L * (sat_f - 1.0) + rgb) / sat_f;

    // last gamut compress for bottom end
    rgb.x = compress_powerptoe(rgb.x, 0.05, 1.0, 1.0, 0);
    rgb.y = compress_powerptoe(rgb.y, 0.05, 1.0, 1.0, 0);
    rgb.z = compress_powerptoe(rgb.z, 0.05, 1.0, 1.0, 0);

    // Apply tonescale to RGB Ratios
    rgb = rgb * norm;

    // Apply display scale
    rgb *= ds;

    // Clamp
    rgb = clampf3(rgb, 0.0, ds);

    // Apply inverse Display EOTF
    float eotf_p = 2.0 + eotf * 0.2;
    if ((eotf > 0) && (eotf < 4)) {
    	rgb = spowf3(rgb, 1.0 / eotf_p);
    } else if (eotf == 4) {
    	rgb = eotf_pq(rgb, 1);
    } else if (eotf == 5) {
    	rgb = eotf_hlg(rgb, 1);
    }

    return rgb;
}

#undef eotf

#endif // MODULES_OPENDT_OPENDRT_v03
