#ifndef MODULES_OPENDT_TESSERACT
#define MODULES_OPENDT_TESSERACT



/*------------------------------------------------
    Tesseract View Transform
    v0.0.3
    
    Written by Jed Smith
    https://github.com/jedypod/open-display-transform

    License: GPL v3
-------------------------------------------------*/

//
// Converted to GLSL by sw-52
//

#include "/modules/opendt/lib.glsl"


//////////////////////////////////////////////////////////////////////////////

//DEFINE_UI_PARAMS(in_gamut, input gamut, DCTLUI_COMBO_BOX, 13, {i_xyz, i_ap0, i_ap1, i_p3d65, i_rec2020, i_rec709, i_awg3, i_awg4, i_rwg, i_sgamut3, i_sgamut3cine, i_bmdwg, i_egamut, i_davinciwg}, {XYZ, ACES 2065-1, ACEScg, P3D65, Rec.2020, Rec.709, Arri Wide Gamut 3, Arri Wide Gamut 4, Red Wide Gamut RGB, Sony SGamut3, Sony SGamut3Cine, Blackmagic Wide Gamut, Filmlight E - Gamut, DaVinci Wide Gamut})
// Tonescale options
#define Lp   TESSERACT_LP       // Lp            100.0 [100.0 - 1000.0] 
#define Lg   TESSERACT_LG       // Lg            10.0  [4.0 - 20.0]
#define offs TESSERACT_OFFSET   // offset        0.005 [0.0 - 0.02]
#define con  TESSERACT_CONTRAST // contrast      1.6   [1.0 - 2.0]

// Color options
#define pb   TESSERACT_PB   // Purity      0.6 [0.0 - 1.0]
#define hl_r TESSERACT_HL_R // hue low r   0.3 [-1.0, 1.0]
#define hl_g TESSERACT_HL_G // hue low g   0.2 [-1.0, 1.0]
#define hl_b TESSERACT_HL_B // hue low b   0.2 [-1.0, 1.0]
#define hh_r TESSERACT_HH_R // hue high r  0.5 [-1.0, 1.0]
#define hh_g TESSERACT_HH_G // hue high g  0.0 [-1.0, 1.0]
#define hh_b TESSERACT_HH_B // hue high b -0.5 [-1.0, 1.0]

#define clamp_rgb 0 // clamp 0 [0 1]

#define eotf EOTF_lin // EOTF EOTF_lin [EOTF_lin EOTF_srgb EOTF_rec1886 EOTF_dci EOTF_pq EOTF_hlg]

/*const mat3 in_to_xyz = matrix_rec2020_to_xyz;
mat3 xyz_to_display = inverse(matrix_rec2020_to_xyz);
mat3 xyz_to_in = xyz_to_display;
const mat3 display_to_xyz = matrix_rec2020_to_xyz;*/

/*

// Output options

DEFINE_UI_PARAMS(clamp, clamp, DCTLUI_CHECK_BOX, 1);
DEFINE_UI_PARAMS(display_gamut, display gamut, DCTLUI_COMBO_BOX, 0, {Rec709, P3D65, Rec2020}, {Rec.709, P3 D65, Rec.2020})
DEFINE_UI_PARAMS(EOTF, display eotf, DCTLUI_COMBO_BOX, 2, {lin, srgb, rec1886, dci, pq, hlg}, {Linear, 2.2 Power - sRGB Display, 2.4 Power - Rec .1886, 2.6 Power - DCI, ST 2084 PQ, HLG})

// Gamut Conversion Matrices
#define matrix_ap0_to_xyz mat3(vec3(0.93863094875f, -0.00574192055f, 0.017566898852f), vec3(0.338093594922f, 0.727213902811f, -0.065307497733f), vec3(0.000723121511f, 0.000818441849f, 1.0875161874f))
#define matrix_ap1_to_xyz mat3(vec3(0.652418717672f, 0.127179925538f, 0.170857283842f), vec3(0.268064059194f, 0.672464478993f, 0.059471461813f), vec3(-0.00546992851f, 0.005182799977f, 1.08934487929f))
#define matrix_rec709_to_xyz mat3(vec3(0.412390917540f, 0.357584357262f, 0.180480793118f), vec3(0.212639078498f, 0.715168714523f, 0.072192311287f), vec3(0.019330825657f, 0.119194783270f, 0.950532138348f))
#define matrix_p3d65_to_xyz mat3(vec3(0.486571133137f, 0.265667706728f, 0.198217317462f), vec3(0.228974640369f, 0.691738605499f, 0.079286918044f), vec3(-0.000000000000f, 0.045113388449, 1.043944478035f))
#define matrix_rec2020_to_xyz mat3(vec3(0.636958122253f, 0.144616916776f, 0.168880969286f), vec3(0.262700229883f, 0.677998125553f, 0.059301715344f), vec3(0.000000000000f, 0.028072696179, 1.060985088348f))
#define matrix_arriwg3_to_xyz mat3(vec3(0.638007619284f, 0.214703856337f, 0.097744451431f), vec3(0.291953779f, 0.823841041511f, -0.11579482051f), vec3(0.002798279032f, -0.067034235689f, 1.15329370742f))
#define matrix_arriwg4_to_xyz mat3(vec3(0.704858320407f, 0.12976029517f, 0.115837311474f), vec3(0.254524176404f, 0.781477732712f, -0.036001909116f), vec3(0.0f, 0.0f, 1.08905775076f))
#define matrix_redwg_to_xyz mat3(vec3(0.735275208950f, 0.068609409034f, 0.146571278572f), vec3(0.286694079638f, 0.842979073524f, -0.129673242569f), vec3(-0.079680845141f, -0.347343206406, 1.516081929207f))
#define matrix_sonysgamut3_to_xyz mat3(vec3(0.706482713192f, 0.128801049791f, 0.115172164069f), vec3(0.270979670813f, 0.786606411221f, -0.057586082034f), vec3(-0.009677845386f, 0.004600037493f, 1.09413555865f))
#define matrix_sonysgamut3cine_to_xyz mat3(vec3(0.599083920758f, 0.248925516115f, 0.102446490178f), vec3(0.215075820116f, 0.885068501744f, -0.100144321859f), vec3(-0.032065849545f, -0.027658390679f, 1.14878199098f))
#define matrix_bmdwg_to_xyz mat3(vec3(0.606538414955f, 0.220412746072f, 0.123504832387f), vec3(0.267992943525f, 0.832748472691f, -0.100741356611f), vec3(-0.029442556202f, -0.086612440646, 1.205112814903f))
#define matrix_egamut_to_xyz mat3(vec3(0.705396831036f, 0.164041340351f, 0.081017754972f), vec3(0.280130714178f, 0.820206701756f, -0.100337378681f), vec3(-0.103781513870f, -0.072907261550, 1.265746593475f))
#define matrix_davinciwg_to_xyz mat3(vec3(0.700622320175f, 0.148774802685f, 0.101058728993f), vec3(0.274118483067f, 0.873631775379f, -0.147750422359f), vec3(-0.098962903023f, -0.137895315886, 1.325916051865f))

#define matrix_xyz_to_rec709 mat3(vec3(3.2409699419f, -1.53738317757f, -0.498610760293f), vec3(-0.969243636281f, 1.87596750151f, 0.041555057407f), vec3(0.055630079697f, -0.203976958889f, 1.05697151424f))
#define matrix_xyz_to_p3d65 mat3(vec3(2.49349691194f, -0.931383617919f, -0.402710784451f), vec3(-0.829488969562f, 1.76266406032f, 0.023624685842f), vec3(0.035845830244f, -0.076172389268f, 0.956884524008f))
#define matrix_xyz_to_rec2020 mat3(vec3(1.71665118797f, -0.355670783776f, -0.253366281374f), vec3(-0.666684351832f, 1.61648123664f, 0.015768545814f), vec3(0.017639857445f, -0.042770613258f, 0.942103121235f))
*/


////////////////////////////////////////////////////////////////////////////////////////////
float compress_quadratic_lift(float x, float t0, float x0, int inv) {
  /* Quadratic toe compression function: https://www.desmos.com/calculator/tkwkrrbyy5
      t0: Threshold splice point below which values are compressed
      x0: How much to compress: -x0 will be compressed to 0.0
  */

    float s0 = sqrt(max(1e-6, x0)) / -t0;
    float o0 = 1.0 / (4.0 * s0 * s0) + t0;
    float o1 = t0 - sqrt(o0 - t0) / s0;
    
    if (inv == 1) return x > t0 ? x : o0 - spowf(s0, 2.0) * spowf(o1 - x, 2.0);
    else          return x > t0 ? x : o1 + sqrt(o0 - x) / s0;
}



float tonescale(float x, float Lp0, float Lg0, float off, float p, int inv) {
  /* Michaelis-Menten Contrast Tonescale Function
      Has 2 components: a contrast curve with linear extension, and a hyperbolic compression.
      Always outputs 0-1 despite Lp/Lg normalize
      https://www.desmos.com/calculator/de7q5mpmju

      Lp0: peak display luminance in nits
      Lg0: luminance of middle grey in nits
      off: black offset (this is part of the tonescale constraint but applied upstream)
      p: contrast
  
    Piecewise pivoted power function: https://www.desmos.com/calculator/401xfuag0q
      p: power
      x0: x coordinate for intersection constraint, and pivot point where the linear extension begins
      y0: y coordinate for intersection constraint
      s0: scale which allows the pivot to be adjusted
      m0: slope for the linear extension
      o0: offset for the linear extension

    Super simple Michaelis-Menten tonescale function with intersection constraints at (y0, y0) and (w0, w1).
      There is no scene-linear scale since we do that earlier in the contrast function.
      https://www.desmos.com/calculator/gesdthnhd3
      m: output y scale
      s: input x scale
  */

    // Amount of stops per stop of Lp increase to increase middle grey; this puts Lg at about 14 nits @ Lp=1000 nits
    const float Lgb = 0.12;
    float x0 = 0.18 + off;
    float y0 = Lg0 / Lp0 * (1.0 + Lgb * log(Lp0 / 100.0) / log(2.0)); // Lg / 100 would normalize middle grey to be at display lin val

    float s0 = y0 * spowf(x0, -p);
    float o0 = s0 * spowf(x0, p);
    float m0 = x0 == 0.0 ? 0.0 : p * o0 / x0;

    float w0 = 256.0 * log(Lp0) / log(100.0) - 128.0;
    const float w1 = 1.0; // Lp/100.0 would normalize 1.0 to 100 nits and 10.0 to 1000 nits
    float m = (w1 * (m0 * w0 - y0)) / (m0 * w0 - w1);
    float s = m - y0;
    
    if (inv == 1) {
        x = x > m ? x : s * x / (m - x);
        return x < o0 ? spowf(x / s0, 1.0 / p) : (x - o0) / m0 + x0;
    } else {
        x = x < x0 ? s0 * spowf(x, p) : m0 * (x - x0) + o0;
        return m * x / (x + s);
    }
}





vec3 tesseract_transform(vec3 rgb) {


    // const float ds = eotf == 4 ? 0.01 : eotf == 5 ? 0.1 : 100.0/Lp;
    const float ds = eotf == 4 ? Lp / 10000.0 : eotf == 5 ? Lp / 1000.0 : 1.0;



    /* --------------------------- */

    /* Forward Rendering Code ------------------------------------------ */

    // Convert into display gamut
    rgb = vdot(in_to_xyz, rgb);
    rgb = vdot(xyz_to_display, rgb);

    // Offset
    rgb += offs;

    float mx = maxf3(rgb);
    float mn = minf3(rgb);

    // RGB Ratios
    rgb = sdivf3f(rgb, mx);

    // Preserve sanity in shadow grain by limiting crazy ratios
    rgb = clampf3(rgb, -2.0, 2.0);

    /* Purity Compression --------------------------------------- */
    // rgb purity compression strength, tuned for common camera gamuts
    vec3 pc_rats = rgb;

    // minrgb with out of gamut values compressed
    pc_rats.x = compress_quadratic_lift(pc_rats.x, 0.2, 0.3, 0);
    pc_rats.y = compress_quadratic_lift(pc_rats.y, 0.2, 0.3, 0);
    pc_rats.z = compress_quadratic_lift(pc_rats.z, 0.2, 0.03, 0);

    // minrgb with out of gamut values compressed
    float pc_rats_mn = minf3(pc_rats);

    // remove achromatic
    pc_rats = pc_rats - pc_rats_mn;
    // restore full range for tonescale curve
    pc_rats_mn *= mx;

    // the factor to mix back to maxrgb
    vec3 pc_s = vec3(0.03, 0.16, 0.06);
    float pc_f = dot(pc_s, pc_rats);

    // minrgb mix to maxrgb by pc_f
    mn = max(0.0, pc_rats_mn * (1.0 - pc_f) + mx * pc_f);

    // Apply tonescale to our modified minrgb to use for purity compression
    mn = tonescale(mn, Lp, Lg, offs, con, 0);
    // Apply tonescale to maxrgb
    mx = tonescale(mx, Lp, Lg, offs, con, 0);

    // Apply purity boost
    float pb_m0 = 1.0 + pb;
    float pb_m1 = 2.0 - pb_m0;
    float pb_f = mx * (pb_m1 - pb_m0) + pb_m0;
    float rats_mn = max(0.0, minf3(rgb));
    rgb = (rgb * pb_f + 1.0 - pb_f) * rats_mn + rgb * (1.0 - rats_mn);

    // Apply purity compress using mn by lerping to 1.0 in rgb ratios (peak achromatic)
    rgb = rgb * (1.0 - mn) + mn;


    float mn_pc = minf3(rgb);
    vec3 ha_rgb = rgb - mn_pc;
    // Narrow hue angles
    ha_rgb = vec3(
        min(1.0, max(0.0, ha_rgb.x - (ha_rgb.y + ha_rgb.z))),
        min(1.0, max(0.0, ha_rgb.y - (ha_rgb.x + ha_rgb.z))),
        min(1.0, max(0.0, ha_rgb.z - (ha_rgb.x + ha_rgb.y)))
    );


    // Apply HueContrast Opposing
    vec3 ha_rgb_hl = spowf3(ha_rgb, 2.0); // only for hue contrast not hue shift
    vec3 hc_m0 = vec3(1.0 - hl_r, 1.0 - hl_g, 1.0 - hl_b);
    vec3 hc_m1 = vec3(1.0 + hl_r, 1.0 + hl_g, 1.0 + hl_b);

    // Remap with a lift/multiply function: https://www.desmos.com/calculator/eaxwxrhrev
    vec3 f = vec3(
        mn * (hc_m1.x - hc_m0.x) + hc_m0.x,
        mn * (hc_m1.y - hc_m0.y) + hc_m0.y,
        mn * (hc_m1.z - hc_m0.z) + hc_m0.z
    );

    // Apply hue contrast
    rgb = vec3(
        rgb.x * f.y * ha_rgb_hl.y + rgb.x * f.z * ha_rgb_hl.z + rgb.x * (1.0 - (ha_rgb_hl.y + ha_rgb_hl.z)),
        rgb.y * f.x * ha_rgb_hl.x + rgb.y * f.z * ha_rgb_hl.z + rgb.y * (1.0 - (ha_rgb_hl.x + ha_rgb_hl.z)),
        rgb.z * f.x * ha_rgb_hl.x + rgb.z * f.y * ha_rgb_hl.y + rgb.z * (1.0 - (ha_rgb_hl.x + ha_rgb_hl.y))
    );

    // Apply hue shift to RGB Ratios mixed by mn
    rgb = rgb * (1.0 - mn) + mn * vec3(
        rgb.x + ha_rgb.z * hh_b - ha_rgb.y * hh_g, 
        rgb.y + ha_rgb.x * hh_r - ha_rgb.z * hh_b, 
        rgb.z + ha_rgb.y * hh_g - ha_rgb.x * hh_r
    );

    /* Final gamut compress to smooth discontinuities on the low end at the expense of saturation */
    rgb.x = compress_quadratic_lift(rgb.x, 0.06, 0.02, 0);
    rgb.y = compress_quadratic_lift(rgb.y, 0.06, 0.06, 0);
    rgb.z = compress_quadratic_lift(rgb.z, 0.06, 0.02, 0);

    // Depart from RGB Ratios
    rgb = mx * rgb;


    /* Display scale
        Tonescale outputs nits/100. 10 nits = 0.1, 1000 nits = 10.0
        This scale either normalizes into a 0-1 range for SDR
        or into the 10,000 nit PQ container, or 1,000 nit HLG container.
    */
    rgb *= ds;

    // Clamp
    if (clamp_rgb == 1)
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

#undef Lp
#undef Lg
#undef offs
#undef con

#undef pb
#undef hl_r
#undef hl_g
#undef hl_b
#undef hh_r
#undef hh_g
#undef hh_b

#undef clamp_rgb

#undef eotf

//#undef spowf

#endif // MODULES_OPENDT_TESSERACT