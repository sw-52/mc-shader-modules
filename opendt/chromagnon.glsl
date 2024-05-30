#ifndef MODULES_OPENDT_CHROMAGNON
#define MODULES_OPENDT_CHROMAGNON



/*------------------------------------------------
    The Chromagnon View Transform
    v0.0.7

    Written by Jed Smith
    https://github.com/jedypod/open-display-transform

    License: GPL v3
-------------------------------------------------*/

//
// Converted to GLSL by sw-52
//

#include "/modules/opendt/lib.glsl"


//////////////////////////////////////////////////////////////////////////////
/* Parameters */
// DEFINE_UI_PARAMS(in_gamut, input gamut, DCTLUI_COMBO_BOX, 13, {i_xyz, i_ap0, i_ap1, i_p3d65, i_rec2020, i_rec709, i_awg3, i_awg4, i_rwg, i_sgamut3, i_sgamut3cine, i_bmdwg, i_egamut, i_davinciwg}, {XYZ, ACES 2065-1, ACEScg, P3D65, Rec.2020, Rec.709, Arri Wide Gamut 3, Arri Wide Gamut 4, Red Wide Gamut RGB, Sony SGamut3, Sony SGamut3Cine, Blackmagic Wide Gamut, Filmlight E - Gamut, DaVinci Wide Gamut})

//  // Title  |  Default Value  | Range

// Tonescale
#define Lp   CHROMAGNON_LP            // Lp            100.0 [100.0 - 1000.0] 
#define Lg   CHROMAGNON_LG            // Lg            10.0  [4.0 - 20.0]
#define Lgb  CHROMAGNON_LG_BOOST      // Lg boost      0.12  [0.0 - 0.5]
#define offs CHROMAGNON_OFFSET        // offset        0.005 [0.0 - 0.02]
#define con  CHROMAGNON_CONTRAST      // contrast      1.8   [1.0 - 3.0]
#define conh CHROMAGNON_HIGH_CONTRAST // high contrast 0.0   [0.0 - 0.5]
#define w0   CHROMAGNON_HIGH_CLIP     // high clip     128.0 [16.0 - 256.0]

// "Split-toning" - adjust neutral axis
#define wtm CHROMAGNON_WHITE_TEMP // white temp 0.0 [-0.6 - 0.6]
#define wtn CHROMAGNON_WHITE_TINT // white tint 0.0 [-1.0 - 1.0]
#define gtm CHROMAGNON_GREY_TEMP  // grey temp  0.0 [-0.6 - 0.6]
#define gtn CHROMAGNON_GREY_TINT  // grey tint  0.0 [-1.0 - 1.0]

// Compress "Chroma" RGB
#define pc_p  CHROMAGNON_PC_P  // purity compress: slope    0.5 [0.0 - 1.0]
#define pc_x0 CHROMAGNON_PC_X0 //                amount     0.1 [0.0 - 0.5]
#define pc_t0 CHROMAGNON_PC_T0 //                threshold  0.1 [0.0 - 1.0]

// Rendering Matrix
#define rs_rnd CHROMAGNON_RND_S_R // render matrix: scale R 0.010 [0.0 - 0.3]
#define gs_rnd CHROMAGNON_RND_S_G //                      G 0.030 [0.0 - 0.3]
#define bs_rnd CHROMAGNON_RND_S_B //                      B 0.025 [0.0 - 0.3]
#define rr_rnd CHROMAGNON_RND_R_R //               rotate R 0.010 [-0.1 - 0.1]
#define gr_rnd CHROMAGNON_RND_R_G //                      G 0.000 [-0.1 - 0.1]
#define br_rnd CHROMAGNON_RND_R_B //                      B 0.015 [-0.1 - 0.1]

// Look Matrix
#define rs_lk CHROMAGNON_LK_R // Look MTX: R dist 0.03 [0.0 - 0.3]
#define gs_lk CHROMAGNON_LK_G // Look MTX: G dist 0.12 [0.0 - 0.3]
#define bs_lk CHROMAGNON_LK_B // Look MTX: B dist 0.02 [0.0 - 0.3]

#define clamp_rgb 0 // clamp [0, 1]
#define invert 0 // invert [0, 1]

#define eotf EOTF_lin

/*const mat3 in_to_xyz = matrix_rec2020_to_xyz;
mat3 xyz_to_display = inverse(matrix_rec2020_to_xyz);
mat3 xyz_to_in = xyz_to_display;
const mat3 display_to_xyz = matrix_rec2020_to_xyz;*/

/*

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

mat3 scrot_mtx(vec3 s, vec3 a) {
/* Create a 3x3 matrix that scales and rotates each primary,
      while preserving the neutral axis.
      s: amount to scale each primary
      a: amount to "rotate" each primary off axis
  */
    return mat3(
        vec3(1.0 - s.y - a.y - s.z + a.z, s.y + a.y, s.z - a.z),
        vec3(s.x + a.x, 1.0 - s.x - a.x - s.z - a.z, s.z + a.z),
        vec3(s.x - a.x, s.y - a.y, 1.0 - s.x + a.x - s.y + a.y)
    );
}


vec3 xlog(vec3 rgb, int inv) {
/* XLog - a simple log curve with no linear extension,
      intersection constraints at (w0, 1), and (0, y1), asymptote at o0.
      https://www.desmos.com/calculator/n64uurfttk
      https://colab.research.google.com/drive/1NwjaD0NzNMzQeNQqZECj33PdcYGkeBM4
  */
    //const float w0 = 128.0;
    const float k = -4.0;
    const float y1 = 0.0;
    const float o0 = -spowf(2.0, k);
    const float s0 = (1.0 - y1) / log(1.0 - w0 / o0);
    const float o1 = 1.0 - s0 * log(w0 - o0);
    if (inv == 1) {
        rgb.x = exp((rgb.x - o1) / s0) + o0;
        rgb.y = exp((rgb.y - o1) / s0) + o0;
        rgb.z = exp((rgb.z - o1) / s0) + o0;
    } else {
        rgb.x = s0 * log(rgb.x - o0) + o1;
        rgb.y = s0 * log(rgb.y - o0) + o1;
        rgb.z = s0 * log(rgb.z - o0) + o1;
    }
    return rgb;
}


float compress_powerptoe(float x, float p, float m, float t0, int inv) {
/* Variable slope compression function.
      p: Slope of the compression curve. Controls how compressed values are distributed.
         p=0.0 is a clip. p=1.0 is a hyperbolic curve.
      m: Compression amount. How far to reach outside of the gamut boundary to pull values in. m=0.0 is a clip.
      t0: Threshold point within gamut to start compression. t0=0.0 is a clip.
      https://www.desmos.com/calculator/igy3az7maq
  */

    float i = inv == 1 ? -1.0 : 1.0;
    return (x > t0) ? x : (x - t0) * spowf(1.0 + i * spowf((t0 - x) / (t0 - m), 1.0 / p), -p) + t0;
}


float contrast(float x, float p, float s0, float t0, float o0, float m0, int inv) {
/* Piecewise pivoted power function.
      p: power
      s0: scale which allows the pivot to be adjusted
      t0: threshold point where the high linear extension begins
      m0: slope for the linear extension
      o0: offset for the linear extension
  */
    if (inv == 1) return x < o0 ? spowf(x / s0, 1.0 / p) : (x - o0) / m0 + t0;
    else return (x < t0) ? s0 * spowf(x, p) : m0 * (x - t0) + o0;
}


float tonescale(float x, float m, float s, int inv) {
/*  Basic hyperbolic compression function.
        m: output y scale
        s: input x scale
        https://www.desmos.com/calculator/fkqf9jrju7
  */
    if (inv == 1) return (x > m) ? x : s * x / (m - x);
    else return m * x / (x + s);
}





vec3 chromagnon_transform(vec3 rgb) {
    // Input gamut conversion matrix (CAT02 chromatic adaptation to D65)


    const float ds = eotf == 4 ? 0.01 : eotf == 5 ? 0.1 : 100.0 / Lp;
    const float clamp_max = ds * Lp / 100.0;

    // Assemble scrot matrices
    mat3 look_mtx = scrot_mtx(vec3(rs_lk, gs_lk, bs_lk), vec3(0.0));
    mat3 look_mtx_inv = inverse(look_mtx);
    mat3 rnd_mtx = scrot_mtx(vec3(rs_rnd, gs_rnd, bs_rnd), vec3(rr_rnd, gr_rnd, br_rnd));
    mat3 rnd_mtx_inv = inverse(rnd_mtx);



    /* --------------------------- */
    /* Pre-Calculations for parameter space */

    // Precalculations for Purity Compress intersection constraint at (-x0, 0)
    const float pc_m0 = spowf((pc_t0 + max(1e-6, pc_x0)) / pc_t0, 1.0 / pc_p) - 1.0;
    const float pc_m = spowf(pc_m0, -pc_p)*(pc_t0 * spowf(pc_m0, pc_p) - pc_t0 - max(1e-6, pc_x0));

    /* Michaelis-Menten Contrast Tonescale Function
          Has 2 components: a contrast curve with linear extension,
          and a hyperbolic compression.
          https://www.desmos.com/calculator/de7q5mpmju
      */

    /* A contrast curve with an intersection constraint at (x0, y0), and
          a linear extension above t0.
          https://www.desmos.com/calculator/57lthhguq7
      */
    // Contrast precalculations
    const float x0 = 0.18 + offs;
    const float y0 = Lg / 100.0 * (1.0 + Lgb * log2(Lp / 100.0));
    const float t0 = x0 + conh;
    const float s0 = y0 * spowf(x0, -con);
    const float o0 = s0 * spowf(t0, con);
    const float m0 = t0 == 0.0 ? 0.0 : con * o0 / t0;

    /* Super simple Michaelis-Menten tonescale function with
          intersection constraints at (w0, w1) and (y0, y0).
          There is no scene-linear scale since we do that earlier in the contrast function.
          https://www.desmos.com/calculator/gesdthnhd3
      */

    // Tonescale precalculations
    const float w1 = Lp / 100.0;
    const float m1 = (w1 * (m0 * w0 - y0)) / (m0 * w0 - w1);
    const float wtm_abs = wtm > 0.0 ? wtm : -wtm;
    const vec3 m = vec3(
        m1*spowf(2.0, wtm - wtm_abs),
        m1*spowf(2.0, wtm * wtn - wtm_abs),
        m1*spowf(2.0, -wtm - wtm_abs)
    );
    const vec3 s1 = vec3(m.x - y0, m.y - y0, m.z - y0);
    const vec3 s = vec3(
        s1.x * spowf(2.0, -gtm),
        s1.y * spowf(2.0, -gtm * gtn),
        s1.z * spowf(2.0, gtm)
    );





    if (invert == 0) {
        /* Forward Rendering Code ------------------------------------------ */

        // Convert into display gamut
        rgb = vdot(in_to_xyz, rgb);
        rgb = vdot(xyz_to_display, rgb);

        // Offset
        rgb += offs;

        /* Purity Compress
                A low complexity per-channel purity / "chroma" / "saturation" compression.
                Note we could expose the individual compression parameters to achieve hue
                control, but we skip this since we will use the scrot matrices for this purpose.
            */
        if (pc_p < 0.01 || pc_t0 < 0.001) {
            // Just clamp min of slope = 0
            rgb = clampminf3(rgb, 0.0);
        } else {
            float mx = maxf3(rgb);
            rgb = sdivf3f(rgb, mx);
            rgb.x = compress_powerptoe(rgb.x, pc_p, pc_m, pc_t0, 0);
            rgb.y = compress_powerptoe(rgb.y, pc_p, pc_m, pc_t0, 0);
            rgb.z = compress_powerptoe(rgb.z, pc_p, pc_m, pc_t0, 0);
            rgb = rgb * mx;
        }

        // Rendering matrix
        rgb = vdot(rnd_mtx, rgb);

        // Clamp negatives
        rgb = clampminf3(rgb, 0.0);

        // Look matrix in log
        rgb = xlog(rgb, 0);
        rgb = vdot(look_mtx, rgb);
        rgb = xlog(rgb, 1);

        // Apply tonescale
        rgb.x = contrast(rgb.x, con, s0, t0, o0, m0, 0);
        rgb.y = contrast(rgb.y, con, s0, t0, o0, m0, 0);
        rgb.z = contrast(rgb.z, con, s0, t0, o0, m0, 0);

        rgb.x = tonescale(rgb.x, m.x, s.x, 0);
        rgb.y = tonescale(rgb.y, m.y, s.y, 0);
        rgb.z = tonescale(rgb.z, m.z, s.z, 0);

        /* Display scale
                Tonescale outputs nits/100. 10 nits = 0.1, 1000 nits = 10.0
                This scale either normalizes into a 0-1 range for SDR
                or into the 10,000 nit PQ container, or 1,000 nit HLG container.
            */
        rgb *= ds;

        // Clamp
        if (clamp_rgb == 1) rgb = clampf3(rgb, 0.0, clamp_max);

        // Apply inverse Display EOTF
        float eotf_p = 2.0 + eotf * 0.2;
        if ((eotf > 0) && (eotf < 4)) {
            rgb = spowf3(rgb, 1.0 / eotf_p);
        } else if (eotf == 4) {
            rgb = eotf_pq(rgb, 1);
        } /*else if (eotf == 5) {
            rgb = eotf_hlg(rgb, 1);
        }*/

        return rgb;

    } else {
        /* Inverse Rendering Code ------------------------------------------ */

        rgb = clampf3(rgb, 0.0, 1.0);

        // Apply Display EOTF
        float eotf_p = 2.0 + eotf * 0.2;
        if ((eotf > 0) && (eotf < 4)) {
            rgb = spowf3(rgb, eotf_p);
        } else if (eotf == 4) {
            rgb = eotf_pq(rgb, 0);
        } /*else if (eotf == 5) {
            rgb = eotf_hlg(rgb, 0);
        }*/

        // Inverse Display scale
        rgb /= ds;

        // Invert tonescale
        rgb.x = tonescale(rgb.x, m.x, s.x, 1);
        rgb.y = tonescale(rgb.y, m.y, s.y, 1);
        rgb.z = tonescale(rgb.z, m.z, s.z, 1);
        rgb.x = contrast(rgb.x, con, s0, t0, o0, m0, 1);
        rgb.y = contrast(rgb.y, con, s0, t0, o0, m0, 1);
        rgb.z = contrast(rgb.z, con, s0, t0, o0, m0, 1);

        // Inverse Look matrix in XLog
        rgb = xlog(rgb, 0);
        rgb = vdot(look_mtx_inv, rgb);
        rgb = xlog(rgb, 1);

        // Inverse rendering matrix
        rgb = vdot(rnd_mtx_inv, rgb);

        /* Inverse Low-complexity per-channel "chroma" compression */
        float mx = maxf3(rgb);
        rgb = sdivf3f(rgb, mx);
        rgb.x = compress_powerptoe(rgb.x, pc_p, pc_m, pc_t0, 1);
        rgb.y = compress_powerptoe(rgb.y, pc_p, pc_m, pc_t0, 1);
        rgb.z = compress_powerptoe(rgb.z, pc_p, pc_m, pc_t0, 1);
        rgb = rgb * mx;

        // Inverse Offset
        rgb -= offs;

        // Convert from display gamut back to scene gamut
        rgb = vdot(inverse(xyz_to_display), rgb);
        rgb = vdot(inverse(in_to_xyz), rgb);

        return rgb;

    }
}

#undef Lp
#undef Lg
#undef Lgb
#undef offs
#undef con
#undef conh
#undef w0

#undef wtm
#undef wtn
#undef gtm
#undef gtn

#undef pc_p
#undef pc_x0
#undef pc_t0

#undef rs_rnd
#undef gs_rnd
#undef bs_rnd
#undef rr_rnd
#undef gr_rnd
#undef br_rnd

#undef rs_lk
#undef gs_lk
#undef bs_lk

#undef clamp_rgb
#undef invert

#undef eotf

//#undef spowf

#endif // MODULES_OPENDT_CHROMAGNON
