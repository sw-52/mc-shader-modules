#ifndef MODULES_OPENDT_OPENDRT_COMMON
#define MODULES_OPENDT_OPENDRT_COMMON



/*  OpenDRT -------------------------------------------------/
      Written by Jed Smith
      https://github.com/jedypod/open-display-transform

      License: GPL v3
-------------------------------------------------*/

//
// Converted to GLSL by sw-52
//

#include "/modules/opendt/lib.glsl"



/* Functions for the OpenDRT v0.3.x Transform ---------------------------------------- */

float compress_powerptoe_p(float x, const float p, const float x0, const float t0, int inv) {
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



#endif // MODULES_OPENDT_OPENDRT_COMMON
