/*
 *  Copyright 2016 ARTED developers
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#define ENABLE_STENCIL_CODE_WITH_PADDING

/* Hand-Code Vector processing for Fujitsu SPARC64 XIfx (HPC-ACE2) */

#include <complex.h>
#include <fjintrin.h>
#include "../interop.h"

void hpsi1_rt_stencil_( double         const* restrict A_
                      , double         const           B[restrict NLx][NLy][NLz]
                      , double         const           C[restrict 12]
                      , double         const           D[restrict 12]
                      , double complex const           E[restrict PNLx][PNLy][PNLz]
                      , double complex                 F[restrict PNLx][PNLy][PNLz]
)
{
  const  double                  A = *A_;
  double         const* restrict b;
  double complex const* restrict e;
  double complex      * restrict f;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  _fjsp_simd_mode_4();

  _fjsp_v4r8 c[12], d[12];
  for(n = 0 ; n < 12 ; ++n) {
    double a = C[n] * -0.5;
    double b = D[n];
    c[n] = _fjsp_broadcastload_v4r8(&a);
    d[n] = _fjsp_broadcastload_v4r8(&b);
  }

  const _fjsp_v4r8 at   = _fjsp_broadcastload_v4r8(&A);
  const _fjsp_v4r8 INV  = _fjsp_set_v4r8(-1, 1, -1, 1);
  const _fjsp_v4r8 sel  = _fjsp_set_v4i8(2, 3, 0, 1);
  const _fjsp_v4r8 cast = _fjsp_set_v4i4(1, 1, 0, 0);

  const _fjsp_v4r8 pnly  = _fjsp_broadcastload_v4i4(&PNLy);
  const _fjsp_v4r8 pnlz  = _fjsp_broadcastload_v4i4(&PNLz);
  const _fjsp_v4r8 pnlyz = _fjsp_mul_v4i4(pnly, pnlz);

#ifdef ARTED_DOMAIN_POWER_OF_TWO
  const _fjsp_v4r8 nlx = _fjsp_broadcastload_v4i4(&NLx);
  const _fjsp_v4r8 nly = _fjsp_broadcastload_v4i4(&NLy);

  const int MX = NLx - 1;
  const int MY = NLy - 1;
  const _fjsp_v4r8 mx = _fjsp_broadcastload_v4i4(&MX);
  const _fjsp_v4r8 my = _fjsp_broadcastload_v4i4(&MY);

  const _fjsp_v4r8 mtbl = _fjsp_set_v4i4(-1, -2, -3, -4);
  const _fjsp_v4r8 ptbl = _fjsp_set_v4i4( 4,  3,  2,  1);

  const _fjsp_v4r8 ymtbl = _fjsp_add_v4i4(nly, mtbl);
  const _fjsp_v4r8 yptbl = _fjsp_add_v4i4(nly, ptbl);
  const _fjsp_v4r8 xmtbl = _fjsp_add_v4i4(nlx, mtbl);
  const _fjsp_v4r8 xptbl = _fjsp_add_v4i4(nlx, ptbl);
#endif

  int idx[8], idy[8];

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix)
#else
  for(ix = 0 ; ix < NLx ; ++ix)
#endif
  {
    _fjsp_v4r8 tix = _fjsp_broadcastload_v4i4(&ix);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
    _fjsp_v4r8 txm = _fjsp_sub_v4i4(tix, _fjsp_and_v4i4(_fjsp_add_v4i4(tix, xmtbl), mx));
    _fjsp_v4r8 txp = _fjsp_sub_v4i4(tix, _fjsp_and_v4i4(_fjsp_add_v4i4(tix, xptbl), mx));
#else
    _fjsp_v4r8 txm = _fjsp_sub_v4i4(tix, _fjsp_load_v4i4(modx + (ix - 4 + NLx)));
    _fjsp_v4r8 txp = _fjsp_sub_v4i4(tix, _fjsp_load_v4i4(modx + (ix + 1 + NLx)));
#endif
               txm = _fjsp_mul_v4i4(pnlyz, txm);
               txp = _fjsp_mul_v4i4(pnlyz, txp);
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    _fjsp_v4r8 tiy = _fjsp_broadcastload_v4i4(&iy);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
    _fjsp_v4r8 tym = _fjsp_sub_v4i4(tiy, _fjsp_and_v4i4(_fjsp_add_v4i4(tiy, ymtbl), my));
    _fjsp_v4r8 typ = _fjsp_sub_v4i4(tiy, _fjsp_and_v4i4(_fjsp_add_v4i4(tiy, yptbl), my));
#else
    _fjsp_v4r8 tym = _fjsp_sub_v4i4(tiy, _fjsp_load_v4i4(mody + (iy - 4 + NLy)));
    _fjsp_v4r8 typ = _fjsp_sub_v4i4(tiy, _fjsp_load_v4i4(mody + (iy + 1 + NLy)));
#endif
               tym = _fjsp_mul_v4i4(pnlz, tym);
               typ = _fjsp_mul_v4i4(pnlz, typ);

    b = &B[ix][iy][0];
    e = &E[ix][iy][0];
    f = &F[ix][iy][0];

    const double complex *ee = e + 128;
    asm("prefetch [%0],20": :"r"(ee));

    for(iz = 0 ; iz < NLz ; iz += 2)
    {
      _fjsp_v4r8 tiz = _fjsp_broadcastload_v4i4(&iz);

      _fjsp_store_v4i4(idy+0, _fjsp_sub_v4i4(tiz, tym));
      _fjsp_store_v4i4(idy+4, _fjsp_sub_v4i4(tiz, typ));
      _fjsp_store_v4i4(idx+0, _fjsp_sub_v4i4(tiz, txm));
      _fjsp_store_v4i4(idx+4, _fjsp_sub_v4i4(tiz, txp));

#define STENCIL_CALC(MM,PP,CC,DD,TT,UT) \
      TT = _fjsp_madd_v4r8(CC, _fjsp_add_v4r8(PP, MM), TT); \
      UT = _fjsp_madd_v4r8(DD, _fjsp_sub_v4r8(PP, MM), UT);

#define STENCIL_LOAD(IDM,IDP,MM,PP) \
      MM = _fjsp_load_v4r8((double *) (e + IDM)); \
      PP = _fjsp_load_v4r8((double *) (e + IDP));

      _fjsp_v4r8 tt0 = _fjsp_setzero_v4r8();
      _fjsp_v4r8 ut0 = _fjsp_setzero_v4r8();
      _fjsp_v4r8 tt1 = _fjsp_setzero_v4r8();
      _fjsp_v4r8 ut1 = _fjsp_setzero_v4r8();
      _fjsp_v4r8 tt2 = _fjsp_setzero_v4r8();
      _fjsp_v4r8 ut2 = _fjsp_setzero_v4r8();

      _fjsp_v4r8 ez;
      _fjsp_v4r8 z0, z1, z2, z3, z4, z5, z6, z7;
      _fjsp_v4r8 y0, y1, y2, y3, y4, y5, y6, y7;
      _fjsp_v4r8 x0, x1, x2, x3, x4, x5, x6, x7;

      ez = _fjsp_load_v4r8((double *)(e + iz));

#ifdef ARTED_DOMAIN_POWER_OF_TWO
      z0 = _fjsp_load_v4r8((double *)(e + ((iz - 4 + NLz) & (NLz - 1))));
      z1 = _fjsp_load_v4r8((double *)(e + ((iz - 2 + NLz) & (NLz - 1))));
      z2 = _fjsp_load_v4r8((double *)(e + ((iz + 2 + NLz) & (NLz - 1))));
      z3 = _fjsp_load_v4r8((double *)(e + ((iz + 4 + NLz) & (NLz - 1))));
#else
      z0 = _fjsp_load_v4r8((double *)(e + modz[iz - 4 + NLz]));
      z1 = _fjsp_load_v4r8((double *)(e + modz[iz - 2 + NLz]));
      z2 = _fjsp_load_v4r8((double *)(e + modz[iz + 2 + NLz]));
      z3 = _fjsp_load_v4r8((double *)(e + modz[iz + 4 + NLz]));
#endif
      asm("### MAGIC"); /* !!! don't remove !!! */

      z6 = _fjsp_ecsl_v4(z0, z1, 2);
      z4 = _fjsp_ecsl_v4(z1, ez, 2);
      z5 = _fjsp_ecsl_v4(ez, z2, 2);
      z7 = _fjsp_ecsl_v4(z2, z3, 2);

      STENCIL_LOAD(idy[3], idy[4], y0, y4)
      STENCIL_LOAD(idy[2], idy[5], y1, y5)
      STENCIL_LOAD(idy[1], idy[6], y2, y6)
      STENCIL_LOAD(idy[0], idy[7], y3, y7)

      STENCIL_LOAD(idx[3], idx[4], x0, x4)
      STENCIL_LOAD(idx[2], idx[5], x1, x5)
      STENCIL_LOAD(idx[1], idx[6], x2, x6)
      STENCIL_LOAD(idx[0], idx[7], x3, x7)

      /* z-dimension (unit stride) */
      STENCIL_CALC(z0, z3, c[11], d[11], tt0, ut0);
      STENCIL_CALC(z1, z2, c[ 9], d[ 9], tt1, ut1);
      STENCIL_CALC(z4, z5, c[ 8], d[ 8], tt2, ut2);
      STENCIL_CALC(z6, z7, c[10], d[10], tt0, ut0);

      /* y-dimension (NLz stride) */
      STENCIL_CALC(y0, y4, c[4], d[4], tt1, ut1);
      STENCIL_CALC(y1, y5, c[5], d[5], tt2, ut2);
      STENCIL_CALC(y2, y6, c[6], d[6], tt0, ut0);
      STENCIL_CALC(y3, y7, c[7], d[7], tt1, ut1);

      /* x-dimension (NLy*NLz stride)  */
      STENCIL_CALC(x0, x4, c[0], d[0], tt2, ut2);
      STENCIL_CALC(x1, x5, c[1], d[1], tt0, ut0);
      STENCIL_CALC(x2, x6, c[2], d[2], tt1, ut1);
      STENCIL_CALC(x3, x7, c[3], d[3], tt2, ut2);

      _fjsp_v4r8 tt3 = _fjsp_add_v4r8(tt0, tt1);
      _fjsp_v4r8 ut3 = _fjsp_add_v4r8(ut0, ut1);
      _fjsp_v4r8 tt  = _fjsp_add_v4r8(tt2, tt3);
      _fjsp_v4r8 ut  = _fjsp_add_v4r8(ut2, ut3);

      _fjsp_v4r8 bv = _fjsp_load_v4r8(b + iz);
      _fjsp_v4r8 bt = _fjsp_eperm_v4(bv, cast);
      _fjsp_v4r8 ab = _fjsp_add_v4r8(at, bt);
      _fjsp_v4r8 tu = _fjsp_eperm_v4(ut, sel);

      _fjsp_v4r8 v2 = _fjsp_madd_v4r8(ab, ez, tt);
      _fjsp_v4r8 v1 = _fjsp_madd_v4r8(INV, tu, v2);

      _fjsp_store_v4r8((double *) (f + iz), v1);
    } /* NLz */
  } /* NLy */
  } /* NLx */

  _fjsp_simd_mode_2();
}
