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

/* Hand-Code Vector processing for Fujitsu SPARC64 XIfx (HPC-ACE2) */

#include <complex.h>
#include <fjintrin.h>
#include "../interop.h"

void current_stencil_( double         const           C[restrict 12]
                     , double complex const           E[restrict NLx][NLy][NLz]
                     , double              * restrict F
                     , double              * restrict G
                     , double              * restrict H
)
{
  double complex const* restrict e;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  _fjsp_simd_mode_4();

  _fjsp_v4r8 tt[3];
  for(n = 0 ; n < 3 ; ++n)
    tt[n] = _fjsp_setzero_v4r8();

  _fjsp_v4r8 c[12];
  for(n = 0 ; n < 12 ; ++n)
    c[n] = _fjsp_broadcastload_v4r8(C+n);

  const _fjsp_v4r8 INV  = _fjsp_set_v4u8(1ULL << 63, 0, 1ULL << 63, 0);
  const _fjsp_v4r8 sel  = _fjsp_set_v4i8(2, 3, 0, 1);

  const _fjsp_v4r8 nly  = _fjsp_broadcastload_v4i4(&NLy);
  const _fjsp_v4r8 nlz  = _fjsp_broadcastload_v4i4(&NLz);
  const _fjsp_v4r8 nlyz = _fjsp_mul_v4i4(nly, nlz);

#ifdef ARTED_DOMAIN_POWER_OF_TWO
  const _fjsp_v4r8 nlx = _fjsp_broadcastload_v4i4(&NLx);

  const int MX = NLx - 1;
  const int MY = NLy - 1;
  const _fjsp_v4r8 mx = _fjsp_broadcastload_v4i4(&MX);
  const _fjsp_v4r8 my = _fjsp_broadcastload_v4i4(&MY);

  const _fjsp_v4r8 ptbl  = _fjsp_set_v4i4(4, 3, 2, 1);
  const _fjsp_v4r8 yptbl = _fjsp_add_v4i4(nly, ptbl);
  const _fjsp_v4r8 xptbl = _fjsp_add_v4i4(nlx, ptbl);
#endif

  int idx[4], idy[4];

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
    _fjsp_v4r8 txp = _fjsp_sub_v4i4(tix, _fjsp_and_v4i4(_fjsp_add_v4i4(tix, xptbl), mx));
#else
    _fjsp_v4r8 txp = _fjsp_sub_v4i4(tix, _fjsp_load_v4i4(modx + (ix + 1 + NLx)));
#endif
               txp = _fjsp_mul_v4i4(nlyz, txp);
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    _fjsp_v4r8 tiy = _fjsp_broadcastload_v4i4(&iy);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
    _fjsp_v4r8 typ = _fjsp_sub_v4i4(tiy, _fjsp_and_v4i4(_fjsp_add_v4i4(tiy, yptbl), my));
#else
    _fjsp_v4r8 typ = _fjsp_sub_v4i4(tiy, _fjsp_load_v4i4(mody + (iy + 1 + NLy)));
#endif
               typ = _fjsp_mul_v4i4(nlz, typ);

    e = &E[ix][iy][0];

    const double complex *ee = e + 128;
    asm("prefetch [%0],20": :"r"(ee));

    for(iz = 0 ; iz < NLz ; iz += 2)
    {
      _fjsp_v4r8 tiz = _fjsp_broadcastload_v4i4(&iz);

      _fjsp_store_v4i4(idy, _fjsp_sub_v4i4(tiz, typ));
      _fjsp_store_v4i4(idx, _fjsp_sub_v4i4(tiz, txp));

#define STENCIL_CALC(PP,CC,TT) \
      TT = _fjsp_madd_v4r8(CC, PP, TT);

#define STENCIL_LOAD(IDP,PP) \
      PP = _fjsp_load_v4r8((double *) (e + IDP));

      _fjsp_v4r8 ez, w;
      _fjsp_v4r8 z0, z1, z2, z3;
      _fjsp_v4r8 y0, y1, y2, y3;
      _fjsp_v4r8 x0, x1, x2, x3;
      _fjsp_v4r8 vx, vy, vz;

      // conj(e[iz])
      ez = _fjsp_load_v4r8((double *)(e + iz));
      w  = _fjsp_xor_v4r8(ez, INV);

#ifdef ARTED_DOMAIN_POWER_OF_TWO
      z1 = _fjsp_load_v4r8((double *)(e + ((iz + 2 + NLz) & (NLz - 1))));
      z3 = _fjsp_load_v4r8((double *)(e + ((iz + 4 + NLz) & (NLz - 1))));
#else
      z1 = _fjsp_load_v4r8((double *)(e + modz[iz + 2 + NLz]));
      z3 = _fjsp_load_v4r8((double *)(e + modz[iz + 4 + NLz]));
#endif
      asm("### MAGIC"); /* !!! don't remove !!! */

      z0 = _fjsp_ecsl_v4(ez, z1, 2);
      z2 = _fjsp_ecsl_v4(z1, z3, 2);

      STENCIL_LOAD(idy[0], y0)
      STENCIL_LOAD(idy[1], y1)
      STENCIL_LOAD(idy[2], y2)
      STENCIL_LOAD(idy[3], y3)

      STENCIL_LOAD(idx[0], x0)
      STENCIL_LOAD(idx[1], x1)
      STENCIL_LOAD(idx[2], x2)
      STENCIL_LOAD(idx[3], x3)

      /* z-dimension (unit stride) */
      vz = _fjsp_setzero_v4r8();
      STENCIL_CALC(z3, c[11], vz);
      STENCIL_CALC(z1, c[ 9], vz);
      STENCIL_CALC(z0, c[ 8], vz);
      STENCIL_CALC(z2, c[10], vz);
      tt[2] = _fjsp_madd_v4r8(w, _fjsp_eperm_v4(vz, sel), tt[2]);

      /* y-dimension (NLz stride) */
      vy = _fjsp_setzero_v4r8();
      STENCIL_CALC(y0, c[4], vy);
      STENCIL_CALC(y1, c[5], vy);
      STENCIL_CALC(y2, c[6], vy);
      STENCIL_CALC(y3, c[7], vy);
      tt[1] = _fjsp_madd_v4r8(w, _fjsp_eperm_v4(vy, sel), tt[1]);

      /* x-dimension (NLy*NLz stride) */
      vx = _fjsp_setzero_v4r8();
      STENCIL_CALC(x0, c[0], vx);
      STENCIL_CALC(x1, c[1], vx);
      STENCIL_CALC(x2, c[2], vx);
      STENCIL_CALC(x3, c[3], vx);
      tt[0] = _fjsp_madd_v4r8(w, _fjsp_eperm_v4(vx, sel), tt[0]);
    } /* NLz */
  } /* NLy */
  } /* NLx */

  const _fjsp_v4r8 two = _fjsp_set_v4r8(2, 2, 2, 2);

  tt[0] = _fjsp_mul_v4r8(tt[0], two);
  tt[1] = _fjsp_mul_v4r8(tt[1], two);
  tt[2] = _fjsp_mul_v4r8(tt[2], two);

  double rf[4], rg[4], rh[4];

  _fjsp_store_v4r8(rf, tt[0]);
  _fjsp_store_v4r8(rg, tt[1]);
  _fjsp_store_v4r8(rh, tt[2]);

  *F = rf[0] + rf[1] + rf[2] + rf[3];
  *G = rg[0] + rg[1] + rg[2] + rg[3];
  *H = rh[0] + rh[1] + rh[2] + rh[3];

  _fjsp_simd_mode_2();
}
