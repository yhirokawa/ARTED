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
#include "./interop.h"

int NL, NLx, NLy, NLz;
int PNLx, PNLy, PNLz;

int BX, BY;

int const* modx;
int const* mody;
int const* modz;

void set_domain_size_(int * NL_, int * NLx_, int * NLy_, int * NLz_) {
  NL  = *NL_;
  NLx = *NLx_;
  NLy = *NLy_;
  NLz = *NLz_;
}

void set_padding_domain_size_(int * PNLx_, int * PNLy_, int * PNLz_) {
  PNLx = *PNLx_;
  PNLy = *PNLy_;
  PNLz = *PNLz_;
}

void set_blocking_factor_(int * BX_, int * BY_) {
  BX = *BX_;
  BY = *BY_;
}

void set_mod_table_(int * modx_, int * mody_, int * modz_) {
  modx = modx_;
  mody = mody_;
  modz = modz_;
}

