/* Copyright (C) 2007-2008 Jean-Marc Valin
 * Copyright (C) 2008 Thorvald Natvig
 */
/**
   @file resample_sse.h
   @brief Resampler functions (SSE version)
*/
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <xmmintrin.h>

#define OVERRIDE_INNER_PRODUCT_SINGLE
static inline float inner_product_single(const float *a, const float *b,
                                         unsigned int len)
{
	float ret;
	__m128 sum = _mm_setzero_ps();
	for (unsigned int i = 0; i < len; i += 8) {
		sum = _mm_add_ps(sum,
		                 _mm_mul_ps(_mm_loadu_ps(a + i), _mm_loadu_ps(b + i)));
		sum = _mm_add_ps(
		    sum, _mm_mul_ps(_mm_loadu_ps(a + i + 4), _mm_loadu_ps(b + i + 4)));
	}
	sum = _mm_add_ps(sum, _mm_movehl_ps(sum, sum));
	sum = _mm_add_ss(sum, _mm_shuffle_ps(sum, sum, 0x55));
	_mm_store_ss(&ret, sum);
	return ret;
}

#define OVERRIDE_INTERPOLATE_PRODUCT_SINGLE
static inline float interpolate_product_single(const float *a, const float *b,
                                               unsigned int len,
                                               const uint32_t oversample,
                                               float *frac)
{
	float ret;
	__m128 sum = _mm_setzero_ps();
	__m128 f = _mm_loadu_ps(frac);
	for (unsigned int i = 0; i < len; i += 2) {
		sum = _mm_add_ps(sum, _mm_mul_ps(_mm_load1_ps(a + i),
		                                 _mm_loadu_ps(b + i * oversample)));
		sum =
		    _mm_add_ps(sum, _mm_mul_ps(_mm_load1_ps(a + i + 1),
		                               _mm_loadu_ps(b + (i + 1) * oversample)));
	}
	sum = _mm_mul_ps(f, sum);
	sum = _mm_add_ps(sum, _mm_movehl_ps(sum, sum));
	sum = _mm_add_ss(sum, _mm_shuffle_ps(sum, sum, 0x55));
	_mm_store_ss(&ret, sum);
	return ret;
}

