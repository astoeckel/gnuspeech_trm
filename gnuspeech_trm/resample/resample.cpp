/*
 *  gnuspeech_trm -- Standalone gnuspeech articulatory synthesiser
 *  Copyright (C) 2018  Andreas Stöckel
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*
 *  This code is a modified version of the resampler included with speexdsp,
 *  licensed under the 3-Clause BSD License. You can find the original code at
 *
 *      https://git.xiph.org/?p=speexdsp.git;a=summary
 *
 *  --------------------------------------------------------------------------
 *
 *  Copyright (C) 2007-2008 Jean-Marc Valin
 *  Copyright (C) 2008      Thorvald Natvig
 *
 *  File: resample.c
 *  Arbitrary resampling code
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. The name of the author may not be used to endorse or promote products
 *  derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 *  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 *  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *  The design goals of this code are:
 *     - Very fast algorithm
 *     - SIMD-friendly algorithm
 *     - Low memory requirement
 *     - Good *perceptual* quality (and not best SNR)
 *
 *  Warning: This resampler is relatively new. Although I think I got rid of
 *  all the major bugs and I don't expect the API to change anymore, there
 *  may be something I've missed. So use with caution.
 *
 *  This algorithm is based on this original resampling algorithm:
 *  Smith, Julius O. Digital Audio Resampling Home Page
 *  Center for Computer Research in Music and Acoustics (CCRMA),
 *  Stanford University, 2007.
 *  Web published at http://ccrma.stanford.edu/~jos/resample/.
 *
 *  There is one main difference, though. This resampler uses cubic
 *  interpolation instead of linear interpolation in the above paper. This
 *  makes the table much smaller and makes it possible to compute that table
 *  on a per-stream basis. In turn, being able to tweak the table for each
 *  stream makes it possible to both reduce complexity on simple ratios
 *  (e.g. 2/3), and get rid of the rounding operations in the inner loop.
 *  The latter both reduces CPU time and makes the algorithm more SIMD-friendly.
 */

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <resample/resample.h>
#ifdef __SSE__ /* Use SSE if available on the target architecture */
#include <resample/resample_sse.h>
#endif /* __SSE__ */

/******************************************************************************
 * MEMORY ALIGNMENT UTILITIES                                                 *
 ******************************************************************************/

/**
 * Memory alignment for pointers internally used by Stanchion. Aligning memory
 * and telling the compiler about it allows the compiler to perform better
 * optimization. Furthermore, some platforms (WASM) do not allow unaligned
 * memory access.
 */
#define ALIGN 16

/**
 * Macro telling the compiler that P is aligned with the alignment defined
 * above.
 */
#define ASSUME_ALIGNED(P) P
#ifdef __GNUC__
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || (__GNUC__ > 4)
#undef ASSUME_ALIGNED
#define ASSUME_ALIGNED(P) (__builtin_assume_aligned(P, ALIGN))
#endif /* (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || (__GNUC__ > 4) */
#endif /* __GNUC__ */

/**
 * Forces a pointer to have the alignment defined above.
 */
#define ALIGN_ADDR(P) \
	(ASSUME_ALIGNED((void *)(((uintptr_t)(P) + ALIGN - 1) & (~(ALIGN - 1)))))

/**
 * Function used internally to compute the total size of a datastructure
 * consisting of multiple substructure. Calling this function updates the size
 * of the outer datastructure by adding a substructure of size n_bytes. Assumes
 * that the beginning of the substructure should be aligned.
 *
 * @param size is a pointer at the variable holding the size of the
 * datastructure.
 * @param n_bytes size of the sub-structure that should be added.
 * @return zero if there was an overflow, one otherwise.
 */
static inline bool mem_update_size(uint32_t *size, uint32_t n_bytes) {
	const uint32_t new_size = ((*size + ALIGN - 1) & (~(ALIGN - 1))) + n_bytes;
	if (new_size < *size) {
		return false; /* error, there has been an overflow */
	}
	*size = new_size;
	return true; /* success */
}

/**
 * Computes the aligned pointer pointing at the substructure of the given size.
 *
 * @param mem pointer at the variable holding the pointer at the current
 * pointer. The pointer is advanced by the given size after the return value is
 * computed.
 * @param size is the size of the substructure for which the pointer should be
 * returned.
 * @return an aligned pointer pointing at the beginning of the substructure.
 */
static inline void* mem_align(void **mem, uint32_t size) {
	void *res = ALIGN_ADDR(*mem);
	*mem = (void*)((uintptr_t)res + size);
	return res;
}

/******************************************************************************
 * TYPE DEFINITIONS                                                           *
 ******************************************************************************/

struct trm_resampler_func_def {
	const double *table;
	int oversample;
};

struct trm_quality_mapping {
	int base_length;
	int oversample;
	float downsample_bandwidth;
	float upsample_bandwidth;
	const struct trm_resampler_func_def *window_func;
};

/**
 * Datastructure holding the private data of the resampler instance. Note that
 * the arrays inside this structure are placed in memory after this structure.
 */
struct trm_resampler {
	/**
	 * Input sample rate in samples per second.
	 */
	uint32_t in_rate;

	/**
	 * Output sample rate in samples per second.
	 */
	uint32_t out_rate;

	/**
	 * Numerator of the normalized input/output ratio.
	 */
	uint32_t num_rate;

	/**
	 * Denumerator of the normalized input/output ratio.
	 */
	uint32_t den_rate;

	/**
	 * Quality 0 <= q <= 10.
	 */
	uint8_t quality;

	/**
	 * Number of channels that are being processed by the resampler. Determines
	 * the amount of memory that is allocated for the sample buffer.
	 */
	uint32_t nb_channels;

	/**
	 * Length of the filter in samples. The length depends on the input/output
	 * ratio as well as the selected quality.
	 */
	uint32_t filt_len;

	/**
	 * Length of the Sinc function table.
	 */
	uint32_t sinc_table_length;

	/**
	 * Size of the sample buffer per channel.
	 */
	uint32_t mem_alloc_size;

	/**
	 * Integer part of the ratio between num_rate and den_rate. Used to
	 * increment the input sample pointer.
	 */
	uint32_t int_advance;

	/**
	 * Fractional part of the ratio between num_rate and den_rate. Used to
	 * increment the fractional part of the input sample pointer.
	 */
	uint32_t frac_advance;

	/**
	 * Cutoff frequency of the filter.
	 */
	float cutoff;

	/**
	 * Number of taps on the input between two samples. (???)
	 */
	uint32_t oversample;

	/**
	 * Index of the last sample that was processed for each channel.
	 */
	uint32_t *last_sample;

	/**
	 * Fractional index of the last sample that was processed for each channel.
	 */
	uint32_t *samp_frac_num;

	/**
	 * Pointer at the sample buffer.
	 */
	float *mem;

	/**
	 * Raw filter; depending on the "use_direct" flag, this may either directly
	 * correspond to the filter, or the filter may be interpolated from this
	 * table.
	 */
	float *sinc_table;

	/**
	 * If true, does not interpolate the sinc function table (the filter is the
	 * sinc_table), if false, interpolates the sinc function table (the filter
	 * is interpolated from the sinc_table).
	 */
	bool use_direct;
};

/******************************************************************************
 * CONSTANTS AND TABLES                                                       *
 ******************************************************************************/

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/**
 * Additional padding added to the sample buffer for each channel.
 *
 * TODO: Figure out what this is actually used for.
 */
#define BUFFER_SIZE 160

static const double kaiser12_table[68] = {
    0.99859849, 1.00000000, 0.99859849, 0.99440475, 0.98745105, 0.97779076,
    0.96549770, 0.95066529, 0.93340547, 0.91384741, 0.89213598, 0.86843014,
    0.84290116, 0.81573067, 0.78710866, 0.75723148, 0.72629970, 0.69451601,
    0.66208321, 0.62920216, 0.59606986, 0.56287762, 0.52980938, 0.49704014,
    0.46473455, 0.43304576, 0.40211431, 0.37206735, 0.34301800, 0.31506490,
    0.28829195, 0.26276832, 0.23854851, 0.21567274, 0.19416736, 0.17404546,
    0.15530766, 0.13794294, 0.12192957, 0.10723616, 0.09382272, 0.08164178,
    0.07063950, 0.06075685, 0.05193064, 0.04409466, 0.03718069, 0.03111947,
    0.02584161, 0.02127838, 0.01736250, 0.01402878, 0.01121463, 0.00886058,
    0.00691064, 0.00531256, 0.00401805, 0.00298291, 0.00216702, 0.00153438,
    0.00105297, 0.00069463, 0.00043489, 0.00025272, 0.00013031, 0.0000527734,
    0.00001000, 0.00000000};

static const double kaiser10_table[36] = {
    0.99537781, 1.00000000, 0.99537781, 0.98162644, 0.95908712, 0.92831446,
    0.89005583, 0.84522401, 0.79486424, 0.74011713, 0.68217934, 0.62226347,
    0.56155915, 0.50119680, 0.44221549, 0.38553619, 0.33194107, 0.28205962,
    0.23636152, 0.19515633, 0.15859932, 0.12670280, 0.09935205, 0.07632451,
    0.05731132, 0.04193980, 0.02979584, 0.02044510, 0.01345224, 0.00839739,
    0.00488951, 0.00257636, 0.00115101, 0.00035515, 0.00000000, 0.00000000};

static const double kaiser8_table[36] = {
    0.99635258, 1.00000000, 0.99635258, 0.98548012, 0.96759014, 0.94302200,
    0.91223751, 0.87580811, 0.83439927, 0.78875245, 0.73966538, 0.68797126,
    0.63451750, 0.58014482, 0.52566725, 0.47185369, 0.41941150, 0.36897272,
    0.32108304, 0.27619388, 0.23465776, 0.19672670, 0.16255380, 0.13219758,
    0.10562887, 0.08273982, 0.06335451, 0.04724088, 0.03412321, 0.02369490,
    0.01563093, 0.00959968, 0.00527363, 0.00233883, 0.00050000, 0.00000000};

static const double kaiser6_table[36] = {
    0.99733006, 1.00000000, 0.99733006, 0.98935595, 0.97618418, 0.95799003,
    0.93501423, 0.90755855, 0.87598009, 0.84068475, 0.80211977, 0.76076565,
    0.71712752, 0.67172623, 0.62508937, 0.57774224, 0.53019925, 0.48295561,
    0.43647969, 0.39120616, 0.34752997, 0.30580127, 0.26632152, 0.22934058,
    0.19505503, 0.16360756, 0.13508755, 0.10953262, 0.08693120, 0.06722600,
    0.05031820, 0.03607231, 0.02432151, 0.01487334, 0.00752000, 0.00000000};

static const struct trm_resampler_func_def _KAISER12 = {kaiser12_table, 64};
#define KAISER12 (&_KAISER12)
static const struct trm_resampler_func_def _KAISER10 = {kaiser10_table, 32};
#define KAISER10 (&_KAISER10)
static const struct trm_resampler_func_def _KAISER8 = {kaiser8_table, 32};
#define KAISER8 (&_KAISER8)
static const struct trm_resampler_func_def _KAISER6 = {kaiser6_table, 32};
#define KAISER6 (&_KAISER6)

/* This table maps conversion quality to internal parameters. There are two
   reasons that explain why the up-sampling bandwidth is larger than the
   down-sampling bandwidth:
   1) When up-sampling, we can assume that the spectrum is already attenuated
      close to the Nyquist rate (from an A/D or a previous resampling filter)
   2) Any aliasing that occurs very close to the Nyquist rate will be masked
      by the sinusoids/noise just below the Nyquist rate (guaranteed only for
      up-sampling).
*/
static const struct trm_quality_mapping quality_map[11] = {
    {8, 4, 0.830f, 0.860f, KAISER6},  /* Q0 */
    {16, 4, 0.850f, 0.880f, KAISER6}, /* Q1 */
    {32, 4, 0.882f, 0.910f, KAISER6},
    /* Q2 */ /* 82.3% cutoff ( ~60 dB stop) 6  */
    {48, 8, 0.895f, 0.917f, KAISER8},
    /* Q3 */ /* 84.9% cutoff ( ~80 dB stop) 8  */
    {64, 8, 0.921f, 0.940f, KAISER8},
    /* Q4 */ /* 88.7% cutoff ( ~80 dB stop) 8  */
    {80, 16, 0.922f, 0.940f, KAISER10},
    /* Q5 */ /* 89.1% cutoff (~100 dB stop) 10 */
    {96, 16, 0.940f, 0.945f, KAISER10},
    /* Q6 */ /* 91.5% cutoff (~100 dB stop) 10 */
    {128, 16, 0.950f, 0.950f, KAISER10},
    /* Q7 */ /* 93.1% cutoff (~100 dB stop) 10 */
    {160, 16, 0.960f, 0.960f, KAISER10},
    /* Q8 */ /* 94.5% cutoff (~100 dB stop) 10 */
    {192, 32, 0.968f, 0.968f, KAISER12},
    /* Q9 */ /* 95.5% cutoff (~100 dB stop) 10 */
    {256, 32, 0.975f, 0.975f, KAISER12},
    /* Q10 */ /* 96.6% cutoff (~100 dB stop) 10 */
};

/******************************************************************************
 * MATH UTILITIES                                                             *
 ******************************************************************************/

/*8,24,40,56,80,104,128,160,200,256,320*/
static double _compute_func(float x, const struct trm_resampler_func_def *func) {
	float y, frac;
	double interp[4];
	int ind;
	y = x * func->oversample;
	ind = (int)floor(y);
	frac = (y - ind);
	/* CSE with handle the repeated powers */
	interp[3] = -0.1666666667 * frac + 0.1666666667 * (frac * frac * frac);
	interp[2] = frac + 0.5 * (frac * frac) - 0.5 * (frac * frac * frac);
	/*interp[2] = 1.f - 0.5f*frac - frac*frac + 0.5f*frac*frac*frac;*/
	interp[0] = -0.3333333333 * frac + 0.5 * (frac * frac) -
	            0.1666666667 * (frac * frac * frac);
	/* Just to make sure we don't have rounding problems */
	interp[1] = 1.f - interp[3] - interp[2] - interp[0];

	/*sum = frac*accum[1] + (1-frac)*accum[2];*/
	return interp[0] * func->table[ind] + interp[1] * func->table[ind + 1] +
	       interp[2] * func->table[ind + 2] + interp[3] * func->table[ind + 3];
}

/* The slow way of computing a sinc for the table. Should improve that some day
 */
static float _sinc(float cutoff, float x, int N,
                   const struct trm_resampler_func_def *window_func) {
	float xx = x * cutoff;
	if (fabs(x) < 1e-6)
		return cutoff;
	else if (fabs(x) > .5 * N)
		return 0;
	/*FIXME: Can it really be any slower than this? */
	return cutoff * sin(M_PI * xx) / (M_PI * xx) *
	       _compute_func(fabs(2. * x / N), window_func);
}

static void _cubic_coef(float frac, float interp[4]) {
	/* Compute interpolation coefficients. I'm not sure whether this corresponds
	to cubic interpolation but I know it's MMSE-optimal on a sinc */
	interp[0] = -0.16667f * frac + 0.16667f * frac * frac * frac;
	interp[1] = frac + 0.5f * frac * frac - 0.5f * frac * frac * frac;
	/*interp[2] = 1.f - 0.5f*frac - frac*frac + 0.5f*frac*frac*frac;*/
	interp[3] =
	    -0.33333f * frac + 0.5f * frac * frac - 0.16667f * frac * frac * frac;
	/* Just to make sure we don't have rounding problems */
	interp[2] = 1. - interp[0] - interp[1] - interp[3];
}

/**
 * Computes the greatest common divisor between a, b.
 */
static inline uint32_t _gcd(uint32_t a, uint32_t b) {
	while (b != 0) {
		uint32_t temp = a;
		a = b;
		b = temp % b;
	}
	return a;
}

/**
 * Multiplies the given value with the multiplicator and divides by the given
 * divisor; checks for overflows.
 *
 * @param result is a pointer at the variable in which the result should be
 * stored.
 * @param value is the value that should be multiplied with the given
 * multiplier.
 * @param div is the divisor that the result of the multiplication should be
 * divided by.
 */
static bool _muldiv(uint32_t *result, uint32_t value, uint32_t mul,
                    uint32_t div) {
	const uint64_t res64 = (value * mul) / div;
	*result = res64;
	return res64 <= UINT32_MAX;
}

/******************************************************************************
 * PRIVATE IMPLEMENTATION DETAILS                                             *
 ******************************************************************************/

static uint32_t trm_resampler_direct_single(trm_resampler_t *st,
                                           uint32_t channel_index,
                                           const float *in, uint32_t *in_len,
                                           float *out, uint32_t *out_len,
                                           uint32_t out_stride) {
	const uint32_t N = st->filt_len;
	uint32_t out_sample = 0;
	uint32_t last_sample = st->last_sample[channel_index];
	uint32_t samp_frac_num = st->samp_frac_num[channel_index];
	const float *sinc_table = st->sinc_table;
	const uint32_t int_advance = st->int_advance;
	const uint32_t frac_advance = st->frac_advance;
	const uint32_t den_rate = st->den_rate;
	float sum;

	while (!(last_sample >= *in_len || out_sample >= *out_len)) {
		const float *sinct = &sinc_table[samp_frac_num * N];
		const float *iptr = &in[last_sample];

#ifndef OVERRIDE_INNER_PRODUCT_SINGLE
		sum = 0;
		for (uint32_t j = 0; j < N; j++) {
			sum += sinct[j] * iptr[j];
		}
#else
		sum = inner_product_single(sinct, iptr, N);
#endif

		out[out_stride * out_sample++] = sum;
		last_sample += int_advance;
		samp_frac_num += frac_advance;
		if (samp_frac_num >= den_rate) {
			samp_frac_num -= den_rate;
			last_sample++;
		}
	}

	st->last_sample[channel_index] = last_sample;
	st->samp_frac_num[channel_index] = samp_frac_num;
	return out_sample;
}

static uint32_t trm_resampler_interpolate_single(
    trm_resampler_t *st, uint32_t channel_index, const float *in,
    uint32_t *in_len, float *out, uint32_t *out_len, uint32_t out_stride) {
	const uint32_t N = st->filt_len;
	uint32_t out_sample = 0;
	uint32_t last_sample = st->last_sample[channel_index];
	uint32_t samp_frac_num = st->samp_frac_num[channel_index];
	const uint32_t int_advance = st->int_advance;
	const uint32_t frac_advance = st->frac_advance;
	const uint32_t den_rate = st->den_rate;
	float sum;

	while (!(last_sample >= *in_len || out_sample >= *out_len)) {
		const float *iptr = &in[last_sample];

		const int offset = samp_frac_num * st->oversample / st->den_rate;
		const float frac =
		    ((float)((samp_frac_num * st->oversample) % st->den_rate)) /
		    st->den_rate;
		float interp[4];

#ifndef OVERRIDE_INTERPOLATE_PRODUCT_SINGLE
		float accum[4] = {0, 0, 0, 0};
		for (uint32_t j = 0; j < N; j++) {
			const float curr_in = iptr[j];
			accum[0] +=
			    curr_in *
			    st->sinc_table[4 + (j + 1) * st->oversample - offset - 2];
			accum[1] +=
			    curr_in *
			    st->sinc_table[4 + (j + 1) * st->oversample - offset - 1];
			accum[2] +=
			    curr_in * st->sinc_table[4 + (j + 1) * st->oversample - offset];
			accum[3] +=
			    curr_in *
			    st->sinc_table[4 + (j + 1) * st->oversample - offset + 1];
		}

		_cubic_coef(frac, interp);
		sum = interp[0] * accum[0] + interp[1] * accum[1] +
		      interp[2] * accum[2] + interp[3] * accum[3];
#else
		_cubic_coef(frac, interp);
		sum = interpolate_product_single(
		    iptr, st->sinc_table + st->oversample + 4 - offset - 2, N,
		    st->oversample, interp);
#endif

		out[out_stride * out_sample++] = sum;
		last_sample += int_advance;
		samp_frac_num += frac_advance;
		if (samp_frac_num >= den_rate) {
			samp_frac_num -= den_rate;
			last_sample++;
		}
	}

	st->last_sample[channel_index] = last_sample;
	st->samp_frac_num[channel_index] = samp_frac_num;
	return out_sample;
}

/**
 * Computes the filter lengths and depdentend values from the input/output rate,
 * number of channels, and quality stored in the given instance.
 */
static bool trm_resampler_compute_filter_size(trm_resampler_t *st) {
	/* Normalize the sample rates */
	uint32_t divisor = _gcd(st->num_rate, st->den_rate);
	st->num_rate = st->num_rate / divisor;
	st->den_rate = st->den_rate / divisor;

	/* Compute the integer and fractional per of the ration between num_rate
	   and den_rate */
	st->int_advance = st->num_rate / st->den_rate;
	st->frac_advance = st->num_rate % st->den_rate;

	/* Compute the oversampling factor and the filter length */
	st->oversample = quality_map[st->quality].oversample;
	st->filt_len = quality_map[st->quality].base_length;
	if (st->num_rate > st->den_rate) {
		/* Down-sampling */
		st->cutoff = quality_map[st->quality].downsample_bandwidth *
		             st->den_rate / st->num_rate;

		/* Compute the filter length */
		if (!_muldiv(&st->filt_len, st->filt_len, st->num_rate, st->den_rate)) {
			return false;
		}

		/* Round up to make sure we have a multiple of 8 for SSE */
		st->filt_len = ((st->filt_len - 1) & (~0x7)) + 8;
		if (2 * st->den_rate < st->num_rate) {
			st->oversample >>= 1;
		}
		if (4 * st->den_rate < st->num_rate) {
			st->oversample >>= 1;
		}
		if (8 * st->den_rate < st->num_rate) {
			st->oversample >>= 1;
		}
		if (16 * st->den_rate < st->num_rate) {
			st->oversample >>= 1;
		}
		if (st->oversample < 1) {
			st->oversample = 1;
		}
	} else {
		/* Up-sampling */
		st->cutoff = quality_map[st->quality].upsample_bandwidth;
	}

	/* Choose the resampling type that requires the least amount of memory */
	st->use_direct =
	    st->filt_len * st->den_rate <= st->filt_len * st->oversample + 8 &&
	    UINT32_MAX / sizeof(float) / st->den_rate >= st->filt_len;
	if (st->use_direct) {
		st->sinc_table_length = st->den_rate * st->filt_len;
	} else {
		/* Make sure the required amount of memory is not larger than 2GiB */
		if ((UINT32_MAX / sizeof(float) - 8) / st->oversample < st->filt_len) {
			return false;
		}
		st->sinc_table_length = st->oversample * st->filt_len + 8;
	}

	/* Compute the sample buffer length */
	st->mem_alloc_size =
	    (st->filt_len - 1 + BUFFER_SIZE + ALIGN - 1) & ~(ALIGN - 1);

	return true;
}

static void trm_resampler_compute_filter(trm_resampler_t *st) {
	if (st->use_direct) {
		for (uint32_t i = 0; i < st->den_rate; i++) {
			for (uint32_t j = 0; j < st->filt_len; j++) {
				st->sinc_table[i * st->filt_len + j] =
				    _sinc(st->cutoff,
				          ((j - (int32_t)st->filt_len / 2 + 1) -
				           ((float)i) / st->den_rate),
				          st->filt_len, quality_map[st->quality].window_func);
			}
		}
	} else {
		for (int32_t i = -4; i < (int32_t)(st->oversample * st->filt_len + 4);
		     i++) {
			st->sinc_table[i + 4] = _sinc(
			    st->cutoff, (i / (float)st->oversample - st->filt_len / 2),
			    st->filt_len, quality_map[st->quality].window_func);
		}
	}
}

static trm_resampler_t *trm_resampler_init_frac(void *mem, uint32_t nb_channels,
                                              uint32_t ratio_num,
                                              uint32_t ratio_den,
                                              uint32_t in_rate,
                                              uint32_t out_rate, int quality) {
	/* Make sure the arguments are sane */
	if (!mem || nb_channels == 0 || ratio_num == 0 || ratio_den == 0 || quality > 10 ||
	    quality < 0) {
		return NULL;
	}

	/* Initialize the trm_resampler struct */
	trm_resampler_t *st =
	    (trm_resampler_t *)mem_align(&mem, sizeof(trm_resampler_t));
	st->in_rate = in_rate;
	st->out_rate = out_rate;
	st->num_rate = ratio_num;
	st->den_rate = ratio_den;
	st->quality = quality;
	st->nb_channels = nb_channels;

	/* Compute the filter lengths */
	trm_resampler_compute_filter_size(st);

	/* Per channel data */
	st->last_sample =
	    (uint32_t *)mem_align(&mem, nb_channels * sizeof(uint32_t));
	st->samp_frac_num =
	    (uint32_t *)mem_align(&mem, nb_channels * sizeof(uint32_t));

	/* Filter and sample buffer */
	st->sinc_table =
	    (float *)mem_align(&mem, st->sinc_table_length * sizeof(float));
	st->mem = (float *)mem_align(
	    &mem, nb_channels * st->mem_alloc_size * sizeof(float));

	/* Compute the filter itself */
	trm_resampler_compute_filter(st);

	/* Reset the memory */
	trm_resampler_reset_mem(st);

	return st;
}

static void trm_resampler_process_native(trm_resampler_t *st,
                                        uint32_t channel_index,
                                        uint32_t *in_len, float *out,
                                        uint32_t *out_len,
                                        uint32_t out_stride) {
	const uint32_t N = st->filt_len;
	uint32_t out_sample = 0;
	float *mem =
	    (float *)ASSUME_ALIGNED(st->mem + channel_index * st->mem_alloc_size);

	/* Call the right resampler */
	if (st->use_direct) {
		out_sample = trm_resampler_direct_single(st, channel_index, mem, in_len,
		                                        out, out_len, out_stride);
	} else {
		out_sample = trm_resampler_interpolate_single(
		    st, channel_index, mem, in_len, out, out_len, out_stride);
	}

	if (st->last_sample[channel_index] < *in_len) {
		*in_len = st->last_sample[channel_index];
	}
	*out_len = out_sample;
	st->last_sample[channel_index] -= *in_len;

	const uint32_t ilen = *in_len;
	for (uint32_t j = 0; j < N - 1; ++j) {
		mem[j] = mem[j + ilen];
	}
}

/******************************************************************************
 * PUBLIC API                                                                 *
 ******************************************************************************/

uint32_t trm_resampler_size(uint32_t nb_channels, uint32_t in_rate,
                           uint32_t out_rate, int quality) {
	uint32_t size = ALIGN - 1; /* Ensure space to align the entire structure */

	/* Make sure the arguments passed to the resampler are correct */
	if (nb_channels == 0 || in_rate == 0 || out_rate == 0 || quality > 10 ||
	    quality < 0) {
		return 0;
	}

	/* Allocate a temporary trm_resampler_t on the stack and fill in the given
	   parameters */
	trm_resampler_t st;
	st.in_rate = st.num_rate = in_rate;
	st.out_rate = st.den_rate = out_rate;
	st.quality = quality;

	/* Compute the size of the sinc filter and the buffer */
	if (!trm_resampler_compute_filter_size(&st)) {
		return 0;
	}

	/* Accumulate the sizes of the individual substructures while taking
	   alignment into account. */
	const bool ok =
	    mem_update_size(&size, sizeof(trm_resampler_t)) &&
	    mem_update_size(&size, nb_channels * sizeof(float)) &&
	    mem_update_size(&size, nb_channels * sizeof(float)) &&
	    mem_update_size(&size, nb_channels * sizeof(float)) &&
	    mem_update_size(&size, st.sinc_table_length * sizeof(float)) &&
	    mem_update_size(&size, nb_channels * st.mem_alloc_size * sizeof(float));
	return ok ? size : 0;
}

trm_resampler_t *trm_resampler_init(void *mem, uint32_t nb_channels,
                                  uint32_t in_rate, uint32_t out_rate,
                                  int quality) {
	return trm_resampler_init_frac(mem, nb_channels, in_rate, out_rate, in_rate,
	                              out_rate, quality);
}

void trm_resampler_process_float(trm_resampler_t *st, uint32_t channel_index,
                                const float *in, uint32_t *in_len,
                                uint32_t in_stride, float *out,
                                uint32_t *out_len, uint32_t out_stride) {
	uint32_t ilen = *in_len;
	uint32_t olen = *out_len;
	float *x =
	    (float *)ASSUME_ALIGNED(st->mem + channel_index * st->mem_alloc_size);
	const int filt_offs = st->filt_len - 1;
	const uint32_t xlen = st->mem_alloc_size - filt_offs;

	while (ilen && olen) {
		uint32_t ichunk = (ilen > xlen) ? xlen : ilen;
		uint32_t ochunk = olen;

		if (in) {
			for (uint32_t j = 0; j < ichunk; ++j) {
				x[j + filt_offs] = in[j * in_stride];
			}
		} else {
			for (uint32_t j = 0; j < ichunk; ++j) {
				x[j + filt_offs] = 0.0f;
			}
		}
		trm_resampler_process_native(st, channel_index, &ichunk, out, &ochunk,
		                            out_stride);
		ilen -= ichunk;
		olen -= ochunk;
		out += ochunk * out_stride;
		if (in) {
			in += ichunk * in_stride;
		}
	}
	*in_len -= ilen;
	*out_len -= olen;
}

void trm_resampler_process_interleaved_float(trm_resampler_t *st, const float *in,
                                            uint32_t *in_len, float *out,
                                            uint32_t *out_len) {
	uint32_t orig_out_len = *out_len, orig_in_len = *in_len;
	uint32_t in_stride = st->nb_channels, out_stride = st->nb_channels;
	for (uint32_t i = 0; i < st->nb_channels; i++) {
		*out_len = orig_out_len, *in_len = orig_in_len;
		if (in != NULL) {
			trm_resampler_process_float(st, i, in + i, in_len, in_stride,
			                           out + i, out_len, out_stride);
		} else {
			trm_resampler_process_float(st, i, NULL, in_len, in_stride, out + i,
			                           out_len, out_stride);
		}
	}
}

uint32_t trm_resampler_get_input_latency(trm_resampler_t *st) {
	return st->filt_len / 2;
}

uint32_t trm_resampler_get_output_latency(trm_resampler_t *st) {
	return ((st->filt_len / 2) * st->den_rate + (st->num_rate >> 1)) /
	       st->num_rate;
}

void trm_resampler_skip_zeros(trm_resampler_t *st) {
	for (uint32_t i = 0; i < st->nb_channels; i++) {
		st->last_sample[i] = st->filt_len / 2;
	}
}

void trm_resampler_reset_mem(trm_resampler_t *st) {
	for (uint32_t i = 0; i < st->nb_channels; i++) {
		st->last_sample[i] = 0;
		st->samp_frac_num[i] = 0;
	}
	for (uint32_t i = 0; i < st->nb_channels * st->mem_alloc_size; i++) {
		st->mem[i] = 0.0f;
	}
}

