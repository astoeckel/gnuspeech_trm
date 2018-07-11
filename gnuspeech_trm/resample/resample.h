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

/**
 * @file resample.h
 *
 * This file provides an interface to the Speex resampler. In contrast to the
 * orignal Speex resampler, the code has been modified to run on a fixed-size
 * memory region. There are no memory allocations, and no dependency on stdlib.
 * All code paths that are not used have been pruned. Furthermore, higher
 * quality settings no longer use double arithmetic, which reduces the code size
 * a little.
 *
 * @author Jean-Marc Valin, Thorvald Natvig, Andreas Stöckel
 */

#ifndef GNUSPEECH_TRM_RESAMPLE_RESAMPLE_H
#define GNUSPEECH_TRM_RESAMPLE_RESAMPLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Opaque struct representing a resampler instance.
 */
struct trm_resampler;

/** 
 * Typedef for the trm_resampler struct.
 */
typedef struct trm_resampler trm_resampler_t;

/**
 * Returns the minimum size of the memory area that should be reserved for a
 * resampler instance with the given properties.
 *
 * @param nb_channels Number of channels to be processed.
 * @param in_rate Input sampling rate (in Hz)
 * @param out_rate Output sampling rate (in Hz)
 * @param quality Resampling quality
 * @return minimum number of bytes to allocate for the resampler, or zero if an
 * error occurred, see the error written to err for more details.
 */
uint32_t trm_resampler_size(uint32_t nb_channels, uint32_t in_rate,
                           uint32_t out_rate, int quality);

/**
 * Create a new resampler with integer input and output rates.
 *
 * @param inst pointer at the memory region the instance should be placed at.
 * Must be at least as large the value returned by trm_resampler_size with the
 * same parameters. If NULL, this function returns NULL.
 * @param nb_channels Number of channels to be processed
 * @param in_rate Input sampling rate (integer number of Hz).
 * @param out_rate Output sampling rate (integer number of Hz).
 * @param quality Resampling quality between 0 and 10, where 0 has poor quality
 * and 10 has very high quality.
 * @param err pointer at an integer the error code is being written to.
 * @return a pointer within the memory region pointer passed to this function
 * or NULL if the given parameters are invalid. The returned value is aligned
 * must be passed to all resampler functions instead of directly passing mem.
 */
trm_resampler_t *trm_resampler_init(void *mem, uint32_t nb_channels,
                                  uint32_t in_rate, uint32_t out_rate,
                                  int quality);

/**
 * Resample a float array with the specified memory layout. The input and output
 * buffers must *not* overlap.
 *
 * @param st resampler state.
 * @param channel_index index of the channel to process for the multi-channel
 * base (0 otherwise)
 * @param in input buffer
 * @param in_len number of input samples in the input buffer. Returns the
 * number of samples processed
 * @param in_stride number of floats between two consecutive samples for the
 * same channel, should be one for a dense (planar) input array.
 * @param out output buffer
 * @param out_len size of the output buffer. Returns the number of samples
 * written
 * @param out_stride number of floats between two consecutive samples for the
 * same channel, should be one for a dense (planar) output array.
 */
void trm_resampler_process_float(trm_resampler_t *st, uint32_t channel_index,
                                const float *in, uint32_t *in_len,
                                uint32_t in_stride, float *out,
                                uint32_t *out_len, uint32_t out_stride);

/**
 * Resample an interleaved float array. The input and output buffers must *not*
 * overlap. This is essentially the same as trm_resampler_process_float, but with
 * input and output stride set to the number of channels.
 *
 * @param st Resampler state
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number
 * of samples processed. This is all per-channel.
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples
 * written. This is all per-channel.
 */
void trm_resampler_process_interleaved_float(trm_resampler_t *st, const float *in,
                                            uint32_t *in_len, float *out,
                                            uint32_t *out_len);

/**
 * Get the latency introduced by the resampler measured in input samples.
 *
 * @param st Resampler state
 */
uint32_t trm_resampler_get_input_latency(trm_resampler_t *st);

/**
 * Get the latency introduced by the resampler measured in output samples.
 *
 * @param st Resampler state
 */
uint32_t trm_resampler_get_output_latency(trm_resampler_t *st);

/**
 * Make sure that the first samples to go out of the resamplers don't have
 * leading zeros. This is only useful before starting to use a newly created
 * resampler. It is recommended to use that when resampling an audio file, as
 * it will generate a file with the same length. For real-time processing,
 * it is probably easier not to use this call (so that the output duration
 * is the same for the first frame).
 *
 * @param st Resampler state
 */
void trm_resampler_skip_zeros(trm_resampler_t *st);

/**
 * Reset a resampler so a new (unrelated) stream can be processed.
 *
 * @param st Resampler state
 */
void trm_resampler_reset_mem(trm_resampler_t *st);

#ifdef __cplusplus
}
#endif

#endif /* GNUSPEECH_TRM_RESAMPLE_RESAMPLE_H */
