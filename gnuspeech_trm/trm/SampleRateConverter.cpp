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

#include <cmath>
#include <cstdlib>

#include <resample/resample.h>
#include <trm/SampleRateConverter.h>

namespace GS {
namespace TRM {

/******************************************************************************
 * Class SampleRateConverter                                                  *
 ******************************************************************************/

static const int RESAMPLER_QUALITY = 3; /* Do not use a high quality setting */

SampleRateConverter::SampleRateConverter(int sample_rate, float output_rate,
                                         std::vector<float> &output_data)
    : m_output_data(output_data)
{
	m_mem = malloc(
	    trm_resampler_size(1, sample_rate, output_rate, RESAMPLER_QUALITY));
	m_resampler = trm_resampler_init(m_mem, 1, sample_rate, output_rate,
	                                 RESAMPLER_QUALITY);
	reset();
}

SampleRateConverter::~SampleRateConverter() { free(m_mem); }

void SampleRateConverter::reset()
{
	m_maximum_sample_value = 0.0002;
	m_number_samples = 0;
	trm_resampler_reset_mem(m_resampler);
	trm_resampler_skip_zeros(m_resampler);
}

void SampleRateConverter::dataFill(float data)
{
	float out[16];
	uint32_t in_len = 1, out_len = 16;
	trm_resampler_process_interleaved_float(m_resampler, &data, &in_len, out,
	                                        &out_len);
	for (uint32_t i = 0; i < out_len; i++) {
		/* Track the maximum value */
		if (std::fabs(out[i]) > m_maximum_sample_value) {
			m_maximum_sample_value = fabs(out[i]);
		}

		/* Push the output sample onto the target list */
		m_output_data.emplace_back(out[i]);
		m_number_samples++;
	}
}
}  // namespace TRM
}  // namespace GS

