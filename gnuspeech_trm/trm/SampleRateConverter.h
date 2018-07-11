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

#ifndef TRM_TRM_SAMPLE_RATE_CONVERTER_H
#define TRM_TRM_SAMPLE_RATE_CONVERTER_H

#include <vector>

/* Forward declaration */
struct trm_resampler;

namespace GS {
namespace TRM {
/**
 * Class used to resample from the frequency the tube update is performed
 * at to the desired output rate.
 */
class SampleRateConverter {
private:
	void *m_mem;
	struct trm_resampler *m_resampler;
	double m_maximum_sample_value;
	long m_number_samples;
	std::vector<float> &m_output_data;

public:
	SampleRateConverter(int sample_rate, float output_rate,
	                    std::vector<float> &output_data);
	~SampleRateConverter();

	void reset();
	void dataFill(float data);

	double maximumSampleValue() const { return m_maximum_sample_value; }
	long numberSamples() const { return m_number_samples; }
};

}  // namespace TRM
}  // namespace GS

#endif /* TRM_TRM_SAMPLE_RATE_CONVERTER_H */
