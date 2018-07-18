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

#include <iostream>

#include "gnuspeech_trm.h"
#include "trm/Tube.h"

/******************************************************************************
 * C++ Private Implementation                                                 *
 ******************************************************************************/

namespace {
/**
 * Structure holding the configuration data for a Tube object. The external API
 * allows to read/write the parameters in this struct
 */
struct TrmConfig {
	/**
	 * Output sample rate in samples per second.
	 */
	double output_rate = 44100.0;

	/**
	 * Update the parameters this many times per second.
	 */
	double control_rate = 250.0;

	/**
	 * Calculate a moving average over the input values with the given period in
	 * seconds.
	 */
	double filter_period = 50e-3;

	/**
	 * Master volume (0 - 60 dB)
	 */
	double volume = 60.0;

	/**
	 * GS waveform type (0=PULSE, 1=SINE)
	 */
	int waveform = 0;

	/**
	 * Tube temperature (25 - 40 C)
	 */
	double temperature = 32.0;

	/**
	 * Junction loss factor in (0 - 5 %)
	 */
	double loss_factor = 0.8;

	/**
	 * Mouth aperture coefficient
	 */
	double mouth_coef = 5000.0;

	/**
	 * Nose aperture coefficient
	 */
	double nose_coef = 5000.0;

	/**
	 * Throat lp cutoff (>= 50)
	 */
	double throat_cutoff = 1500.0;

	/**
	 * Throat volume (0 - 48 dB)
	 */
	double throat_vol = 6.0;

	/**
	 * Pulse mod. of noise (0=OFF, 1=ON)
	 */
	int modulation = 1;

	/**
	 * Noise crossmix offset (30 - 60 dB)
	 */
	double mix_offset = 48.0;

	// Parameters that depend on the voice.

	/**
	 * % glottal pulse rise time
	 */
	double glottal_pulse_tp = 40.0;

	/**
	 * % glottal pulse fall time minimum
	 */
	double glottal_pulse_tn_min = 32.0;

	/**
	 * % glottal pulse fall time maximum
	 */
	double glottal_pulse_tn_max = 32.0;

	/**
	 * % glottal source breathiness
	 */
	double breathiness = 1.5;

	/**
	 * Length of the vocal tract.
	 */
	double vocal_tract_length = 15.0;

	/**
	 * Reference glottal pitch as a musical note (0 = C)
	 */
	double reference_glottal_pitch = 0.0;

	/**
	 * Aperture scl. radius (3.05 - 12 cm)
	 */
	double aperture_radius = 3.05;

	/**
	 * Fixed nose radii (0 - 3 cm)
	 */
	double nose_radius[TRM_CONF_NOSE_RADIUS_COUNT] = {1.35, 1.96, 1.91, 1.3,
	                                                  0.73};

	/**
	 * Radius coefficients.
	 */
	double radius_coef[TRM_CONF_RADIUS_COEF_COUNT] = {1.0, 1.0, 1.0, 1.0,
	                                                  1.0, 1.0, 1.0, 1.0};

	/**
	 * Set to true if the synthesizer needs to be re-initialized due to the
	 * update.
	 */
	bool updated = true;
};

/**
 * Structure holding all parameters as well as their differential.
 */
struct TrmParameters {
	double glot_pitch = 0.0;
	double glot_vol = 0.0;
	double asp_vol = 0.0;
	double fric_vol = 0.0;
	double fric_pos = 0.0;
	double fric_cf = 0.0;
	double fric_bw = 0.0;
	double radius[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double velum = 0.0;

	double dglot_pitch = 0.0;
	double dglot_vol = 0.0;
	double dasp_vol = 0.0;
	double dfric_vol = 0.0;
	double dfric_pos = 0.0;
	double dfric_cf = 0.0;
	double dfric_bw = 0.0;
	double dradius[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double dvelum = 0.0;
};

/**
 * Used to clip the state variables to a sane range.
 */
static inline double clip(double x, double min, double max) {
	return std::max(min, std::min(max, x));
}

/**
 * Actual object behing gnuspeech_trm_t. Holds the configuration, current
 * dynamics data and parameters, as well as
 */
struct TrmInstance {
	/**
	 * Structure holding the configuration options.
	 */
	TrmConfig config;

	/**
	 * Structure holding the current synthesizer parameters.
	 */
	TrmParameters params;

	/**
	 * Actual synthesizer instance.
	 */
	GS::TRM::Tube tube;

	/**
	 * Data position pointer.
	 */
	unsigned int data_pos = 0;

	/**
	 * Remainder used to accumulate very small requests to
	 * gnuspeech_trm_synthesize()
	 */
	double smpls_rem = 0.0;
};
}  // namespace

/******************************************************************************
 * C API                                                                      *
 ******************************************************************************/

gnuspeech_trm_t gnuspeech_trm_create() { return new TrmInstance(); }

void gnuspeech_trm_free(gnuspeech_trm_t inst_)
{
	if (inst_) {
		delete static_cast<TrmInstance *>(inst_);
	}
}

void gnuspeech_trm_reset(gnuspeech_trm_t inst_)
{
	if (inst_) {
		TrmInstance &inst = *(static_cast<TrmInstance *>(inst_));
		inst.config = TrmConfig();     /* Reset the configuration */
		inst.params = TrmParameters(); /* Reset the parameters */
		inst.tube.reset();             /* Reset the synthesizer */
		inst.data_pos = 0;
		inst.smpls_rem = 0.0;
	}
}

bool gnuspeech_trm_set_config(gnuspeech_trm_t inst_, int conf, double value)
{
#define SET_CONFIG_AND_RETURN(NAME) \
	if (c.NAME != value) {          \
		c.updated = true;           \
		c.NAME = value;             \
	}                               \
	return true;

	if (!inst_) {
		return false;
	}

	TrmConfig &c = static_cast<TrmInstance *>(inst_)->config;
	switch (conf) {
		case TRM_CONF_OUTPUT_RATE:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(output_rate);
			}
			return false;

		case TRM_CONF_CONTROL_RATE:
			if (value >= 1.0 && value <= 1000.0) {
				SET_CONFIG_AND_RETURN(control_rate);
			}
			return false;

		case TRM_CONF_FILTER_PERIOD:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(filter_period);
			}
			return false;

		case TRM_CONF_VOLUME:
			if (value >= 0.0 && value <= 60.0) {
				SET_CONFIG_AND_RETURN(volume);
			}
			return false;

		case TRM_CONF_WAVEFORM:
			if (value == 0.0 || value == 1.0) {
				SET_CONFIG_AND_RETURN(waveform)
			}
			return false;

		case TRM_CONF_TEMPERATURE:
			if (value >= 25.0 && value <= 40.0) {
				SET_CONFIG_AND_RETURN(temperature);
			}
			return false;

		case TRM_CONF_LOSS_FACTOR:
			if (value >= 0.0 && value <= 5.0) {
				SET_CONFIG_AND_RETURN(loss_factor);
			}
			return false;

		case TRM_CONF_MOUTH_COEF:
			if (value >= 0.0) {
				SET_CONFIG_AND_RETURN(mouth_coef);
			}
			return false;

		case TRM_CONF_NOSE_COEF:
			if (value >= 0.0) {
				SET_CONFIG_AND_RETURN(nose_coef);
			}
			return false;

		case TRM_CONF_THROAT_CUTOFF:
			if (value >= 50.0) {
				SET_CONFIG_AND_RETURN(throat_cutoff);
			}
			return false;

		case TRM_CONF_THROAT_VOL:
			if (value >= 0.0 && value <= 48.0) {
				SET_CONFIG_AND_RETURN(throat_vol);
			}
			return false;

		case TRM_CONF_MODULATION:
			if (value == 0.0 || value == 1.0) {
				SET_CONFIG_AND_RETURN(modulation);
			}
			return false;

		case TRM_CONF_MIX_OFFSET:
			if (value >= 30.0 && value <= 60.0) {
				SET_CONFIG_AND_RETURN(mix_offset);
			}
			return false;

		case TRM_CONF_GLOTTAL_PULSE_TP:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(glottal_pulse_tp);
			}
			return false;

		case TRM_CONF_GLOTTAL_PULSE_TN_MIN:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(glottal_pulse_tn_min);
			}
			return false;

		case TRM_CONF_GLOTTAL_PULSE_TN_MAX:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(glottal_pulse_tn_max);
			}
			return false;

		case TRM_CONF_BREATHINESS:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(breathiness);
			}
			return false;

		case TRM_CONF_VOCAL_TRACT_LENGTH:
			if (value > 0.0) {
				SET_CONFIG_AND_RETURN(vocal_tract_length);
			}
			return false;

		case TRM_CONF_REFERENCE_GLOTTAL_PITCH:
			SET_CONFIG_AND_RETURN(reference_glottal_pitch);

		case TRM_CONF_APERTURE_RADIUS:
			if (value >= 3.05 && value <= 12.0) {
				SET_CONFIG_AND_RETURN(aperture_radius);
			}
			return false;

		case TRM_CONF_NOSE_RADIUS_1:
		case TRM_CONF_NOSE_RADIUS_2:
		case TRM_CONF_NOSE_RADIUS_3:
		case TRM_CONF_NOSE_RADIUS_4:
		case TRM_CONF_NOSE_RADIUS_5:
		case TRM_CONF_NOSE_RADIUS_6: {
			if (value >= GS_TRM_TUBE_MIN_RADIUS && value <= 3.0) {
				const int i = conf - TRM_CONF_NOSE_RADIUS_MIN;
				if (c.nose_radius[i] != value) {
					c.nose_radius[i] = value;
					c.updated = true;
				}
				return true;
			}
			return false;
		}

		case TRM_CONF_RADIUS_COEF_1:
		case TRM_CONF_RADIUS_COEF_2:
		case TRM_CONF_RADIUS_COEF_3:
		case TRM_CONF_RADIUS_COEF_4:
		case TRM_CONF_RADIUS_COEF_5:
		case TRM_CONF_RADIUS_COEF_6:
		case TRM_CONF_RADIUS_COEF_7:
		case TRM_CONF_RADIUS_COEF_8: {
			if (value > GS_TRM_TUBE_MIN_RADIUS) {
				const int i = conf - TRM_CONF_RADIUS_COEF_1;
				if (c.radius_coef[i] != value) {
					c.radius_coef[i] = value;
					c.updated = true;
				}
				return true;
			}
			return false;
		}

		default:
			return false;
	}
#undef SET_CONFIG_AND_RETURN
}

bool gnuspeech_trm_get_config(const gnuspeech_trm_t inst_, int conf,
                              double *value)
{
#define RETURN_VALUE(NAME) \
	*value = c.NAME;       \
	return true;

	/* Abort if the given value is invalid */
	if (!inst_ || !value) {
		return false;
	}

	const TrmConfig &c = static_cast<TrmInstance *>(inst_)->config;
	switch (conf) {
		case TRM_CONF_OUTPUT_RATE:
			RETURN_VALUE(output_rate);

		case TRM_CONF_CONTROL_RATE:
			RETURN_VALUE(control_rate);

		case TRM_CONF_FILTER_PERIOD:
			RETURN_VALUE(filter_period);

		case TRM_CONF_VOLUME:
			RETURN_VALUE(volume);

		case TRM_CONF_WAVEFORM:
			RETURN_VALUE(waveform);

		case TRM_CONF_TEMPERATURE:
			RETURN_VALUE(temperature);

		case TRM_CONF_LOSS_FACTOR:
			RETURN_VALUE(loss_factor);

		case TRM_CONF_MOUTH_COEF:
			RETURN_VALUE(mouth_coef);

		case TRM_CONF_NOSE_COEF:
			RETURN_VALUE(nose_coef);

		case TRM_CONF_THROAT_CUTOFF:
			RETURN_VALUE(throat_cutoff);

		case TRM_CONF_THROAT_VOL:
			RETURN_VALUE(throat_vol);

		case TRM_CONF_MODULATION:
			RETURN_VALUE(modulation);

		case TRM_CONF_MIX_OFFSET:
			RETURN_VALUE(mix_offset);

		case TRM_CONF_GLOTTAL_PULSE_TP:
			RETURN_VALUE(glottal_pulse_tp);

		case TRM_CONF_GLOTTAL_PULSE_TN_MIN:
			RETURN_VALUE(glottal_pulse_tn_min);

		case TRM_CONF_GLOTTAL_PULSE_TN_MAX:
			RETURN_VALUE(glottal_pulse_tn_max);

		case TRM_CONF_BREATHINESS:
			RETURN_VALUE(breathiness);

		case TRM_CONF_VOCAL_TRACT_LENGTH:
			RETURN_VALUE(vocal_tract_length);

		case TRM_CONF_REFERENCE_GLOTTAL_PITCH:
			RETURN_VALUE(reference_glottal_pitch);

		case TRM_CONF_APERTURE_RADIUS:
			RETURN_VALUE(aperture_radius);

		case TRM_CONF_NOSE_RADIUS_1:
		case TRM_CONF_NOSE_RADIUS_2:
		case TRM_CONF_NOSE_RADIUS_3:
		case TRM_CONF_NOSE_RADIUS_4:
		case TRM_CONF_NOSE_RADIUS_5:
		case TRM_CONF_NOSE_RADIUS_6: {
			const int i = conf - TRM_CONF_NOSE_RADIUS_MIN;
			*value = c.nose_radius[i];
			return true;
		}

		case TRM_CONF_RADIUS_COEF_1:
		case TRM_CONF_RADIUS_COEF_2:
		case TRM_CONF_RADIUS_COEF_3:
		case TRM_CONF_RADIUS_COEF_4:
		case TRM_CONF_RADIUS_COEF_5:
		case TRM_CONF_RADIUS_COEF_6:
		case TRM_CONF_RADIUS_COEF_7:
		case TRM_CONF_RADIUS_COEF_8: {
			const int i = conf - TRM_CONF_RADIUS_COEF_MIN;
			*value = c.radius_coef[i];
			return true;
		}

		default:
			return false;
	}
#undef RETURN_VALUE
}

bool gnuspeech_trm_set_parameter(gnuspeech_trm_t inst_, int param, double value)
{
#define SET_PARAM_AND_RETURN(NAME) \
	p.NAME = value;                \
	return true;

	if (!inst_) {
		return false;
	}

	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	switch (param) {
		case TRM_PARAM_GLOT_PITCH:
			SET_PARAM_AND_RETURN(glot_pitch);
		case TRM_PARAM_GLOT_VOL:
			SET_PARAM_AND_RETURN(glot_vol);
		case TRM_PARAM_ASP_VOL:
			SET_PARAM_AND_RETURN(asp_vol);
		case TRM_PARAM_FRIC_VOL:
			SET_PARAM_AND_RETURN(fric_vol);
		case TRM_PARAM_FRIC_POS:
			SET_PARAM_AND_RETURN(fric_pos);
		case TRM_PARAM_FRIC_CF:
			SET_PARAM_AND_RETURN(fric_cf);
		case TRM_PARAM_FRIC_BW:
			SET_PARAM_AND_RETURN(fric_bw);
		case TRM_PARAM_R1:
		case TRM_PARAM_R2:
		case TRM_PARAM_R3:
		case TRM_PARAM_R4:
		case TRM_PARAM_R5:
		case TRM_PARAM_R6:
		case TRM_PARAM_R7:
		case TRM_PARAM_R8: {
			const int i = param - TRM_PARAM_R1;
			p.radius[i] = std::max(value, GS_TRM_TUBE_MIN_RADIUS);
			return true;
		}
		case TRM_PARAM_VELUM:
			SET_PARAM_AND_RETURN(velum);
		default:
			return false;
	}
#undef SET_PARAM_AND_RETURN
}

bool gnuspeech_trm_get_parameter(const gnuspeech_trm_t inst_, int param,
                                 double *value)
{
#define RETURN_VALUE(NAME) \
	*value = p.NAME;       \
	return true;

	if (!inst_ || !value || param < 0 || param >= TRM_PARAM_COUNT) {
		return false;
	}

	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	switch (param) {
		case TRM_PARAM_GLOT_PITCH:
			RETURN_VALUE(glot_pitch);
		case TRM_PARAM_GLOT_VOL:
			RETURN_VALUE(glot_vol);
		case TRM_PARAM_ASP_VOL:
			RETURN_VALUE(asp_vol);
		case TRM_PARAM_FRIC_VOL:
			RETURN_VALUE(fric_vol);
		case TRM_PARAM_FRIC_POS:
			RETURN_VALUE(fric_pos);
		case TRM_PARAM_FRIC_CF:
			RETURN_VALUE(fric_cf);
		case TRM_PARAM_FRIC_BW:
			RETURN_VALUE(fric_bw);
		case TRM_PARAM_R1:
		case TRM_PARAM_R2:
		case TRM_PARAM_R3:
		case TRM_PARAM_R4:
		case TRM_PARAM_R5:
		case TRM_PARAM_R6:
		case TRM_PARAM_R7:
		case TRM_PARAM_R8: {
			const int i = param - TRM_PARAM_R1;
			*value = p.radius[i];
			return true;
		}
		case TRM_PARAM_VELUM:
			RETURN_VALUE(velum);
		default:
			return false;
	}
}

bool gnuspeech_trm_get_parameters(gnuspeech_trm_t inst_, double *dx)
{
	if (!inst_) {
		return false;
	}
	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	dx[TRM_PARAM_GLOT_PITCH] = p.glot_pitch;
	dx[TRM_PARAM_GLOT_VOL] = p.glot_vol;
	dx[TRM_PARAM_ASP_VOL] = p.asp_vol;
	dx[TRM_PARAM_FRIC_VOL] = p.fric_vol;
	dx[TRM_PARAM_FRIC_POS] = p.fric_pos;
	dx[TRM_PARAM_FRIC_CF] = p.fric_cf;
	dx[TRM_PARAM_FRIC_BW] = p.fric_bw;
	for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
		dx[TRM_PARAM_R1 + j] = p.radius[j];
	}
	dx[TRM_PARAM_VELUM] = p.velum;
	return true;
}

bool gnuspeech_trm_set_parameters(gnuspeech_trm_t inst_, const double *dx)
{
	if (!inst_) {
		return false;
	}
	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	p.glot_pitch = dx[TRM_PARAM_GLOT_PITCH];
	p.glot_vol = dx[TRM_PARAM_GLOT_VOL];
	p.asp_vol = dx[TRM_PARAM_ASP_VOL];
	p.fric_vol = dx[TRM_PARAM_FRIC_VOL];
	p.fric_pos = dx[TRM_PARAM_FRIC_POS];
	p.fric_cf = dx[TRM_PARAM_FRIC_CF];
	p.fric_bw = dx[TRM_PARAM_FRIC_BW];
	for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
		p.radius[j] = dx[TRM_PARAM_R1 + j];
	}
	p.velum = dx[TRM_PARAM_VELUM];
	return true;
}

bool gnuspeech_trm_get_parameter_dynamics(gnuspeech_trm_t inst_, double *dx)
{
	if (!inst_) {
		return false;
	}
	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	dx[TRM_PARAM_GLOT_PITCH] = p.dglot_pitch;
	dx[TRM_PARAM_GLOT_VOL] = p.dglot_vol;
	dx[TRM_PARAM_ASP_VOL] = p.dasp_vol;
	dx[TRM_PARAM_FRIC_VOL] = p.dfric_vol;
	dx[TRM_PARAM_FRIC_POS] = p.dfric_pos;
	dx[TRM_PARAM_FRIC_CF] = p.dfric_cf;
	dx[TRM_PARAM_FRIC_BW] = p.dfric_bw;
	for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
		p.dradius[j] = dx[TRM_PARAM_R1 + j];
	}
	dx[TRM_PARAM_VELUM] = p.dvelum;
	return true;
}

bool gnuspeech_trm_set_parameter_dynamics(gnuspeech_trm_t inst_,
                                          const double *dx)
{
	if (!inst_) {
		return false;
	}
	TrmParameters &p = static_cast<TrmInstance *>(inst_)->params;
	p.dglot_pitch = dx[TRM_PARAM_GLOT_PITCH];
	p.dglot_vol = dx[TRM_PARAM_GLOT_VOL];
	p.dasp_vol = dx[TRM_PARAM_ASP_VOL];
	p.dfric_vol = dx[TRM_PARAM_FRIC_VOL];
	p.dfric_pos = dx[TRM_PARAM_FRIC_POS];
	p.dfric_cf = dx[TRM_PARAM_FRIC_CF];
	p.dfric_bw = dx[TRM_PARAM_FRIC_BW];
	for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
		p.dradius[j] = dx[TRM_PARAM_R1 + j];
	}
	p.dvelum = dx[TRM_PARAM_VELUM];
	return true;
}

int gnuspeech_trm_synthesize(gnuspeech_trm_t inst_, float *sample_buf,
                             unsigned int n_samples, bool flush)
{
	/* Check parameters for validity */
	if (!inst_ || !sample_buf) {
		return -1;
	}

	/* Fetch the actual instance */
	TrmInstance &inst = *(static_cast<TrmInstance *>(inst_));

	/* Some handy aliases */
	std::vector<float> &data = inst.tube.outputData();
	TrmConfig &c = inst.config;
	TrmParameters &p = inst.params;
	unsigned int samples_written = 0;

	/* Update the synthesizer if required */
	if (inst.config.updated) {
		/* Make sure to read remaining data */
		if (inst.tube.initialized()) {
			static const double scale = inst.tube.calculateMonoScale();
			while (inst.data_pos < data.size() &&
			       (samples_written < n_samples)) {
				sample_buf[samples_written++] = data[inst.data_pos++] * scale;
			}
			if (samples_written == n_samples) {
				return samples_written;
			}
		}

		/* Try to (re)initialize the synthesizer */
		try {
			inst.tube.loadConfig(inst.config);
			inst.tube.initializeSynthesizer();
			inst.tube.initializeInputFilters(c.filter_period);
			inst.data_pos = 0;
			inst.smpls_rem = 0.0;
			inst.config.updated = false;
		}
		catch (...) {
			return -1;
		}
	}

	/* Calculate the number of control samples */
	double scale = 0.0;
	const double dt = inst.tube.secPerSample();
	const double smpls_per_smpl = inst.tube.outputSamplesPerSample();
	int n_it = (n_samples + inst.smpls_rem) / smpls_per_smpl;
	inst.smpls_rem = (n_samples + inst.smpls_rem) - n_it * smpls_per_smpl;

	while ((samples_written < n_samples) &&
	       ((n_it > 0) || (inst.data_pos < data.size()))) {
		/* Check whether we need to synthesize new samples. If yes, update the
		   control parameters and run the synthesizer. */
		if (inst.data_pos >= data.size()) {
			/* Remove elements from the input queue */
			if (inst.data_pos > 0) {
				data.erase(data.begin(), data.begin() + inst.data_pos);
				inst.data_pos = 0;
			}

			/* Assemble the input data */
			GS::TRM::Tube::InputData input;
			input.glotPitch = p.glot_pitch + c.reference_glottal_pitch;
			input.glotVol = p.glot_vol;
			input.aspVol = p.asp_vol;
			input.fricVol = p.fric_vol;
			input.fricPos = p.fric_pos;
			input.fricCF = p.fric_cf;
			input.fricBW = p.fric_bw;
			for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
				input.radius[j] = std::max(p.radius[j] * c.radius_coef[j],
				                           GS_TRM_TUBE_MIN_RADIUS);
			}
			input.velum = p.velum;

			/* Run the synthesis for the given input value */
			inst.tube.synthesizeForSingleInput(input);

			/* Interpolate the control parameters */
			p.glot_pitch = clip(p.glot_pitch + p.dglot_pitch * dt, -20., 20.);
			p.glot_vol = clip(p.glot_vol + p.dglot_vol * dt, 0., 60.);
			p.asp_vol = clip(p.asp_vol + p.dasp_vol * dt, 0., 60.);
			p.fric_vol = clip(p.fric_vol + p.dfric_vol * dt, 0., 60.);
			p.fric_pos = clip(p.fric_pos + p.dfric_pos * dt, 0., 8.0);
			p.fric_cf = clip(p.fric_cf + p.dfric_cf * dt, 0., 4000.);
			p.fric_bw = clip(p.fric_bw + p.dfric_bw * dt, 0., 4000.);
			for (int j = 0; j < GS::TRM::Tube::TOTAL_REGIONS; j++) {
				p.radius[j] = clip(p.radius[j] + p.dradius[j] * dt,
				                           GS_TRM_TUBE_MIN_RADIUS, 1.0);
			}
			p.velum = clip(p.velum + p.dvelum * dt, 0.0, 1.0);
			n_it--;
		}
		else {
			/* Fetch the value used for volume normalisation */
			if (scale == 0.0) {
				scale = inst.tube.calculateMonoScale();
			}

			/* Write data from the output buffer to the target buffer */
			sample_buf[samples_written++] = data[inst.data_pos++] * scale;
		}
	}

	/* We've flushed the buffer, reset the sample remainder accumulator */
	if (flush) {
		inst.smpls_rem = 0.0;
	}

	return samples_written;
}

