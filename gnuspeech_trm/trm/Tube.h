/***************************************************************************
 *  Copyright 1991, 1992, 1993, 1994, 1995, 1996, 2001, 2002               *
 *    David R. Hill, Leonard Manzara, Craig Schock                         *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/
// 2014-09
// This file was copied from Gnuspeech and modified by Marcelo Y. Matuda.

#ifndef TRM_TUBE_H_
#define TRM_TUBE_H_

#include <algorithm> /* max, min */
#include <memory>
#include <vector>

#include "BandpassFilter.h"
#include "MovingAverageFilter.h"
#include "NoiseFilter.h"
#include "NoiseSource.h"
#include "RadiationFilter.h"
#include "ReflectionFilter.h"
#include "SampleRateConverter.h"
#include "Throat.h"
#include "WavetableGlottalSource.h"

#define GS_TRM_TUBE_MIN_RADIUS (0.001)

namespace GS {
namespace TRM {

struct VocalTractModelParameterValue {
	int index;
	double value;
};

class Tube {
public:
	enum {         /*  OROPHARYNX REGIONS  */
		   R1 = 0, /*  S1  */
		   R2 = 1, /*  S2  */
		   R3 = 2, /*  S3  */
		   R4 = 3, /*  S4 & S5  */
		   R5 = 4, /*  S6 & S7  */
		   R6 = 5, /*  S8  */
		   R7 = 6, /*  S9  */
		   R8 = 7, /*  S10  */
		   TOTAL_REGIONS = 8
	};
	enum { /*  NASAL TRACT SECTIONS  */
		   N1 = 0,
		   N2 = 1,
		   N3 = 2,
		   N4 = 3,
		   N5 = 4,
		   N6 = 5,
		   TOTAL_NASAL_SECTIONS = 6
	};
	enum ParameterIndex {
		PARAM_GLOT_PITCH = 0,
		PARAM_GLOT_VOL = 1,
		PARAM_ASP_VOL = 2,
		PARAM_FRIC_VOL = 3,
		PARAM_FRIC_POS = 4,
		PARAM_FRIC_CF = 5,
		PARAM_FRIC_BW = 6,
		PARAM_R1 = 7,
		PARAM_R2 = 8,
		PARAM_R3 = 9,
		PARAM_R4 = 10,
		PARAM_R5 = 11,
		PARAM_R6 = 12,
		PARAM_R7 = 13,
		PARAM_R8 = 14,
		PARAM_VELUM = 15
	};

	struct InputData {
		double glotPitch;
		double glotVol;
		double aspVol;
		double fricVol;
		double fricPos;
		double fricCF;
		double fricBW;
		double radius[TOTAL_REGIONS];
		double velum;
	};

	Tube();
	~Tube();

	template <typename T>
	void loadConfig(const T &config);
	bool initialized() const { return bool(srConv_); }
	void initializeSynthesizer();
	void initializeInputFilters(double period);

	void synthesizeForSingleInput(const InputData &inputData);
	float calculateMonoScale();

	void reset();
	void flush();

	int outputSamplesPerSample() const {
		return double(outputRate_) / double(sampleRate_);
	}

	double secPerSample() const {
		return 1.0 / double(sampleRate_);
	}

	std::vector<float>& outputData() { return outputData_; }

private:
	enum { VELUM = N1 };
	enum { /*  OROPHARYNX SCATTERING JUNCTION COEFFICIENTS (BETWEEN EACH REGION)
		    */
		   C1 = R1, /*  R1-R2 (S1-S2)  */
		   C2 = R2, /*  R2-R3 (S2-S3)  */
		   C3 = R3, /*  R3-R4 (S3-S4)  */
		   C4 = R4, /*  R4-R5 (S5-S6)  */
		   C5 = R5, /*  R5-R6 (S7-S8)  */
		   C6 = R6, /*  R6-R7 (S8-S9)  */
		   C7 = R7, /*  R7-R8 (S9-S10)  */
		   C8 = R8, /*  R8-AIR (S10-AIR)  */
		   TOTAL_COEFFICIENTS = TOTAL_REGIONS
	};
	enum {          /*  OROPHARYNX SECTIONS  */
		   S1 = 0,  /*  R1  */
		   S2 = 1,  /*  R2  */
		   S3 = 2,  /*  R3  */
		   S4 = 3,  /*  R4  */
		   S5 = 4,  /*  R4  */
		   S6 = 5,  /*  R5  */
		   S7 = 6,  /*  R5  */
		   S8 = 7,  /*  R6  */
		   S9 = 8,  /*  R7  */
		   S10 = 9, /*  R8  */
		   TOTAL_SECTIONS = 10
	};
	enum {           /*  NASAL TRACT COEFFICIENTS  */
		   NC1 = N1, /*  N1-N2  */
		   NC2 = N2, /*  N2-N3  */
		   NC3 = N3, /*  N3-N4  */
		   NC4 = N4, /*  N4-N5  */
		   NC5 = N5, /*  N5-N6  */
		   NC6 = N6, /*  N6-AIR  */
		   TOTAL_NASAL_COEFFICIENTS = TOTAL_NASAL_SECTIONS
	};
	enum { /*  THREE-WAY JUNCTION ALPHA COEFFICIENTS  */
		   LEFT = 0,
		   RIGHT = 1,
		   UPPER = 2,
		   TOTAL_ALPHA_COEFFICIENTS = 3
	};
	enum {          /*  FRICATION INJECTION COEFFICIENTS  */
		   FC1 = 0, /*  S3  */
		   FC2 = 1, /*  S4  */
		   FC3 = 2, /*  S5  */
		   FC4 = 3, /*  S6  */
		   FC5 = 4, /*  S7  */
		   FC6 = 5, /*  S8  */
		   FC7 = 6, /*  S9  */
		   FC8 = 7, /*  S10  */
		   TOTAL_FRIC_COEFFICIENTS = 8
	};

	struct InputFilters {
		MovingAverageFilter<double> glotPitchFilter;
		MovingAverageFilter<double> glotVolFilter;
		MovingAverageFilter<double> aspVolFilter;
		MovingAverageFilter<double> fricVolFilter;
		MovingAverageFilter<double> fricPosFilter;
		MovingAverageFilter<double> fricCFFilter;
		MovingAverageFilter<double> fricBWFilter;
		MovingAverageFilter<double> radius0Filter;
		MovingAverageFilter<double> radius1Filter;
		MovingAverageFilter<double> radius2Filter;
		MovingAverageFilter<double> radius3Filter;
		MovingAverageFilter<double> radius4Filter;
		MovingAverageFilter<double> radius5Filter;
		MovingAverageFilter<double> radius6Filter;
		MovingAverageFilter<double> radius7Filter;
		MovingAverageFilter<double> velumFilter;
		InputFilters(double sampleRate, double period)
		    : glotPitchFilter(sampleRate, period),
		      glotVolFilter(sampleRate, period),
		      aspVolFilter(sampleRate, period),
		      fricVolFilter(sampleRate, period),
		      fricPosFilter(sampleRate, period),
		      fricCFFilter(sampleRate, period),
		      fricBWFilter(sampleRate, period),
		      radius0Filter(sampleRate, period),
		      radius1Filter(sampleRate, period),
		      radius2Filter(sampleRate, period),
		      radius3Filter(sampleRate, period),
		      radius4Filter(sampleRate, period),
		      radius5Filter(sampleRate, period),
		      radius6Filter(sampleRate, period),
		      radius7Filter(sampleRate, period),
		      velumFilter(sampleRate, period)
		{
		}
		void reset()
		{
			glotPitchFilter.reset();
			glotVolFilter.reset();
			aspVolFilter.reset();
			fricVolFilter.reset();
			fricPosFilter.reset();
			fricCFFilter.reset();
			fricBWFilter.reset();
			radius0Filter.reset();
			radius1Filter.reset();
			radius2Filter.reset();
			radius3Filter.reset();
			radius4Filter.reset();
			radius5Filter.reset();
			radius6Filter.reset();
			radius7Filter.reset();
			velumFilter.reset();
		}
	};

	Tube(const Tube &) = delete;
	Tube &operator=(const Tube &) = delete;

	void synthesize();
	void calculateTubeCoefficients();
	void initializeNasalCavity();
	void setControlRateParameters(int pos);
	void setFricationTaps();
	double vocalTract(double input, double frication);

	static double amplitude(double decibelLevel);
	static double frequency(double pitch);
	static double speedOfSound(double temperature);

	float outputRate_;  /*  output sample rate (22.05, 44.1)  */
	float controlRate_; /*  1.0-1000.0 input tables/second (Hz)  */

	double volume_; /*  master volume (0 - 60 dB)  */

	int waveform_;       /*  GS waveform type (0=PULSE, 1=SINE  */
	double tp_;          /*  % glottal pulse rise time  */
	double tnMin_;       /*  % glottal pulse fall time minimum  */
	double tnMax_;       /*  % glottal pulse fall time maximum  */
	double breathiness_; /*  % glottal source breathiness  */

	double length_;      /*  nominal tube length (10 - 20 cm)  */
	double temperature_; /*  tube temperature (25 - 40 C)  */
	double lossFactor_;  /*  junction loss factor in (0 - 5 %)  */

	double apertureRadius_; /*  aperture scl. radius (3.05 - 12 cm)  */
	double mouthCoef_;      /*  mouth aperture coefficient  */
	double noseCoef_;       /*  nose aperture coefficient  */

	double
	    noseRadius_[TOTAL_NASAL_SECTIONS]; /*  fixed nose radii (0 - 3 cm)  */

	double throatCutoff_; /*  throat lp cutoff (50 - nyquist Hz)  */
	double throatVol_;    /*  throat volume (0 - 48 dB) */

	int modulation_;   /*  pulse mod. of noise (0=OFF, 1=ON)  */
	double mixOffset_; /*  noise crossmix offset (30 - 60 dB)  */

	/*  DERIVED VALUES  */
	int controlPeriod_;
	int sampleRate_;
	double actualTubeLength_; /*  actual length in cm  */

	/*  MEMORY FOR TUBE AND TUBE COEFFICIENTS  */
	double oropharynx_[TOTAL_SECTIONS][2][2];
	double oropharynxCoeff_[TOTAL_COEFFICIENTS];

	double nasal_[TOTAL_NASAL_SECTIONS][2][2];
	double nasalCoeff_[TOTAL_NASAL_COEFFICIENTS];

	double alpha_[TOTAL_ALPHA_COEFFICIENTS];
	int currentPtr_;
	int prevPtr_;

	/*  MEMORY FOR FRICATION TAPS  */
	double fricationTap_[TOTAL_FRIC_COEFFICIENTS];

	double dampingFactor_;  /*  calculated damping factor  */
	double crossmixFactor_; /*  calculated crossmix factor  */
	double breathinessFactor_;

	double prevGlotAmplitude_;

	InputData currentData_;
	std::vector<float> outputData_;
	std::unique_ptr<SampleRateConverter> srConv_;
	std::unique_ptr<RadiationFilter> mouthRadiationFilter_;
	std::unique_ptr<ReflectionFilter> mouthReflectionFilter_;
	std::unique_ptr<RadiationFilter> nasalRadiationFilter_;
	std::unique_ptr<ReflectionFilter> nasalReflectionFilter_;
	std::unique_ptr<Throat> throat_;
	std::unique_ptr<WavetableGlottalSource> glottalSource_;
	std::unique_ptr<BandpassFilter> bandpassFilter_;
	std::unique_ptr<NoiseFilter> noiseFilter_;
	std::unique_ptr<NoiseSource> noiseSource_;
	std::unique_ptr<InputFilters> inputFilters_;
};

template <typename T>
void Tube::loadConfig(const T &config)
{
	outputRate_ = config.output_rate;
	controlRate_ = config.control_rate;
	volume_ = config.volume;
	waveform_ = config.waveform;
	temperature_ = config.temperature;
	lossFactor_ = config.loss_factor;
	mouthCoef_ = config.mouth_coef;
	noseCoef_ = config.nose_coef;
	throatCutoff_ = config.throat_cutoff;
	throatVol_ = config.throat_vol;
	modulation_ = config.modulation;
	mixOffset_ = config.mix_offset;
	tp_ = config.glottal_pulse_tp;
	tnMin_ = config.glottal_pulse_tn_min;
	tnMax_ = config.glottal_pulse_tn_max;
	length_ = config.vocal_tract_length;
	breathiness_ = config.breathiness;
	apertureRadius_ = config.aperture_radius;

	for (int i = 0; i < Tube::TOTAL_NASAL_COEFFICIENTS; i++) {
		noseRadius_[i] = config.nose_radius[i];
	}
}
} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_TUBE_H_ */
