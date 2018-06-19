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

/******************************************************************************
 *
 *     Program:       tube
 *
 *     Description:   Software (non-real-time) implementation of the Tube
 *                    Resonance Model for speech production.
 *
 *     Author:        Leonard Manzara
 *
 *     Date:          July 5th, 1994
 *
 ******************************************************************************/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility> /* move */

#include "Tube.h"

#define OUTPUT_VECTOR_RESERVE 1024

#define GLOTTAL_SOURCE_PULSE 0
#define GLOTTAL_SOURCE_SINE 1

/*  PITCH VARIABLES  */
#define PITCH_BASE 220.0
#define PITCH_OFFSET 3 /*  MIDDLE C = 0  */

/*  RANGE OF ALL VOLUME CONTROLS  */
#define VOL_MAX 60

/*  SCALING CONSTANT FOR INPUT TO VOCAL TRACT & THROAT (MATCHES DSP)  */
//#define VT_SCALE                  0.03125     /*  2^(-5)  */
// this is a temporary fix only, to try to match dsp synthesizer
#define VT_SCALE 0.125 /*  2^(-3)  */

/*  FINAL OUTPUT SCALING, SO THAT .SND FILES APPROX. MATCH DSP OUTPUT  */
#define OUTPUT_SCALE 0.95

/*  BI-DIRECTIONAL TRANSMISSION LINE POINTERS  */
#define TOP 0
#define BOTTOM 1

namespace GS {
namespace TRM {

Tube::Tube()
{
	reset();

	outputData_.reserve(OUTPUT_VECTOR_RESERVE);
}

Tube::~Tube() {}

void Tube::reset()
{
	outputRate_ = 0.0;
	controlRate_ = 0.0;
	volume_ = 0.0;
	waveform_ = 0;
	tp_ = 0.0;
	tnMin_ = 0.0;
	tnMax_ = 0.0;
	breathiness_ = 0.0;
	length_ = 0.0;
	temperature_ = 0.0;
	lossFactor_ = 0.0;
	apertureRadius_ = 0.0;
	mouthCoef_ = 0.0;
	noseCoef_ = 0.0;
	memset(noseRadius_, 0, sizeof(double) * TOTAL_NASAL_SECTIONS);
	throatCutoff_ = 0.0;
	throatVol_ = 0.0;
	modulation_ = 0;
	mixOffset_ = 0.0;
	controlPeriod_ = 0;
	sampleRate_ = 0;
	actualTubeLength_ = 0.0;
	memset(&oropharynx_[0][0][0], 0, sizeof(double) * TOTAL_SECTIONS * 2 * 2);
	memset(oropharynxCoeff_, 0, sizeof(double) * TOTAL_COEFFICIENTS);
	memset(&nasal_[0][0][0], 0, sizeof(double) * TOTAL_NASAL_SECTIONS * 2 * 2);
	memset(nasalCoeff_, 0, sizeof(double) * TOTAL_NASAL_COEFFICIENTS);
	memset(alpha_, 0, sizeof(double) * TOTAL_ALPHA_COEFFICIENTS);
	currentPtr_ = 1;
	prevPtr_ = 0;
	memset(fricationTap_, 0, sizeof(double) * TOTAL_FRIC_COEFFICIENTS);
	dampingFactor_ = 0.0;
	crossmixFactor_ = 0.0;
	breathinessFactor_ = 0.0;
	prevGlotAmplitude_ = -1.0;
	memset(&currentData_, 0, sizeof(InputData));
	outputData_.resize(0);

	if (srConv_)
		srConv_->reset();
	if (mouthRadiationFilter_)
		mouthRadiationFilter_->reset();
	if (mouthReflectionFilter_)
		mouthReflectionFilter_->reset();
	if (nasalRadiationFilter_)
		nasalRadiationFilter_->reset();
	if (nasalReflectionFilter_)
		nasalReflectionFilter_->reset();
	if (throat_)
		throat_->reset();
	if (glottalSource_)
		glottalSource_->reset();
	if (bandpassFilter_)
		bandpassFilter_->reset();
	if (noiseFilter_)
		noiseFilter_->reset();
	if (noiseSource_)
		noiseSource_->reset();
	if (inputFilters_)
		inputFilters_->reset();
}

/******************************************************************************
 *
 *  function:  speedOfSound
 *
 *  purpose:   Returns the speed of sound according to the value of
 *             the temperature (in Celsius degrees).
 *
 ******************************************************************************/
double Tube::speedOfSound(double temperature)
{
	return 331.4 + (0.6 * temperature);
}

/******************************************************************************
 *
 *  function:  initializeSynthesizer
 *
 *  purpose:   Initializes all variables so that the synthesis can
 *             be run.
 *
 ******************************************************************************/
void Tube::initializeSynthesizer()
{
	double nyquist;

	/*  CALCULATE THE SAMPLE RATE, BASED ON NOMINAL TUBE LENGTH AND SPEED OF
	 * SOUND  */
	if (length_ > 0.0) {
		double c = speedOfSound(temperature_);
		controlPeriod_ = static_cast<int>(
		    rint((c * TOTAL_SECTIONS * 100.0) / (length_ * controlRate_)));
		sampleRate_ = static_cast<int>(controlRate_ * controlPeriod_);
		actualTubeLength_ = (c * TOTAL_SECTIONS * 100.0) / sampleRate_;
		nyquist = sampleRate_ / 2.0;
	}
	else {
		throw std::runtime_error("Illegal tube length.");
	}

	/*  CALCULATE THE BREATHINESS FACTOR  */
	breathinessFactor_ = breathiness_ / 100.0;

	/*  CALCULATE CROSSMIX FACTOR  */
	crossmixFactor_ = 1.0 / amplitude(mixOffset_);

	/*  CALCULATE THE DAMPING FACTOR  */
	dampingFactor_ = (1.0 - (lossFactor_ / 100.0));

	/*  INITIALIZE THE WAVE TABLE  */
	glottalSource_.reset(new WavetableGlottalSource(
	    waveform_ == GLOTTAL_SOURCE_PULSE ? WavetableGlottalSource::TYPE_PULSE
	                                      : WavetableGlottalSource::TYPE_SINE,
	    sampleRate_, tp_, tnMin_, tnMax_));

	/*  INITIALIZE REFLECTION AND RADIATION FILTER COEFFICIENTS FOR MOUTH  */
	double mouthApertureCoeff = (nyquist - mouthCoef_) / nyquist;
	mouthRadiationFilter_.reset(new RadiationFilter(mouthApertureCoeff));
	mouthReflectionFilter_.reset(new ReflectionFilter(mouthApertureCoeff));

	/*  INITIALIZE REFLECTION AND RADIATION FILTER COEFFICIENTS FOR NOSE  */
	double nasalApertureCoeff = (nyquist - noseCoef_) / nyquist;
	nasalRadiationFilter_.reset(new RadiationFilter(nasalApertureCoeff));
	nasalReflectionFilter_.reset(new ReflectionFilter(nasalApertureCoeff));

	/*  INITIALIZE NASAL CAVITY FIXED SCATTERING COEFFICIENTS  */
	initializeNasalCavity();

	/*  INITIALIZE THE THROAT LOWPASS FILTER  */
	throat_.reset(
	    new Throat(sampleRate_, throatCutoff_, amplitude(throatVol_)));

	/*  INITIALIZE THE SAMPLE RATE CONVERSION ROUTINES  */
	srConv_.reset(
	    new SampleRateConverter(sampleRate_, outputRate_, outputData_));

	/*  INITIALIZE THE OUTPUT VECTOR  */
	outputData_.clear();

	bandpassFilter_.reset(new BandpassFilter());
	noiseFilter_.reset(new NoiseFilter());
	noiseSource_.reset(new NoiseSource());
}

void Tube::initializeInputFilters(double period)
{
	inputFilters_.reset(new InputFilters(sampleRate_, period));
}

void Tube::synthesizeForSingleInput(const InputData &inputData)
{
	currentData_.glotPitch =
	    inputFilters_->glotPitchFilter.filter(inputData.glotPitch);
	currentData_.glotVol =
	    inputFilters_->glotVolFilter.filter(inputData.glotVol);
	currentData_.aspVol = inputFilters_->aspVolFilter.filter(inputData.aspVol);
	currentData_.fricVol =
	    inputFilters_->fricVolFilter.filter(inputData.fricVol);
	currentData_.fricPos =
	    inputFilters_->fricPosFilter.filter(inputData.fricPos);
	currentData_.fricCF = inputFilters_->fricCFFilter.filter(inputData.fricCF);
	currentData_.fricBW = inputFilters_->fricBWFilter.filter(inputData.fricBW);
	currentData_.radius[0] =
	    inputFilters_->radius0Filter.filter(inputData.radius[0]);
	currentData_.radius[1] =
	    inputFilters_->radius1Filter.filter(inputData.radius[1]);
	currentData_.radius[2] =
	    inputFilters_->radius2Filter.filter(inputData.radius[2]);
	currentData_.radius[3] =
	    inputFilters_->radius3Filter.filter(inputData.radius[3]);
	currentData_.radius[4] =
	    inputFilters_->radius4Filter.filter(inputData.radius[4]);
	currentData_.radius[5] =
	    inputFilters_->radius5Filter.filter(inputData.radius[5]);
	currentData_.radius[6] =
	    inputFilters_->radius6Filter.filter(inputData.radius[6]);
	currentData_.radius[7] =
	    inputFilters_->radius7Filter.filter(inputData.radius[7]);
	currentData_.velum = inputFilters_->velumFilter.filter(inputData.velum);

	synthesize();
}

void Tube::flush() { srConv_->flushBuffer(); }

/******************************************************************************
 *
 *  function:  synthesize
 *
 *  purpose:   Performs the actual synthesis of sound samples.
 *
 ******************************************************************************/

void Tube::synthesize()
{
	/*  CONVERT PARAMETERS HERE  */
	double f0 = frequency(currentData_.glotPitch);
	double ax = amplitude(currentData_.glotVol);
	double ah1 = amplitude(currentData_.aspVol);
	calculateTubeCoefficients();
	setFricationTaps();
	bandpassFilter_->update(sampleRate_, currentData_.fricBW,
	                        currentData_.fricCF);

	/*  DO SYNTHESIS HERE  */
	/*  CREATE LOW-PASS FILTERED NOISE  */
	double lpNoise = noiseFilter_->filter(noiseSource_->getSample());

	/*  UPDATE THE SHAPE OF THE GLOTTAL PULSE, IF NECESSARY  */
	if (waveform_ == GLOTTAL_SOURCE_PULSE) {
		if (ax != prevGlotAmplitude_) {
			glottalSource_->updateWavetable(ax);
		}
	}

	/*  CREATE GLOTTAL PULSE (OR SINE TONE)  */
	double pulse = glottalSource_->getSample(f0);

	/*  CREATE PULSED NOISE  */
	double pulsedNoise = lpNoise * pulse;

	/*  CREATE NOISY GLOTTAL PULSE  */
	pulse = ax * ((pulse * (1.0 - breathinessFactor_)) +
	              (pulsedNoise * breathinessFactor_));

	double signal;
	/*  CROSS-MIX PURE NOISE WITH PULSED NOISE  */
	if (modulation_) {
		double crossmix = ax * crossmixFactor_;
		crossmix = (crossmix < 1.0) ? crossmix : 1.0;
		signal = (pulsedNoise * crossmix) + (lpNoise * (1.0 - crossmix));
	}
	else {
		signal = lpNoise;
	}

	/*  PUT SIGNAL THROUGH VOCAL TRACT  */
	signal = vocalTract(((pulse + (ah1 * signal)) * VT_SCALE),
	                    bandpassFilter_->filter(signal));

	/*  PUT PULSE THROUGH THROAT  */
	signal += throat_->process(pulse * VT_SCALE);

	/*  OUTPUT SAMPLE HERE  */
	srConv_->dataFill(signal);

	prevGlotAmplitude_ = ax;
}

/******************************************************************************
 *
 *  function:  initializeNasalCavity
 *
 *  purpose:   Calculates the scattering coefficients for the fixed
 *             sections of the nasal cavity.
 *
 ******************************************************************************/
void Tube::initializeNasalCavity()
{
	double radA2, radB2;

	/*  CALCULATE COEFFICIENTS FOR INTERNAL FIXED SECTIONS OF NASAL CAVITY  */
	for (int i = N2, j = NC2; i < N6; i++, j++) {
		radA2 = noseRadius_[i] * noseRadius_[i];
		radB2 = noseRadius_[i + 1] * noseRadius_[i + 1];
		nasalCoeff_[j] = (radA2 - radB2) / (radA2 + radB2);
	}

	/*  CALCULATE THE FIXED COEFFICIENT FOR THE NOSE APERTURE  */
	radA2 = noseRadius_[N6] * noseRadius_[N6];
	radB2 = apertureRadius_ * apertureRadius_;
	nasalCoeff_[NC6] = (radA2 - radB2) / (radA2 + radB2);
}

/******************************************************************************
 *
 *  function:  calculateTubeCoefficients
 *
 *  purpose:   Calculates the scattering coefficients for the vocal
 *             ract according to the current radii.  Also calculates
 *             the coefficients for the reflection/radiation filter
 *             pair for the mouth and nose.
 *
 ******************************************************************************/
void Tube::calculateTubeCoefficients()
{
	double radA2, radB2, r0_2, r1_2, r2_2, sum;

	/*  CALCULATE COEFFICIENTS FOR THE OROPHARYNX  */
	for (int i = 0; i < (TOTAL_REGIONS - 1); i++) {
		radA2 = currentData_.radius[i] * currentData_.radius[i];
		radB2 = currentData_.radius[i + 1] * currentData_.radius[i + 1];
		oropharynxCoeff_[i] = (radA2 - radB2) / (radA2 + radB2);
	}

	/*  CALCULATE THE COEFFICIENT FOR THE MOUTH APERTURE  */
	radA2 = currentData_.radius[R8] * currentData_.radius[R8];
	radB2 = apertureRadius_ * apertureRadius_;
	oropharynxCoeff_[C8] = (radA2 - radB2) / (radA2 + radB2);

	/*  CALCULATE ALPHA COEFFICIENTS FOR 3-WAY JUNCTION  */
	/*  NOTE:  SINCE JUNCTION IS IN MIDDLE OF REGION 4, r0_2 = r1_2  */
	r0_2 = r1_2 = currentData_.radius[R4] * currentData_.radius[R4];
	r2_2 = currentData_.velum * currentData_.velum;
	sum = 2.0 / (r0_2 + r1_2 + r2_2);
	alpha_[LEFT] = sum * r0_2;
	alpha_[RIGHT] = sum * r1_2;
	alpha_[UPPER] = sum * r2_2;

	/*  AND 1ST NASAL PASSAGE COEFFICIENT  */
	radA2 = currentData_.velum * currentData_.velum;
	radB2 = noseRadius_[N2] * noseRadius_[N2];
	nasalCoeff_[NC1] = (radA2 - radB2) / (radA2 + radB2);
}

/******************************************************************************
 *
 *  function:  setFricationTaps
 *
 *  purpose:   Sets the frication taps according to the current
 *             position and amplitude of frication.
 *
 ******************************************************************************/
void Tube::setFricationTaps()
{
	int integerPart;
	double complement, remainder;
	double fricationAmplitude = amplitude(currentData_.fricVol);

	/*  CALCULATE POSITION REMAINDER AND COMPLEMENT  */
	integerPart = (int)currentData_.fricPos;
	complement = currentData_.fricPos - (double)integerPart;
	remainder = 1.0 - complement;

	/*  SET THE FRICATION TAPS  */
	for (int i = FC1; i < TOTAL_FRIC_COEFFICIENTS; i++) {
		if (i == integerPart) {
			fricationTap_[i] = remainder * fricationAmplitude;
			if ((i + 1) < TOTAL_FRIC_COEFFICIENTS) {
				fricationTap_[++i] = complement * fricationAmplitude;
			}
		}
		else {
			fricationTap_[i] = 0.0;
		}
	}
}

/******************************************************************************
 *
 *  function:  vocalTract
 *
 *  purpose:   Updates the pressure wave throughout the vocal tract,
 *             and returns the summed output of the oral and nasal
 *             cavities.  Also injects frication appropriately.
 *
 ******************************************************************************/
double Tube::vocalTract(double input, double frication)
{
	int i, j, k;
	double delta, output, junctionPressure;

	/*  INCREMENT CURRENT AND PREVIOUS POINTERS  */
	if (++currentPtr_ > 1) {
		currentPtr_ = 0;
	}
	if (++prevPtr_ > 1) {
		prevPtr_ = 0;
	}

	/*  UPDATE OROPHARYNX  */
	/*  INPUT TO TOP OF TUBE  */
	oropharynx_[S1][TOP][currentPtr_] =
	    (oropharynx_[S1][BOTTOM][prevPtr_] * dampingFactor_) + input;

	/*  CALCULATE THE SCATTERING JUNCTIONS FOR S1-S2  */
	delta = oropharynxCoeff_[C1] * (oropharynx_[S1][TOP][prevPtr_] -
	                                oropharynx_[S2][BOTTOM][prevPtr_]);
	oropharynx_[S2][TOP][currentPtr_] =
	    (oropharynx_[S1][TOP][prevPtr_] + delta) * dampingFactor_;
	oropharynx_[S1][BOTTOM][currentPtr_] =
	    (oropharynx_[S2][BOTTOM][prevPtr_] + delta) * dampingFactor_;

	/*  CALCULATE THE SCATTERING JUNCTIONS FOR S2-S3 AND S3-S4  */
	for (i = S2, j = C2, k = FC1; i < S4; i++, j++, k++) {
		delta = oropharynxCoeff_[j] * (oropharynx_[i][TOP][prevPtr_] -
		                               oropharynx_[i + 1][BOTTOM][prevPtr_]);
		oropharynx_[i + 1][TOP][currentPtr_] =
		    ((oropharynx_[i][TOP][prevPtr_] + delta) * dampingFactor_) +
		    (fricationTap_[k] * frication);
		oropharynx_[i][BOTTOM][currentPtr_] =
		    (oropharynx_[i + 1][BOTTOM][prevPtr_] + delta) * dampingFactor_;
	}

	/*  UPDATE 3-WAY JUNCTION BETWEEN THE MIDDLE OF R4 AND NASAL CAVITY  */
	junctionPressure = (alpha_[LEFT] * oropharynx_[S4][TOP][prevPtr_]) +
	                   (alpha_[RIGHT] * oropharynx_[S5][BOTTOM][prevPtr_]) +
	                   (alpha_[UPPER] * nasal_[VELUM][BOTTOM][prevPtr_]);
	oropharynx_[S4][BOTTOM][currentPtr_] =
	    (junctionPressure - oropharynx_[S4][TOP][prevPtr_]) * dampingFactor_;
	oropharynx_[S5][TOP][currentPtr_] =
	    ((junctionPressure - oropharynx_[S5][BOTTOM][prevPtr_]) *
	     dampingFactor_) +
	    (fricationTap_[FC3] * frication);
	nasal_[VELUM][TOP][currentPtr_] =
	    (junctionPressure - nasal_[VELUM][BOTTOM][prevPtr_]) * dampingFactor_;

	/*  CALCULATE JUNCTION BETWEEN R4 AND R5 (S5-S6)  */
	delta = oropharynxCoeff_[C4] * (oropharynx_[S5][TOP][prevPtr_] -
	                                oropharynx_[S6][BOTTOM][prevPtr_]);
	oropharynx_[S6][TOP][currentPtr_] =
	    ((oropharynx_[S5][TOP][prevPtr_] + delta) * dampingFactor_) +
	    (fricationTap_[FC4] * frication);
	oropharynx_[S5][BOTTOM][currentPtr_] =
	    (oropharynx_[S6][BOTTOM][prevPtr_] + delta) * dampingFactor_;

	/*  CALCULATE JUNCTION INSIDE R5 (S6-S7) (PURE DELAY WITH DAMPING)  */
	oropharynx_[S7][TOP][currentPtr_] =
	    (oropharynx_[S6][TOP][prevPtr_] * dampingFactor_) +
	    (fricationTap_[FC5] * frication);
	oropharynx_[S6][BOTTOM][currentPtr_] =
	    oropharynx_[S7][BOTTOM][prevPtr_] * dampingFactor_;

	/*  CALCULATE LAST 3 INTERNAL JUNCTIONS (S7-S8, S8-S9, S9-S10)  */
	for (i = S7, j = C5, k = FC6; i < S10; i++, j++, k++) {
		delta = oropharynxCoeff_[j] * (oropharynx_[i][TOP][prevPtr_] -
		                               oropharynx_[i + 1][BOTTOM][prevPtr_]);
		oropharynx_[i + 1][TOP][currentPtr_] =
		    ((oropharynx_[i][TOP][prevPtr_] + delta) * dampingFactor_) +
		    (fricationTap_[k] * frication);
		oropharynx_[i][BOTTOM][currentPtr_] =
		    (oropharynx_[i + 1][BOTTOM][prevPtr_] + delta) * dampingFactor_;
	}

	/*  REFLECTED SIGNAL AT MOUTH GOES THROUGH A LOWPASS FILTER  */
	oropharynx_[S10][BOTTOM][currentPtr_] =
	    dampingFactor_ *
	    mouthReflectionFilter_->filter(oropharynxCoeff_[C8] *
	                                   oropharynx_[S10][TOP][prevPtr_]);

	/*  OUTPUT FROM MOUTH GOES THROUGH A HIGHPASS FILTER  */
	output = mouthRadiationFilter_->filter((1.0 + oropharynxCoeff_[C8]) *
	                                       oropharynx_[S10][TOP][prevPtr_]);

	/*  UPDATE NASAL CAVITY  */
	for (i = VELUM, j = NC1; i < N6; i++, j++) {
		delta = nasalCoeff_[j] *
		        (nasal_[i][TOP][prevPtr_] - nasal_[i + 1][BOTTOM][prevPtr_]);
		nasal_[i + 1][TOP][currentPtr_] =
		    (nasal_[i][TOP][prevPtr_] + delta) * dampingFactor_;
		nasal_[i][BOTTOM][currentPtr_] =
		    (nasal_[i + 1][BOTTOM][prevPtr_] + delta) * dampingFactor_;
	}

	/*  REFLECTED SIGNAL AT NOSE GOES THROUGH A LOWPASS FILTER  */
	nasal_[N6][BOTTOM][currentPtr_] =
	    dampingFactor_ * nasalReflectionFilter_->filter(
	                         nasalCoeff_[NC6] * nasal_[N6][TOP][prevPtr_]);

	/*  OUTPUT FROM NOSE GOES THROUGH A HIGHPASS FILTER  */
	output += nasalRadiationFilter_->filter((1.0 + nasalCoeff_[NC6]) *
	                                        nasal_[N6][TOP][prevPtr_]);
	/*  RETURN SUMMED OUTPUT FROM MOUTH AND NOSE  */
	return output;
}

float Tube::calculateMonoScale()
{
	float scale = static_cast<float>(
	    (OUTPUT_SCALE / std::max(1e-4, srConv_->maximumSampleValue())) *
	    amplitude(volume_));
	return scale;
}

/******************************************************************************
 *
 *  function:  amplitude
 *
 *  purpose:   Converts dB value to amplitude value.
 *
 ******************************************************************************/
double Tube::amplitude(double decibelLevel)
{
	/*  CONVERT 0-60 RANGE TO -60-0 RANGE  */
	decibelLevel -= VOL_MAX;

	/*  IF -60 OR LESS, RETURN AMPLITUDE OF 0  */
	if (decibelLevel <= (-VOL_MAX)) {
		return 0.0;
	}

	/*  IF 0 OR GREATER, RETURN AMPLITUDE OF 1  */
	if (decibelLevel >= 0.0) {
		return 1.0;
	}

	/*  ELSE RETURN INVERSE LOG VALUE  */
	return pow(10.0, decibelLevel / 20.0);
}

/******************************************************************************
 *
 *  function:  frequency
 *
 *  purpose:   Converts a given pitch (0 = middle C) to the
 *             corresponding frequency.
 *
 ******************************************************************************/
double Tube::frequency(double pitch)
{
	return PITCH_BASE * pow(2.0, (pitch + PITCH_OFFSET) / 12.0);
}
} /* namespace TRM */
} /* namespace GS */
