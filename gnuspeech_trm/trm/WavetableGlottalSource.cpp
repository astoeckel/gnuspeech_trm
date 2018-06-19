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

#include <cmath>

#include "FIRFilter.h"
#include "WavetableGlottalSource.h"

/*  COMPILE WITH OVERSAMPLING OR PLAIN OSCILLATOR  */
#define OVERSAMPLING_OSCILLATOR 1

/*  GLOTTAL SOURCE OSCILLATOR TABLE VARIABLES  */
#define TABLE_LENGTH 512
#define TABLE_MODULUS (TABLE_LENGTH - 1)

/*  OVERSAMPLING FIR FILTER CHARACTERISTICS  */
#define FIR_BETA .2
#define FIR_GAMMA .1
#define FIR_CUTOFF .00000001

namespace GS {
namespace TRM {

WavetableGlottalSource::WavetableGlottalSource(Type type, double sampleRate,
                                               double tp, double tnMin,
                                               double tnMax)
    : wavetable_(TABLE_LENGTH)
{
	// Calculates the initial glottal pulse and stores it
	// in the wavetable, for use in the oscillator.

	/*  CALCULATE WAVE TABLE PARAMETERS  */
	tableDiv1_ = static_cast<int>(rint(TABLE_LENGTH * (tp / 100.0)));
	tableDiv2_ = static_cast<int>(rint(TABLE_LENGTH * ((tp + tnMax) / 100.0)));
	tnLength_ = tableDiv2_ - tableDiv1_;
	tnDelta_ = rint(TABLE_LENGTH * ((tnMax - tnMin) / 100.0));
	basicIncrement_ = (double)TABLE_LENGTH / (double)sampleRate;
	currentPosition_ = 0;

	/*  INITIALIZE THE WAVETABLE WITH EITHER A GLOTTAL PULSE OR SINE TONE  */
	if (type == TYPE_PULSE) {
		/*  CALCULATE RISE PORTION OF WAVE TABLE  */
		for (int i = 0; i < tableDiv1_; i++) {
			double x = (double)i / (double)tableDiv1_;
			double x2 = x * x;
			double x3 = x2 * x;
			wavetable_[i] = (3.0 * x2) - (2.0 * x3);
		}

		/*  CALCULATE FALL PORTION OF WAVE TABLE  */
		for (int i = tableDiv1_, j = 0; i < tableDiv2_; i++, j++) {
			double x = (double)j / tnLength_;
			wavetable_[i] = 1.0 - (x * x);
		}

		/*  SET CLOSED PORTION OF WAVE TABLE  */
		for (int i = tableDiv2_; i < TABLE_LENGTH; i++) {
			wavetable_[i] = 0.0;
		}
	}
	else {
		/*  SINE WAVE  */
		for (int i = 0; i < TABLE_LENGTH; i++) {
			wavetable_[i] =
			    sin(((double)i / (double)TABLE_LENGTH) * 2.0 * M_PI);
		}
	}

#if OVERSAMPLING_OSCILLATOR
	firFilter_.reset(new FIRFilter(FIR_BETA, FIR_GAMMA, FIR_CUTOFF));
#endif
}

WavetableGlottalSource::~WavetableGlottalSource() {}

void WavetableGlottalSource::reset()
{
	currentPosition_ = 0;
	firFilter_->reset();
}

/******************************************************************************
 *
 *  function:  updateWavetable
 *
 *  purpose:   Rewrites the changeable part of the glottal pulse
 *             according to the amplitude.
 *
 ******************************************************************************/
void WavetableGlottalSource::updateWavetable(double amplitude)
{
	/*  CALCULATE NEW CLOSURE POINT, BASED ON AMPLITUDE  */
	double newDiv2 = tableDiv2_ - rint(amplitude * tnDelta_);
	double invNewTnLength = 1.0 / (newDiv2 - tableDiv1_);

	/*  RECALCULATE THE FALLING PORTION OF THE GLOTTAL PULSE  */
	double x = 0.0;
	for (int i = tableDiv1_, end = static_cast<int>(newDiv2); i < end;
	     ++i, x += invNewTnLength) {
		wavetable_[i] = 1.0 - (x * x);
	}

	/*  FILL IN WITH CLOSED PORTION OF GLOTTAL PULSE  */
	for (int i = static_cast<int>(newDiv2); i < tableDiv2_; i++) {
		wavetable_[i] = 0.0;
	}
}

/******************************************************************************
 *
 *  function:  incrementTablePosition
 *
 *  purpose:   Increments the position in the wavetable according to
 *             the desired frequency.
 *
 ******************************************************************************/
void WavetableGlottalSource::incrementTablePosition(double frequency)
{
	currentPosition_ = mod0(currentPosition_ + (frequency * basicIncrement_));
}

/******************************************************************************
 *
 *  function:  oscillator
 *
 *  purpose:   Is a 2X oversampling interpolating wavetable
 *             oscillator.
 *
 ******************************************************************************/
#if OVERSAMPLING_OSCILLATOR
double WavetableGlottalSource::getSample(
    double frequency) /*  2X OVERSAMPLING OSCILLATOR  */
{
	int lowerPosition, upperPosition;
	double interpolatedValue, output;

	for (int i = 0; i < 2; i++) {
		/*  FIRST INCREMENT THE TABLE POSITION, DEPENDING ON FREQUENCY  */
		incrementTablePosition(frequency / 2.0);

		/*  FIND SURROUNDING INTEGER TABLE POSITIONS  */
		lowerPosition = static_cast<int>(currentPosition_);
		upperPosition = static_cast<int>(mod0(lowerPosition + 1));

		/*  CALCULATE INTERPOLATED TABLE VALUE  */
		interpolatedValue =
		    wavetable_[lowerPosition] +
		    ((currentPosition_ - lowerPosition) *
		     (wavetable_[upperPosition] - wavetable_[lowerPosition]));

		/*  PUT VALUE THROUGH FIR FILTER  */
		output = firFilter_->filter(interpolatedValue, i);
	}

	/*  SINCE WE DECIMATE, TAKE ONLY THE SECOND OUTPUT VALUE  */
	return output;
}
#else
double WavetableGlottalSource::getSample(
    double frequency) /*  PLAIN OSCILLATOR  */
{
	int lowerPosition, upperPosition;

	/*  FIRST INCREMENT THE TABLE POSITION, DEPENDING ON FREQUENCY  */
	incrementTablePosition(frequency);

	/*  FIND SURROUNDING INTEGER TABLE POSITIONS  */
	lowerPosition = (int)currentPosition;
	upperPosition = mod0(lowerPosition + 1);

	/*  RETURN INTERPOLATED TABLE VALUE  */
	return wavetable_[lowerPosition] +
	       ((currentPosition - lowerPosition) *
	        (wavetable_[upperPosition] - wavetable_[lowerPosition]));
}
#endif

/******************************************************************************
 *
 *  function:  mod0
 *
 *  purpose:   Returns the modulus of 'value', keeping it in the
 *             range 0 -> TABLE_MODULUS.
 *
 ******************************************************************************/
double WavetableGlottalSource::mod0(double value)
{
	if (value > TABLE_MODULUS) {
		value -= TABLE_LENGTH;
	}

	return value;
}

} /* namespace TRM */
} /* namespace GS */
