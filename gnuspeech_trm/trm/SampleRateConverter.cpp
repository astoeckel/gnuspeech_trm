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

#include "SampleRateConverter.h"

#include <cmath>

#define BETA                      5.658        /*  kaiser window parameters  */
#define IzeroEPSILON              1E-21

/*  SAMPLE RATE CONVERSION CONSTANTS  */
#define ZERO_CROSSINGS            13                 /*  SRC CUTOFF FRQ      */
#define LP_CUTOFF                 (11.0/13.0)        /*  (0.846 OF NYQUIST)  */
#define FILTER_LENGTH             (ZERO_CROSSINGS * L_RANGE)

//#define N_BITS                    16
#define L_BITS                    8
#define L_RANGE                   256                  /*  must be 2^L_BITS  */
#define M_BITS                    8
#define M_RANGE                   256                  /*  must be 2^M_BITS  */
#define FRACTION_BITS             (L_BITS + M_BITS)
#define FRACTION_RANGE            65536         /*  must be 2^FRACTION_BITS  */
#define FILTER_LIMIT              (FILTER_LENGTH - 1)

#define N_MASK                    0xFFFF0000
#define L_MASK                    0x0000FF00
#define M_MASK                    0x000000FF
#define FRACTION_MASK             0x0000FFFF

#define nValue(x)                 (((x) & N_MASK) >> FRACTION_BITS)
#define lValue(x)                 (((x) & L_MASK) >> M_BITS)
#define mValue(x)                 ((x) & M_MASK)
#define fractionValue(x)          ((x) & FRACTION_MASK)

#define BUFFER_SIZE               1024                 /*  ring buffer size  */



namespace GS {
namespace TRM {

SampleRateConverter::SampleRateConverter(int sampleRate, float outputRate, std::vector<float>& outputData)
		: sampleRateRatio_(0.0)
		, fillPtr_(0)
		, emptyPtr_(0)
		, padSize_(0)
		, fillSize_(0)
		, timeRegisterIncrement_(0)
		, filterIncrement_(0)
		, phaseIncrement_(0)
		, timeRegister_(0)
		, fillCounter_(0)
		, maximumSampleValue_(0.0)
		, numberSamples_(0)
		, h_(FILTER_LENGTH)
		, deltaH_(FILTER_LENGTH)
		, buffer_(BUFFER_SIZE)
		, outputData_(outputData)
{
	initializeConversion(sampleRate, outputRate);
}

SampleRateConverter::~SampleRateConverter()
{
}

void
SampleRateConverter::reset()
{
	emptyPtr_ = 0;
	timeRegister_ = 0;
	fillCounter_ = 0;
	maximumSampleValue_ = 0.0;
	numberSamples_ = 0;
	initializeBuffer();
}

/******************************************************************************
*
*  function:  initializeConversion
*
*  purpose:   Initializes all the sample rate conversion functions.
*
******************************************************************************/
void
SampleRateConverter::initializeConversion(int sampleRate, float outputRate)
{
	/*  INITIALIZE FILTER IMPULSE RESPONSE  */
	initializeFilter();

	/*  CALCULATE SAMPLE RATE RATIO  */
	sampleRateRatio_ = (double) outputRate / (double) sampleRate;

	/*  CALCULATE TIME REGISTER INCREMENT  */
	timeRegisterIncrement_ =
			(int) rint(pow(2.0, FRACTION_BITS) / sampleRateRatio_);

	/*  CALCULATE ROUNDED SAMPLE RATE RATIO  */
	double roundedSampleRateRatio =
			pow(2.0, FRACTION_BITS) / (double) timeRegisterIncrement_;

	/*  CALCULATE PHASE OR FILTER INCREMENT  */
	if (sampleRateRatio_ >= 1.0) {
		filterIncrement_ = L_RANGE;
	} else {
		phaseIncrement_ = (unsigned int) rint(sampleRateRatio_ * (double) FRACTION_RANGE);
	}

	/*  CALCULATE PAD SIZE  */
	padSize_ = (sampleRateRatio_ >= 1.0) ?
				ZERO_CROSSINGS :
				(int) ((float) ZERO_CROSSINGS / roundedSampleRateRatio) + 1;

	/*  INITIALIZE THE RING BUFFER  */
	initializeBuffer();
}

/******************************************************************************
*
*  function:  Izero
*
*  purpose:   Returns the value for the modified Bessel function of
*             the first kind, order 0, as a double.
*
******************************************************************************/
double
SampleRateConverter::Izero(double x)
{
	double sum, u, halfx, temp;
	int n;

	sum = u = n = 1;
	halfx = x / 2.0;

	do {
		temp = halfx / (double) n;
		n += 1;
		temp *= temp;
		u *= temp;
		sum += u;
	} while (u >= (IzeroEPSILON * sum));

	return sum;
}

/******************************************************************************
*
* function:  initializeBuffer
*
*  purpose:  Initializes the ring buffer used for sample rate
*            conversion.
*
******************************************************************************/
void
SampleRateConverter::initializeBuffer()
{
	/*  FILL THE RING BUFFER WITH ALL ZEROS  */
	for (int i = 0; i < BUFFER_SIZE; i++) {
		buffer_[i] = 0.0;
	}

	/*  INITIALIZE FILL POINTER  */
	fillPtr_ = padSize_;

	/*  CALCULATE FILL SIZE  */
	fillSize_ = BUFFER_SIZE - (2 * padSize_);
}

/******************************************************************************
*
*  function:  initializeFilter
*
*  purpose:   Initializes filter impulse response and impulse delta
*             values.
*
******************************************************************************/
void
SampleRateConverter::initializeFilter()
{
	/*  INITIALIZE THE FILTER IMPULSE RESPONSE  */
	h_[0] = LP_CUTOFF;
	double x = M_PI / (double) L_RANGE;
	for (int i = 1; i < FILTER_LENGTH; i++) {
		double y = (double) i * x;
		h_[i] = sin(y * LP_CUTOFF) / y;
	}

	/*  APPLY A KAISER WINDOW TO THE IMPULSE RESPONSE  */
	double IBeta = 1.0 / Izero(BETA);
	for (int i = 0; i < FILTER_LENGTH; i++) {
		double temp = (double) i / FILTER_LENGTH;
		h_[i] *= Izero(BETA * sqrt(1.0 - (temp * temp))) * IBeta;
	}

	/*  INITIALIZE THE FILTER IMPULSE RESPONSE DELTA VALUES  */
	for (int i = 0; i < FILTER_LIMIT; i++) {
		deltaH_[i] = h_[i + 1] - h_[i];
	}
	deltaH_[FILTER_LIMIT] = 0.0 - h_[FILTER_LIMIT];
}

/******************************************************************************
*
*  function:  dataFill
*
*  purpose:   Fills the ring buffer with a single sample, increments
*             the counters and pointers, and empties the buffer when
*             full.
*
******************************************************************************/
void
SampleRateConverter::dataFill(double data)
{
	/*  PUT THE DATA INTO THE RING BUFFER  */
	buffer_[fillPtr_] = data;

	/*  INCREMENT THE FILL POINTER, MODULO THE BUFFER SIZE  */
	srIncrement(&fillPtr_, BUFFER_SIZE);

	/*  INCREMENT THE COUNTER, AND EMPTY THE BUFFER IF FULL  */
	if (++fillCounter_ >= fillSize_) {
		dataEmpty();
		/* RESET THE FILL COUNTER  */
		fillCounter_ = 0;
	}
}

/******************************************************************************
*
*  function:  dataEmpty
*
*  purpose:   Converts available portion of the input signal to the
*             new sampling rate, and outputs the samples to the
*             sound struct.
*
******************************************************************************/
void
SampleRateConverter::dataEmpty()
{
	/*  CALCULATE END POINTER  */
	int endPtr = fillPtr_ - padSize_;

	/*  ADJUST THE END POINTER, IF LESS THAN ZERO  */
	if (endPtr < 0) {
		endPtr += BUFFER_SIZE;
	}

	/*  ADJUST THE ENDPOINT, IF LESS THEN THE EMPTY POINTER  */
	if (endPtr < emptyPtr_) {
		endPtr += BUFFER_SIZE;
	}

	/*  UPSAMPLE LOOP (SLIGHTLY MORE EFFICIENT THAN DOWNSAMPLING)  */
	if (sampleRateRatio_ >= 1.0) {
		while (emptyPtr_ < endPtr) {
			/*  RESET ACCUMULATOR TO ZERO  */
			double output = 0.0;

			/*  CALCULATE INTERPOLATION VALUE (STATIC WHEN UPSAMPLING)  */
			double interpolation = (double) mValue(timeRegister_) / (double) M_RANGE;

			/*  COMPUTE THE LEFT SIDE OF THE FILTER CONVOLUTION  */
			int index = emptyPtr_;
			for (int filterIndex = lValue(timeRegister_);
					filterIndex < FILTER_LENGTH;
					srDecrement(&index,BUFFER_SIZE), filterIndex += filterIncrement_) {
				output += (buffer_[index] *
						(h_[filterIndex] + (deltaH_[filterIndex] * interpolation)));
			}

			/*  ADJUST VALUES FOR RIGHT SIDE CALCULATION  */
			timeRegister_ = ~timeRegister_;
			interpolation = (double) mValue(timeRegister_) / (double) M_RANGE;

			/*  COMPUTE THE RIGHT SIDE OF THE FILTER CONVOLUTION  */
			index = emptyPtr_;
			srIncrement(&index,BUFFER_SIZE);
			for (int filterIndex = lValue(timeRegister_);
					filterIndex < FILTER_LENGTH;
					srIncrement(&index,BUFFER_SIZE), filterIndex += filterIncrement_) {
				output += (buffer_[index] *
						(h_[filterIndex] + (deltaH_[filterIndex] * interpolation)));
			}

			/*  RECORD MAXIMUM SAMPLE VALUE  */
			double absoluteSampleValue = fabs(output);
			if (absoluteSampleValue > maximumSampleValue_) {
				maximumSampleValue_ = absoluteSampleValue;
			}

			/*  INCREMENT SAMPLE NUMBER  */
			numberSamples_++;

			/*  SAVE THE SAMPLE  */
			outputData_.push_back(static_cast<float>(output));

			/*  CHANGE TIME REGISTER BACK TO ORIGINAL FORM  */
			timeRegister_ = ~timeRegister_;

			/*  INCREMENT THE TIME REGISTER  */
			timeRegister_ += timeRegisterIncrement_;

			/*  INCREMENT THE EMPTY POINTER, ADJUSTING IT AND END POINTER  */
			emptyPtr_ += nValue(timeRegister_);

			if (emptyPtr_ >= BUFFER_SIZE) {
				emptyPtr_ -= BUFFER_SIZE;
				endPtr -= BUFFER_SIZE;
			}

			/*  CLEAR N PART OF TIME REGISTER  */
			timeRegister_ &= (~N_MASK);
		}
	} else {
		/*  DOWNSAMPLING CONVERSION LOOP  */

		while (emptyPtr_ < endPtr) {

			/*  RESET ACCUMULATOR TO ZERO  */
			double output = 0.0;

			/*  COMPUTE P PRIME  */
			unsigned int phaseIndex = (unsigned int) rint(
						((double) fractionValue(timeRegister_)) * sampleRateRatio_);

			/*  COMPUTE THE LEFT SIDE OF THE FILTER CONVOLUTION  */
			int index = emptyPtr_;
			unsigned int impulseIndex;
			while ((impulseIndex = (phaseIndex >> M_BITS)) < FILTER_LENGTH) {
				double impulse = h_[impulseIndex] + (deltaH_[impulseIndex] *
						(((double) mValue(phaseIndex)) / (double) M_RANGE));
				output += (buffer_[index] * impulse);
				srDecrement(&index, BUFFER_SIZE);
				phaseIndex += phaseIncrement_;
			}

			/*  COMPUTE P PRIME, ADJUSTED FOR RIGHT SIDE  */
			phaseIndex = (unsigned int) rint(
						((double) fractionValue(~timeRegister_)) * sampleRateRatio_);

			/*  COMPUTE THE RIGHT SIDE OF THE FILTER CONVOLUTION  */
			index = emptyPtr_;
			srIncrement(&index, BUFFER_SIZE);
			while ((impulseIndex = (phaseIndex >> M_BITS)) < FILTER_LENGTH) {
				double impulse = h_[impulseIndex] + (deltaH_[impulseIndex] *
						(((double) mValue(phaseIndex)) / (double) M_RANGE));
				output += (buffer_[index] * impulse);
				srIncrement(&index, BUFFER_SIZE);
				phaseIndex += phaseIncrement_;
			}

			/*  RECORD MAXIMUM SAMPLE VALUE  */
			double absoluteSampleValue = fabs(output);
			if (absoluteSampleValue > maximumSampleValue_) {
				maximumSampleValue_ = absoluteSampleValue;
			}

			/*  INCREMENT SAMPLE NUMBER  */
			numberSamples_++;

			/*  SAVE THE SAMPLE  */
			outputData_.push_back(static_cast<float>(output));

			/*  INCREMENT THE TIME REGISTER  */
			timeRegister_ += timeRegisterIncrement_;

			/*  INCREMENT THE EMPTY POINTER, ADJUSTING IT AND END POINTER  */
			emptyPtr_ += nValue(timeRegister_);
			if (emptyPtr_ >= BUFFER_SIZE) {
				emptyPtr_ -= BUFFER_SIZE;
				endPtr -= BUFFER_SIZE;
			}

			/*  CLEAR N PART OF TIME REGISTER  */
			timeRegister_ &= (~N_MASK);
		}
	}
}

/******************************************************************************
*
*  function:  srIncrement
*
*  purpose:   Increments the pointer, keeping it within the range
*             0 to (modulus-1).
*
******************************************************************************/
void
SampleRateConverter::srIncrement(int *pointer, int modulus)
{
	if (++(*pointer) >= modulus) {
		(*pointer) -= modulus;
	}
}

/******************************************************************************
*
*  function:  srDecrement
*
*  purpose:   Decrements the pointer, keeping it within the range
*             0 to (modulus-1).
*
******************************************************************************/
void
SampleRateConverter::srDecrement(int *pointer, int modulus)
{
	if (--(*pointer) < 0) {
		(*pointer) += modulus;
	}
}

/******************************************************************************
*
*  function:  flushBuffer
*
*  purpose:   Pads the buffer with zero samples, and flushes it by
*             converting the remaining samples.
*
******************************************************************************/
void
SampleRateConverter::flushBuffer()
{
	/*  PAD END OF RING BUFFER WITH ZEROS  */
	for (int i = 0; i < padSize_ * 2; i++) {
		dataFill(0.0);
	}

	/*  FLUSH UP TO FILL POINTER - PADSIZE  */
	dataEmpty();
}

} /* namespace TRM */
} /* namespace GS */
