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

#include "BandpassFilter.h"

namespace GS {
namespace TRM {

BandpassFilter::BandpassFilter()
    : bpAlpha_(0.0),
      bpBeta_(0.0),
      bpGamma_(0.0),
      xn1_(0.0),
      xn2_(0.0),
      yn1_(0.0),
      yn2_(0.0)
{
}

BandpassFilter::~BandpassFilter() {}

void BandpassFilter::reset()
{
	xn1_ = 0.0;
	xn2_ = 0.0;
	yn1_ = 0.0;
	yn2_ = 0.0;
}

void BandpassFilter::update(double sampleRate, double bandwidth,
                            double centerFreq)
{
	double tanValue = tan((M_PI * bandwidth) / sampleRate);
	double cosValue = cos((2.0 * M_PI * centerFreq) / sampleRate);
	bpBeta_ = (1.0 - tanValue) / (2.0 * (1.0 + tanValue));
	bpGamma_ = (0.5 + bpBeta_) * cosValue;
	bpAlpha_ = (0.5 - bpBeta_) / 2.0;
}

/******************************************************************************
 *
 *  function:  bandpassFilter
 *
 *  purpose:   Frication bandpass filter, with variable center
 *             frequency and bandwidth.
 *
 ******************************************************************************/
double BandpassFilter::filter(double input)
{
	double output = 2.0 * ((bpAlpha_ * (input - xn2_)) + (bpGamma_ * yn1_) -
	                       (bpBeta_ * yn2_));

	xn2_ = xn1_;
	xn1_ = input;
	yn2_ = yn1_;
	yn1_ = output;

	return output;
}

} /* namespace TRM */
} /* namespace GS */
