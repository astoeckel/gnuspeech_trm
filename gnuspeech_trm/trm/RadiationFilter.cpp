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

#include "RadiationFilter.h"



namespace GS {
namespace TRM {

RadiationFilter::RadiationFilter(double apertureCoeff)
		: radiationX_(0.0)
		, radiationY_(0.0)
{
	a20_ = apertureCoeff;
	a21_ = b21_ = -a20_;
}

RadiationFilter::~RadiationFilter()
{
}

void
RadiationFilter::reset()
{
	radiationX_ = 0.0;
	radiationY_ = 0.0;
}

double
RadiationFilter::filter(double input)
{
	double output = (a20_ * input) + (a21_ * radiationX_) - (b21_ * radiationY_);
	radiationX_ = input;
	radiationY_ = output;
	return output;
}

} /* namespace TRM */
} /* namespace GS */
