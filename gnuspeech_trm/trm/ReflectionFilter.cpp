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

#include "ReflectionFilter.h"

#include <cmath>



namespace GS {
namespace TRM {

ReflectionFilter::ReflectionFilter(double apertureCoeff)
		: reflectionY_(0.0)
{
	b11_ = -apertureCoeff;
	a10_ = 1.0 - fabs(b11_);
}

ReflectionFilter::~ReflectionFilter()
{
}

void
ReflectionFilter::reset()
{
	reflectionY_ = 0.0;
}

double
ReflectionFilter::filter(double input)
{
	double output = (a10_ * input) - (b11_ * reflectionY_);
	reflectionY_ = output;
	return output;
}

} /* namespace TRM */
} /* namespace GS */
