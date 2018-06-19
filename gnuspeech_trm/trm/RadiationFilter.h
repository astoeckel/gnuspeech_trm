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

#ifndef TRM_RADIATION_FILTER_H_
#define TRM_RADIATION_FILTER_H_



namespace GS {
namespace TRM {

// Is a variable, one-zero, one-pole, highpass filter,
// whose cutoff point is determined by the aperture
// coefficient.
class RadiationFilter {
public:
	RadiationFilter(double apertureCoeff);
	~RadiationFilter();

	void reset();
	double filter(double input);
private:
	RadiationFilter(const RadiationFilter&) = delete;
	RadiationFilter& operator=(const RadiationFilter&) = delete;

	double a20_;
	double a21_;
	double b21_;
	double radiationX_;
	double radiationY_;
};

} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_RADIATION_FILTER_H_ */
