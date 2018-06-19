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

#ifndef TRM_REFLECTION_FILTER_H_
#define TRM_REFLECTION_FILTER_H_



namespace GS {
namespace TRM {

// Is a variable, one-pole lowpass filter, whose cutoff
// is determined by the aperture coefficient.
class ReflectionFilter {
public:
	ReflectionFilter(double apertureCoeff);
	~ReflectionFilter();

	void reset();
	double filter(double input);
private:
	ReflectionFilter(const ReflectionFilter&) = delete;
	ReflectionFilter& operator=(const ReflectionFilter&) = delete;

	double a10_;
	double b11_;
	double reflectionY_;
};

} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_REFLECTION_FILTER_H_ */
