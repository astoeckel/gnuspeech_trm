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

#ifndef TRM_FIR_FILTER_H_
#define TRM_FIR_FILTER_H_

#include <vector>

namespace GS {
namespace TRM {

/******************************************************************************
 *
 *  class:    FIRFilter
 *
 *  purpose:  Is the linear phase, lowpass FIR filter.
 *
 ******************************************************************************/
class FIRFilter {
public:
	FIRFilter(double beta, double gamma, double cutoff);
	~FIRFilter();

	void reset();
	double filter(double input, int needOutput);

private:
	FIRFilter(const FIRFilter &) = delete;
	FIRFilter &operator=(const FIRFilter &) = delete;

	static int maximallyFlat(double beta, double gamma, int *np,
	                         double *coefficient);
	static void trim(double cutoff, int *numberCoefficients,
	                 double *coefficient);
	static int increment(int pointer, int modulus);
	static int decrement(int pointer, int modulus);
	static void rationalApproximation(double number, int *order, int *numerator,
	                                  int *denominator);

	std::vector<double> data_;
	std::vector<double> coef_;
	int ptr_;
	int numberTaps_;
};

} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_FIR_FILTER_H_ */
