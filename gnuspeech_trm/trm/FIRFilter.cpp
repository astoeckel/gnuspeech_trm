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
#include <stdexcept>

#include "FIRFilter.h"

#define LIMIT 200

namespace GS {
namespace TRM {

FIRFilter::FIRFilter(double beta, double gamma, double cutoff)
{
	int pointer, increment, numberCoefficients;
	double coefficient[LIMIT + 1];

	/*  DETERMINE IDEAL LOW PASS FILTER COEFFICIENTS  */
	maximallyFlat(beta, gamma, &numberCoefficients, coefficient);

	/*  TRIM LOW-VALUE COEFFICIENTS  */
	trim(cutoff, &numberCoefficients, coefficient);

	/*  DETERMINE THE NUMBER OF TAPS IN THE FILTER  */
	numberTaps_ = (numberCoefficients * 2) - 1;

	/*  ALLOCATE MEMORY FOR DATA AND COEFFICIENTS  */
	data_.resize(numberTaps_);
	coef_.resize(numberTaps_);

	/*  INITIALIZE THE COEFFICIENTS  */
	increment = -1;
	pointer = numberCoefficients;
	for (int i = 0; i < numberTaps_; i++) {
		coef_[i] = coefficient[pointer];
		pointer += increment;
		if (pointer <= 0) {
			pointer = 2;
			increment = 1;
		}
	}

	/*  SET POINTER TO FIRST ELEMENT  */
	ptr_ = 0;

#if 0
	/*  PRINT OUT  */
	printf("\n");
	for (int i = 0; i < numberTaps_; i++) {
		printf("coef_[%-d] = %11.8f\n", i, coef_[i]);
	}
#endif
}

FIRFilter::~FIRFilter() {}

void FIRFilter::reset()
{
	for (auto &item : data_)
		item = 0.0;
	ptr_ = 0;
}

/******************************************************************************
 *
 *  function:  maximallyFlat
 *
 *  purpose:   Calculates coefficients for a linear phase lowpass FIR
 *             filter, with beta being the center frequency of the
 *             transition band (as a fraction of the sampling
 *             frequency), and gamme the width of the transition
 *             band.
 *
 ******************************************************************************/
int FIRFilter::maximallyFlat(double beta, double gamma, int *np,
                             double *coefficient)
{
	double a[LIMIT + 1], c[LIMIT + 1], betaMinimum, ac;
	int nt, numerator, n, ll, i;

	/*  INITIALIZE NUMBER OF POINTS  */
	*np = 0;

	/*  CUT-OFF FREQUENCY MUST BE BETWEEN 0 HZ AND NYQUIST  */
	if ((beta <= 0.0) || (beta >= 0.5)) {
		throw std::runtime_error("Beta out of range.");
	}

	/*  TRANSITION BAND MUST FIT WITH THE STOP BAND  */
	betaMinimum =
	    ((2.0 * beta) < (1.0 - 2.0 * beta)) ? (2.0 * beta) : (1.0 - 2.0 * beta);
	if ((gamma <= 0.0) || (gamma >= betaMinimum)) {
		throw std::runtime_error("Gamma out of range.");
	}

	/*  MAKE SURE TRANSITION BAND NOT TOO SMALL  */
	nt = (int)(1.0 / (4.0 * gamma * gamma));
	if (nt > 160) {
		throw std::runtime_error("Gamma too small.");
	}

	/*  CALCULATE THE RATIONAL APPROXIMATION TO THE CUT-OFF POINT  */
	ac = (1.0 + cos((2.0 * M_PI) * beta)) / 2.0;
	rationalApproximation(ac, &nt, &numerator, np);

	/*  CALCULATE FILTER ORDER  */
	n = (2 * (*np)) - 1;
	if (numerator == 0) {
		numerator = 1;
	}

	/*  COMPUTE MAGNITUDE AT NP POINTS  */
	c[1] = a[1] = 1.0;
	ll = nt - numerator;

	for (i = 2; i <= *np; i++) {
		int j;
		double x, sum = 1.0, y;
		c[i] = cos((2.0 * M_PI) * ((double)(i - 1) / (double)n));
		x = (1.0 - c[i]) / 2.0;
		y = x;

		if (numerator == nt) {
			continue;
		}

		for (j = 1; j <= ll; j++) {
			double z = y;
			if (numerator != 1) {
				int jj;
				for (jj = 1; jj <= (numerator - 1); jj++) {
					z *= 1.0 + ((double)j / (double)jj);
				}
			}
			y *= x;
			sum += z;
		}
		a[i] = sum * pow((1.0 - x), numerator);
	}

	/*  CALCULATE WEIGHTING COEFFICIENTS BY AN N-POINT IDFT  */
	for (i = 1; i <= *np; i++) {
		int j;
		coefficient[i] = a[1] / 2.0;
		for (j = 2; j <= *np; j++) {
			int m = ((i - 1) * (j - 1)) % n;
			if (m > nt) {
				m = n - m;
			}
			coefficient[i] += c[m + 1] * a[j];
		}
		coefficient[i] *= 2.0 / (double)n;
	}

	return 0;
}

/******************************************************************************
 *
 *  function:  trim
 *
 *  purpose:   Trims the higher order coefficients of the FIR filter
 *             which fall below the cutoff value.
 *
 ******************************************************************************/
void FIRFilter::trim(double cutoff, int *numberCoefficients,
                     double *coefficient)
{
	for (int i = *numberCoefficients; i > 0; i--) {
		if (fabs(coefficient[i]) >= fabs(cutoff)) {
			*numberCoefficients = i;
			return;
		}
	}
}

/******************************************************************************
 *
 *  function:  increment
 *
 *  purpose:   Increments the pointer to the circular FIR filter
 *             buffer, keeping it in the range 0 -> modulus-1.
 *
 ******************************************************************************/
int FIRFilter::increment(int pointer, int modulus)
{
	if (++pointer >= modulus) {
		return 0;
	}
	else {
		return pointer;
	}
}

/******************************************************************************
 *
 *  function:  decrement
 *
 *  purpose:   Decrements the pointer to the circular FIR filter
 *             buffer, keeping it in the range 0 -> modulus-1.
 *
 ******************************************************************************/
int FIRFilter::decrement(int pointer, int modulus)
{
	if (--pointer < 0) {
		return modulus - 1;
	}
	else {
		return pointer;
	}
}

double FIRFilter::filter(double input, int needOutput)
{
	if (needOutput) {
		int i;
		double output = 0.0;

		/*  PUT INPUT SAMPLE INTO DATA BUFFER  */
		data_[ptr_] = input;

		/*  SUM THE OUTPUT FROM ALL FILTER TAPS  */
		for (i = 0; i < numberTaps_; i++) {
			output += data_[ptr_] * coef_[i];
			ptr_ = increment(ptr_, numberTaps_);
		}

		/*  DECREMENT THE DATA POINTER READY FOR NEXT CALL  */
		ptr_ = decrement(ptr_, numberTaps_);

		/*  RETURN THE OUTPUT VALUE  */
		return output;
	}
	else {
		/*  PUT INPUT SAMPLE INTO DATA BUFFER  */
		data_[ptr_] = input;

		/*  ADJUST THE DATA POINTER, READY FOR NEXT CALL  */
		ptr_ = decrement(ptr_, numberTaps_);

		return 0.0;
	}
}

/******************************************************************************
 *
 *  function:  rationalApproximation
 *
 *  purpose:   Calculates the best rational approximation to 'number',
 *             given the maximum 'order'.
 *
 ******************************************************************************/
void FIRFilter::rationalApproximation(double number, int *order, int *numerator,
                                      int *denominator)
{
	double fractionalPart, minimumError = 1.0;
	int i, orderMaximum, modulus = 0;

	/*  RETURN IMMEDIATELY IF THE ORDER IS LESS THAN ONE  */
	if (*order <= 0) {
		*numerator = 0;
		*denominator = 0;
		*order = -1;
		return;
	}

	/*  FIND THE ABSOLUTE VALUE OF THE FRACTIONAL PART OF THE NUMBER  */
	fractionalPart = fabs(number - (int)number);

	/*  DETERMINE THE MAXIMUM VALUE OF THE DENOMINATOR  */
	orderMaximum = 2 * (*order);
	orderMaximum = (orderMaximum > LIMIT) ? LIMIT : orderMaximum;

	/*  FIND THE BEST DENOMINATOR VALUE  */
	for (i = (*order); i <= orderMaximum; i++) {
		double ps = i * fractionalPart;
		int ip = (int)(ps + 0.5);
		double error = fabs((ps - (double)ip) / (double)i);
		if (error < minimumError) {
			minimumError = error;
			modulus = ip;
			*denominator = i;
		}
	}

	/*  DETERMINE THE NUMERATOR VALUE, MAKING IT NEGATIVE IF NECESSARY  */
	*numerator = (int)fabs(number) * (*denominator) + modulus;
	if (number < 0) {
		*numerator *= -1;
	}

	/*  SET THE ORDER  */
	*order = *denominator - 1;

	/*  RESET THE NUMERATOR AND DENOMINATOR IF THEY ARE EQUAL  */
	if (*numerator == *denominator) {
		*denominator = orderMaximum;
		*order = *numerator = *denominator - 1;
	}
}

} /* namespace TRM */
} /* namespace GS */
