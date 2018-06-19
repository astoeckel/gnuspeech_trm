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

#ifndef TRM_BANDPASS_FILTER_H_
#define TRM_BANDPASS_FILTER_H_

namespace GS {
namespace TRM {

class BandpassFilter {
public:
	BandpassFilter();
	~BandpassFilter();

	void reset();
	void update(double sampleRate, double bandwidth, double centerFreq);
	double filter(double input);

private:
	BandpassFilter(const BandpassFilter &) = delete;
	BandpassFilter &operator=(const BandpassFilter &) = delete;

	double bpAlpha_;
	double bpBeta_;
	double bpGamma_;
	double xn1_;
	double xn2_;
	double yn1_;
	double yn2_;
};

} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_BANDPASS_FILTER_H_ */
