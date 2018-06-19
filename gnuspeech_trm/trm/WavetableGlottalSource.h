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

#ifndef TRM_WAVETABLE_GLOTTAL_SOURCE_H_
#define TRM_WAVETABLE_GLOTTAL_SOURCE_H_

#include <memory>
#include <vector>

namespace GS {
namespace TRM {

class FIRFilter;

class WavetableGlottalSource {
public:
	enum Type { /*  WAVEFORM TYPES  */
		        TYPE_PULSE,
		        TYPE_SINE
	};

	WavetableGlottalSource(Type type, double sampleRate, double tp = 0.0,
	                       double tnMin = 0.0, double tnMax = 0.0);
	~WavetableGlottalSource();

	void reset();
	double getSample(double frequency);
	void updateWavetable(double amplitude);

private:
	WavetableGlottalSource(const WavetableGlottalSource &) = delete;
	WavetableGlottalSource &operator=(const WavetableGlottalSource &) = delete;

	void incrementTablePosition(double frequency);

	static double mod0(double value);

	int tableDiv1_;
	int tableDiv2_;
	double tnLength_;
	double tnDelta_;
	double basicIncrement_;
	double currentPosition_;
	std::vector<double> wavetable_;
	std::unique_ptr<FIRFilter> firFilter_;
};

} /* namespace TRM */
} /* namespace GS */

#endif /* TRM_WAVETABLE_GLOTTAL_SOURCE_H_ */
