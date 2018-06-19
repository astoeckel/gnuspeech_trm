/***************************************************************************
 *  Copyright 2014 Marcelo Y. Matuda                                       *
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

#ifndef MOVING_AVERAGE_FILTER_H_
#define MOVING_AVERAGE_FILTER_H_

#include <cassert>
#include <cmath>   /* round */
#include <cstddef> /* std::size_t */
#include <vector>

namespace GS {
namespace TRM {

/*******************************************************************************
 * Moving average filter.
 */
template <typename FloatType>
class MovingAverageFilter {
public:
	MovingAverageFilter(FloatType sampleRate, FloatType period /* seconds */);

	void reset();
	FloatType filter(FloatType value);

private:
	std::vector<FloatType> buf_;
	typename std::vector<FloatType>::size_type pos_;
	double sum_;
	double invN_;
};

/*******************************************************************************
 * Constructor.
 */
template <typename FloatType>
MovingAverageFilter<FloatType>::MovingAverageFilter(FloatType sampleRate,
                                                    FloatType period)
    : buf_(static_cast<std::size_t>(std::round(sampleRate * period))),
      pos_(buf_.size()),
      sum_(0.0),
      invN_(1.0 / buf_.size())
{
	assert(!buf_.empty());
}

/*******************************************************************************
 *
 */
template <typename FloatType>
void MovingAverageFilter<FloatType>::reset()
{
	for (auto &item : buf_)
		item = 0.0;
	pos_ = buf_.size();
	sum_ = 0.0;
}

/*******************************************************************************
 *
 */
template <typename FloatType>
FloatType MovingAverageFilter<FloatType>::filter(FloatType value)
{
	if (++pos_ >= buf_.size()) {
		if (pos_ > buf_.size()) {  // first value
			buf_.assign(buf_.size(), value);
			sum_ = value * static_cast<double>(buf_.size());
		}
		pos_ = 0;
	}
	sum_ -= buf_[pos_];
	sum_ += value;
	buf_[pos_] = value;
	return sum_ * invN_;
}

} /* namespace TRM */
} /* namespace GS */

#endif /* MOVING_AVERAGE_FILTER_H_ */
