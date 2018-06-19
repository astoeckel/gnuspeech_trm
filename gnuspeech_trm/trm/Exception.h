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

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <algorithm> /* move */
#include <cassert>
#include <cstdio>    /* fprintf */
#include <cstdlib>   /* free, malloc */
#include <cstring>   /* strcpy, strlen */
#include <exception>
#include <sstream>
#include <string>

// MSVC 2013
#if defined(_MSC_VER) && defined(_NOEXCEPT)
# define GS_EXCEPTION_NOEXCEPT _NOEXCEPT
#else
# define GS_EXCEPTION_NOEXCEPT noexcept
#endif

// __func__ is defined in C99/C++11.
// __PRETTY_FUNCTION__ is a gcc extension.
#ifdef __GNUC__
# define GS_EXCEPTION_FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined (_MSC_VER)
# define GS_EXCEPTION_FUNCTION_NAME __FUNCTION__
#else
# define GS_EXCEPTION_FUNCTION_NAME __func__
#endif

#define THROW_EXCEPTION(E,M) \
	do {\
		E exc;\
		try { \
			GS::ErrorMessage em;\
			em << M << "\n[file: " << __FILE__ << "]\n[function: " << GS_EXCEPTION_FUNCTION_NAME << "]\n[line: " << __LINE__ << "]";\
			exc.setMessage(em);\
		} catch (...) {}\
		throw exc;\
	} while (false)



namespace GS {

// Note: string / vector default constructor may throw bad_alloc in C++11.
class ExceptionString {
public:
	ExceptionString() GS_EXCEPTION_NOEXCEPT : str_(nullptr) {}
	ExceptionString(const ExceptionString& o) GS_EXCEPTION_NOEXCEPT : str_(nullptr) {
		*this = o;
	}
	ExceptionString(ExceptionString&& o) GS_EXCEPTION_NOEXCEPT : str_(nullptr) {
		*this = std::move(o);
	}
	~ExceptionString() GS_EXCEPTION_NOEXCEPT {
		std::free(str_);
	}
	ExceptionString& operator=(const ExceptionString& o) GS_EXCEPTION_NOEXCEPT {
		if (this != &o) {
			if (o.str_ == nullptr) {
				std::free(str_);
				str_ = nullptr;
				return *this;
			}
			std::size_t size = std::strlen(o.str_);
			auto p = static_cast<char*>(std::malloc(size + 1));
			if (p == nullptr) {
				std::fprintf(stderr, "Exception string copy error. String: %s\n", o.str_);
				return *this;
			}
			std::free(str_);
			str_ = p;
			std::strcpy(str_, o.str_);
		}
		return *this;
	}
	ExceptionString& operator=(ExceptionString&& o) GS_EXCEPTION_NOEXCEPT {
		assert(this != &o);
		std::free(str_);
		str_ = o.str_;
		o.str_ = nullptr;
		return *this;
	}
	const char* str() const GS_EXCEPTION_NOEXCEPT {
		return str_ ? str_ : "";
	}
	void setStr(const char* s) GS_EXCEPTION_NOEXCEPT {
		if (s == nullptr) {
			std::free(str_);
			str_ = nullptr;
			return;
		}
		std::size_t size = std::strlen(s);
		auto p = static_cast<char*>(std::malloc(size + 1));
		if (p == nullptr) {
			std::fprintf(stderr, "Exception string assignment error. String: %s\n", s);
			return;
		}
		std::free(str_);
		str_ = p;
		std::strcpy(str_, s);
	}
private:
	char* str_;
};

/*******************************************************************************
 *
 * This class may throw std::bad_alloc.
 */
class ErrorMessage {
public:
	ErrorMessage() {}
	~ErrorMessage() {}

	template<typename T>
	ErrorMessage& operator<<(const T& messagePart) {
		buffer_ << messagePart;
		return *this;
	}

	ErrorMessage& operator<<(const std::exception& e) {
		buffer_ << e.what();
		return *this;
	}

	std::string getString() const {
		return buffer_.str();
	}
private:
	ErrorMessage(const ErrorMessage&) = delete;
	ErrorMessage& operator=(const ErrorMessage&) = delete;

	std::ostringstream buffer_;
};

/*******************************************************************************
 *
 */
class Exception : public std::exception {
public:
	virtual const char* what() const GS_EXCEPTION_NOEXCEPT {
		return message_.str();
	}

	// May throw std::bad_alloc.
	void setMessage(const ErrorMessage& em) {
		std::string msg = em.getString();
		message_.setStr(msg.c_str());
	}
protected:
	ExceptionString message_;
};

class EndOfBufferException              : public Exception {};
class ExternalProgramExecutionException : public Exception {};
class InvalidCallException              : public Exception {};
class InvalidDirectoryException         : public Exception {};
class InvalidFileException              : public Exception {};
class InvalidParameterException         : public Exception {};
class InvalidStateException             : public Exception {};
class InvalidValueException             : public Exception {};
class IOException                       : public Exception {};
class MissingValueException             : public Exception {};
class ParsingException                  : public Exception {};
class TextParserException               : public Exception {};
class TRMControlModelException          : public Exception {};
class TRMException                      : public Exception {};
class UnavailableResourceException      : public Exception {};
class WrongBufferSizeException          : public Exception {};
class XMLException                      : public Exception {};

} /* namespace GS */

#endif /* EXCEPTION_H_ */
