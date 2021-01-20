/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef GT2_EXCEPTION_H
#define GT2_EXCEPTION_H

#include <exception>
#include <string>

#include "macros.h"

namespace GeneTrail
{
	class GT2_EXPORT IOError : public std::exception
	{
		public:
			IOError(const std::string& error) noexcept;
			virtual ~IOError() {}

			virtual const char* what() const noexcept override;

		private:
			std::string msg_;
	};

	class GT2_EXPORT NotImplemented : public std::exception
	{
		public:
			NotImplemented(const char* file, int line, const std::string& method) noexcept;
			const char* what() const noexcept override;

		private:
			std::string msg_;
	};

	class GT2_EXPORT InvalidIndex : public std::exception
	{
		public:
			InvalidIndex(unsigned int i, unsigned int max) noexcept;
			const char* what() const noexcept override;

		private:
			std::string msg_;
	};

	class GT2_EXPORT InvalidKey : public std::exception
	{
	  public:
		InvalidKey(const std::string& str) noexcept;
		const char* what() const noexcept override;

	  private:
		std::string msg_;
	};

	class GT2_EXPORT UnknownEntry : public std::exception
	{
		public:
		UnknownEntry(const std::string& entry_as_string) noexcept;
		const char* what() const noexcept override;

		private:
		std::string msg_;
	};
}

#endif // GT2_EXCEPTION_H

