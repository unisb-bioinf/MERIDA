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

#include "Exception.h"

#include <boost/lexical_cast.hpp>

namespace GeneTrail
{
	IOError::IOError(const std::string& error) noexcept
		: msg_("Input/Output error: " + error)
	{
	}

	const char* IOError::what() const noexcept
	{
		return msg_.c_str();
	}

	NotImplemented::NotImplemented(const char* file, int line, const std::string& method) noexcept
	{
		msg_ = "The method '" + method + "' is currently not implemented in file: "
		+ file  + ", line: ";

		try {
			msg_ += boost::lexical_cast<std::string>(line);
		} catch(boost::bad_lexical_cast& e) {
			msg_ = "Error creating exception. (Bad lexical cast)";
		}
	}

	const char* NotImplemented::what() const noexcept
	{
		return msg_.c_str();
	}

	InvalidIndex::InvalidIndex(unsigned int i, unsigned int max) noexcept
	{
		try {
			msg_ = std::string("Invalid index ") + boost::lexical_cast<std::string>(i)
			     + ". Strict upper bound is " + boost::lexical_cast<std::string>(max) + ".";
		} catch(boost::bad_lexical_cast& e) {
			msg_ = std::string("Error creating exception. (Bad lexical cast)");
		}
	}

	const char* InvalidIndex::what() const noexcept
	{
		return msg_.c_str();
	}

	InvalidKey::InvalidKey(const std::string& key) noexcept
	{
		try {
			msg_ = "Invalid key: '" + key + "'!";
		} catch(...) {
			// This method may not throw, however creating a string may
			// allocate memory and thus throw.
		}
	}

	const char* InvalidKey::what() const noexcept
	{
		return msg_.c_str();
	}

	UnknownEntry::UnknownEntry(const std::string& entry_as_string) noexcept
		: msg_("Error during lookup! Unknown entry \"" + entry_as_string + "\".")
	{
	}

	const char* UnknownEntry::what() const noexcept { return msg_.c_str(); }
}

