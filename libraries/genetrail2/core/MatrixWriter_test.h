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

#ifndef GT2_MATRIX_WRITER_H
#define GT2_MATRIX_WRITER_H

#include "macros.h"

#include <ostream>
#include <vector>

namespace GeneTrail
{
	class Matrix;

	class GT2_EXPORT MatrixWriter
	{
		protected:
			void     writeText_       (std::ostream& output, const Matrix& matrix) const;
			uint64_t writeBinary_     (std::ostream& output, const Matrix& matrix) const;
			uint64_t writeChunkHeader_(std::ostream& output, uint8_t type, uint64_t size) const;
			uint64_t writeNames_      (std::ostream& output, const std::vector<std::string>& names) const;
			uint64_t writeHeader_     (std::ostream& output, const Matrix& matrix) const;
			uint64_t writeRowNames_   (std::ostream& output, const Matrix& matrix) const;
			uint64_t writeColNames_   (std::ostream& output, const Matrix& matrix) const;
	};
}

#endif //GT2_MATRIX_WRITER_H
