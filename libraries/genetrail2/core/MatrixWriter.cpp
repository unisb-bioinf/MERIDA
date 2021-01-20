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

#include "MatrixWriter.h"
#include "Matrix.h"

#include <iterator>
#include <iostream>

namespace GeneTrail
{
	uint64_t MatrixWriter::writeBinary_(std::ostream& output, const Matrix& matrix) const
	{
		output.write("BINARYMATRIX", 12);

		// Write Header Chunk
		uint64_t total = 12;
		total += writeHeader_(output, matrix);
		total += writeRowNames_(output, matrix);
		total += writeColNames_(output, matrix);

		return total;
	}

	void MatrixWriter::writeText_(std::ostream& output, const Matrix& matrix) const
	{
		if(matrix.cols() == 0) {
			return;
		}

		output << matrix.colName(0);
		if(matrix.colName(0) == "") {
			std::cerr << "Warning: empty column name supplied in column 0.\n";
		}

		for(Matrix::index_type j = 1; j < matrix.cols(); ++j) {
			if(matrix.colName(j) == "") {
				std::cerr << "Warning: empty column name supplied in column " << j << ".\n";
			}

			output << "\t" << matrix.colName(j);
		}
	}

	uint64_t MatrixWriter::writeChunkHeader_(std::ostream& output, uint8_t type, uint64_t size) const
	{
		output.write((char*)&type, 1);
		output.write((char*)&size, 8);

		return 9;
	}

	uint64_t MatrixWriter::writeNames_(std::ostream& output, const std::vector<std::string>& names) const
	{
		uint64_t bytes_written = 0;

		for(const auto& s : names) {
			output.write(s.c_str(), s.length() + 1);
			bytes_written += s.length() + 1;
		}

		// Fill in the chunk size
		output.seekp(-bytes_written - 8, std::ios::cur);
		output.write((char*)&bytes_written, 8);
		output.seekp(bytes_written, std::ios::cur);

		return bytes_written;
	}

	uint64_t MatrixWriter::writeHeader_(std::ostream& output, const Matrix& matrix) const
	{
		uint64_t total = 0;

		total += writeChunkHeader_(output, 0x0, 0x9);
		uint32_t row_count = matrix.rows();
		uint32_t col_count = matrix.cols();
		uint8_t storage_order = 0x1;

		output.write((char*)&row_count, 4);
		output.write((char*)&col_count, 4);
		output.write((char*)&storage_order, 1);

		return total + 9;
	}

	uint64_t MatrixWriter::writeRowNames_(std::ostream& output, const Matrix& matrix) const
	{
		uint64_t total = 0;
		total += writeChunkHeader_(output, 0x1, 0x0);
		total += writeNames_(output, matrix.rowNames());

		return total;
	}

	uint64_t MatrixWriter::writeColNames_(std::ostream& output, const Matrix& matrix) const
	{
		uint64_t total = 0;

		// Write the chunk header. We temporarily put in 0x0 for the size
		// this will be rectified by writeNames_
		total += writeChunkHeader_(output, 0x2, 0x0);
		total += writeNames_(output, matrix.colNames());

		return total;
	}
}
