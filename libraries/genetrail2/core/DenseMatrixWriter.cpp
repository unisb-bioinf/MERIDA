/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013-2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include "DenseMatrixWriter.h"

#include "DenseMatrix.h"
#include "DenseMatrixSubset.h"

#include <iterator>

#include <iostream>

namespace GeneTrail
{
	uint64_t DenseMatrixWriter::writeBinary(std::ostream& output, const DenseMatrix& matrix) const
	{
		uint64_t
		total  = writeBinary_(output, matrix);
		total += writeData_(output, matrix);

		return total;
	}

	uint64_t DenseMatrixWriter::writeBinary(std::ostream& output, const DenseMatrixSubset& matrix) const
	{
		uint64_t
		total  = writeBinary_(output, matrix);
		total += writeData_(output, matrix);

		return total;
	}

	void DenseMatrixWriter::writeText(std::ostream& output, const DenseMatrix& matrix) const
	{
		writeText_(output, matrix);

		for(DenseMatrix::index_type i = 0; i < matrix.rows(); ++i) {
			output << "\n" << matrix.rowName(i);

			for(DenseMatrix::index_type j = 0; j < matrix.cols(); ++j) {
				output << "\t" << matrix(i,j);
			}
		}
	}

	void DenseMatrixWriter::writeText(std::ostream& output, const DenseMatrixSubset& matrix) const
	{
		writeText_(output, matrix);

		for(DenseMatrix::index_type i = 0; i < matrix.rows(); ++i) {
			output << "\n" << matrix.rowName(i);

			for(DenseMatrix::index_type j = 0; j < matrix.cols(); ++j) {
				output << "\t" << matrix(i,j);
			}
		}
	}

	uint64_t DenseMatrixWriter::writeData_(std::ostream& output, const DenseMatrix& matrix) const
	{
		const uint64_t n = matrix.rows() * matrix.cols() * sizeof(DenseMatrix::value_type);
		uint64_t total = writeChunkHeader_(output, 0x3, n);

		// As we are writing column major we can just
		// pipe the data to the output
		output.write((const char*)matrix.matrix().data(), n);

		return total + n;
	}

	uint64_t DenseMatrixWriter::writeData_(std::ostream& output, const DenseMatrixSubset& matrix) const
	{
		const uint64_t n = matrix.rows() * matrix.cols() * sizeof(Matrix::value_type);
		uint64_t total = writeChunkHeader_(output, 0x3, n);

		for(Matrix::index_type j = 0; j < matrix.cols(); ++j) {
			for(Matrix::index_type i = 0; i < matrix.rows(); ++i) {
				Matrix::value_type tmp = matrix(i,j);
				output.write((const char*)&tmp, sizeof(Matrix::value_type));
			}
		}

		return total + n;
	}
}

