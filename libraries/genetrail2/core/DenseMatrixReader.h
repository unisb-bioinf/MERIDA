/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013-2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_DENSE_MATRIX_READER_H
#define GT2_DENSE_MATRIX_READER_H

#include "macros.h"

#include <istream>
#include <vector>
#include <set>
#include <initializer_list>

namespace GeneTrail
{
	class DenseMatrix;

	class GT2_EXPORT DenseMatrixReader
	{
		public:
			/**
			 * Options that modify the behaviour of the
			 * reader. The flags can be combined using bitwise or (|)
			 */
			enum ReaderOptions
			{
				NO_OPTIONS          = 0,      /// No special option should be set
				READ_ROW_NAMES      = 1 << 0, /// The file contains row names (only ascii)
				READ_COL_NAMES      = 1 << 1, /// The file contains column names (only ascii)
				TRANSPOSE           = 1 << 2, /// Read a transposed version of the matrix
				ADDITIONAL_COL_NAME = 1 << 3  /// The row names have a column name assigned to them (only ascii)
			};

			/**
			 * Returns a set of default options. These are:
			 * READ_COL_NAMES | READ_ROW_NAMES
			 */
			static unsigned int defaultOptions();

			/**
			 * Virtual destructor
			 */
			virtual ~DenseMatrixReader() {}

			/**
			 * Reads a DenseMatrix from the provided input stream.
			 *
			 * @param input a (seekable) stream of a matrix implementation
			 * @param opts a set of options that manipulate the behaviour of the reader
			 *
			 * @throws IOError an IOError is thrown if the provided stream is invalid or the
			 *                 matrix description is invalid.
			 */
			virtual DenseMatrix read(std::istream& input, unsigned int opts = defaultOptions()) const;

		private:

			/**
			 * This set contains symbols, that might occurr in text files, indication a NaN value.
			 */
			std::set<std::string> nan_like_symbols{"NA","NaN","NAN","nan","null","NULL"};

			DenseMatrix textRead_ (std::istream& input, unsigned int opts) const;

			enum ChunkType {
				HEADER   = 0x00,
				ROWNAMES = 0x01,
				COLNAMES = 0x02,
				DATA     = 0x03
			};

			/**
			 * If isBinary_ returns true, binaryRead_ will attempt to read a binary
			 * matrix from the specified file.
			 *
			 * The file format is based on chunks, that, with the exception of the
			 * header can be stored in an arbitrary order.
			 *
			 * The basic chunk format is as follows:
			 *
			 * CHUNK-ID  : uint8_t   --- The identifier of the chunk
			 * CHUNK-SIZE: uint64_t  --- The size of the CHUNK-DATA field in bytes
			 *                           forwarding the file cursor CHUNK-SIZE bytes
			 *                           will either reach the next CHUNK or EOF
			 * CHUNK-DATA: This field depends on the current chunk
			 *
			 * Chunk Types (ID):
			 *  * HEADER (0x00):
			 *   - ROW-COUNT:     uint32_t --- The number of rows of the matrix
			 *   - COL-COUNT:     uint32_t --- The number of columns of the matrix
			 *   - STORAGE-ORDER: uint8_t  --- Zero if the matrix is stored in row major
			 *                                 non-zero if it is stored in column major
			 *  * ROWNAMES (0x01):
			 *   Contains ROW-COUNT zero terminated names.
			 *
			 *  * COLNAMES (0x02):
			 *   Contains COL-COUNT zero terminated names.
			 *
			 *  * DATA (0x03):
			 *   Contains ROW-COUNT * COL-COUNT double values of entries
			 *   in the storage order specified in the header.
			 */
			DenseMatrix binaryRead_(std::istream& input, unsigned int opts = NO_OPTIONS) const;

			void skipEmptyLines_(std::istream& input, std::string& line) const;

			/**
			 * This method checks the magic number of a stream in order to decide
			 * whether the matrix is stored in binary or text format.
			 *
			 * The magic number searched for is the signed 12 byte sequence
			 * "BINARYMATRIX"
			 *
			 */
			bool isBinary_(std::istream& input) const;

			void readRowNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const;
			void readColNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const;
			void readData_    (std::istream& input, DenseMatrix& result, uint8_t storage_order) const;
			void readNames_   (std::istream& input, uint64_t chunk_size, std::vector<std::string>& names) const;
			void readChunkHeader_(std::istream& input, uint8_t& chunk_type, uint64_t& chunk_size) const;
			DenseMatrix readHeader_(std::istream& input, uint8_t& storage_order) const;
	};
}

#endif //GT2_DENSE_MATRIX_READER_H

