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

#include "DenseMatrixReader.h"

#include <vector>
#include <deque>
#include <iostream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "DenseMatrix.h"
#include "Exception.h"

namespace GeneTrail
{
	unsigned int DenseMatrixReader::defaultOptions()
	{
		return READ_COL_NAMES | READ_ROW_NAMES;
	}

	void DenseMatrixReader::skipEmptyLines_(std::istream& input, std::string& line) const
	{
		while(std::getline(input, line))
		{
			// TODO this is dangerous if there is an empty row name
			boost::trim(line);

			if(!line.empty())
			{
				break;
			}
		}
	}

	bool DenseMatrixReader::isBinary_(std::istream& input) const
	{
		char magic[12];
		input.read(magic, 12);
		return strncmp(magic, "BINARYMATRIX", 12) == 0;
	}

	DenseMatrix DenseMatrixReader::read(std::istream& input, unsigned int opts) const
	{
		if(!input) {
			throw IOError("Invalid input stream!");
		}

		if(isBinary_(input)) {
			return binaryRead_(input);
		}

		input.seekg(0, std::ios::beg);
		return textRead_(input, opts);
	}

	void DenseMatrixReader::readChunkHeader_(std::istream& input, uint8_t& chunk_type, uint64_t& size) const
	{
		input.read((char*)&chunk_type, 1);
		input.read((char*)&size, 8);
	}

	DenseMatrix DenseMatrixReader::readHeader_(std::istream& input, uint8_t& storage_order) const
	{
		uint32_t row_count;
		uint32_t col_count;

		input.read((char*)&row_count, 4);
		input.read((char*)&col_count, 4);
		input.read((char*)&storage_order, 1);

		return DenseMatrix(row_count, col_count);
	}

	void DenseMatrixReader::readRowNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const
	{
		// Reserve some storage
		std::vector<std::string> names(result.rows());

		readNames_(input, chunk_size, names);

		result.setRowNames(names);
	}

	void DenseMatrixReader::readColNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const
	{
		// Reserve some storage
		std::vector<std::string> names(result.cols());

		readNames_(input, chunk_size, names);

		result.setColNames(names);
	}

	void DenseMatrixReader::readNames_(std::istream& input, uint64_t chunk_size, std::vector<std::string>& names) const
	{
		unsigned int i = 0;

		uint64_t bytes_read = 0;

		const int BUFFER_SIZE = 1024;

		char buffer[BUFFER_SIZE];
		while(input.good() && i < names.size() && bytes_read < chunk_size) {
			input.getline(buffer, BUFFER_SIZE, '\0');

			bytes_read += input.gcount(); // Keep track of the number of bytes read

			names[i] += buffer;

			if((input.rdstate() & std::istream::failbit) != 0) {
				input.clear();
			} else {
				++i;
			}
		}

		if(bytes_read != chunk_size) {
			std::cerr << "Parsing error" << std::endl;
		}
	}

	void DenseMatrixReader::readData_(std::istream& input, DenseMatrix& result, uint8_t storage_order) const
	{
		uint64_t bytes_read = 0;

		const uint64_t n = sizeof(DenseMatrix::value_type) * result.rows() * result.cols();

		if(storage_order == 0) {
			for(unsigned int i = 0; i < result.rows(); ++i) {
				for(unsigned int j = 0; j < result.cols(); ++j) {
					input.read((char*)&result(i,j), sizeof(DenseMatrix::value_type));
					bytes_read += input.gcount();
				}
			}
		} else {
			// As the internal storage format of matrix is column major this is quite efficient...
			input.read((char*)result.matrix().data(), n);
			bytes_read += input.gcount();
		}

		if(bytes_read != n) {
			std::cerr << "Parsing error" << std::endl;
		}
	}

	DenseMatrix DenseMatrixReader::binaryRead_(std::istream& input, unsigned int opts) const
	{
		uint8_t chunk_type = 0;
		uint64_t chunk_size = 0;

		// Read the header
		readChunkHeader_(input, chunk_type, chunk_size);

		if(chunk_type != 0x0)
		{
			throw new IOError("Unexpected chunk: expected 0 (matrix header), got " + boost::lexical_cast<std::string>(chunk_type));
		}

		if(chunk_size != 9) {
			throw IOError("Inconsistent header size: expected 9 got " + boost::lexical_cast<std::string>(chunk_size));
		}

		uint8_t  storage_order;
		DenseMatrix result = readHeader_(input, storage_order);

		while(input.good()) {
			readChunkHeader_(input, chunk_type, chunk_size);

			if(input.eof()) {
				break;
			}

			switch(chunk_type) {
				case DenseMatrixReader::HEADER:
					throw IOError("Unexpected chunk: did not expect header chunk!");
				case DenseMatrixReader::ROWNAMES:
					readRowNames_(input, result, chunk_size);
					break;
				case DenseMatrixReader::COLNAMES:
					readColNames_(input, result, chunk_size);
					break;
				case DenseMatrixReader::DATA:
					if(sizeof(DenseMatrix::value_type) * result.cols() * result.rows() != chunk_size) {
						std::cerr << "Inconsistent data chunk size!" << std::endl;
						return DenseMatrix(0,0);
					}

					readData_(input, result, storage_order);
					break;
				default:
					std::cout << "Unknown chunk " << chunk_type << " Skipping!" << std::endl;
					input.seekg(chunk_size, std::ios::cur);
			}
		}

		return result;
	}

	//TODO this is a bloody mess! Fix asap!
	DenseMatrix DenseMatrixReader::textRead_(std::istream& input, unsigned int opts) const
	{
		std::string line;
		// Get the first interesting line
		skipEmptyLines_(input, line);

		std::vector<std::string> fields;
		boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

		std::vector<std::string> row_names;
		std::vector<std::string> col_names;

		size_t colname_offset = ((opts & ADDITIONAL_COL_NAME) ? 1 : 0);

		if(opts & READ_COL_NAMES)
		{
			col_names.resize(fields.size() - colname_offset);

			for(unsigned int i = colname_offset; i < fields.size(); ++i)
			{
				col_names[i - colname_offset] = fields[i];
			}

			skipEmptyLines_(input, line);
			boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);
		}

		const size_t num_fields = fields.size();

		std::deque<double> data;

		const int start = (opts & READ_ROW_NAMES) ? 1 : 0;

		if(start && fields.size() > 0)
		{
			row_names.push_back(fields[0]);
		}
		for(size_t i = start; i < fields.size(); ++i)
		{
			try
			{
				data.push_back(boost::lexical_cast<double>(fields[i]));
			}
			catch(boost::bad_lexical_cast& e)
			{
				if(nan_like_symbols.find(fields[i]) == nan_like_symbols.end()) {
					throw IOError(e.what());
				} else {
					data.push_back(std::numeric_limits<double>::quiet_NaN());
				}
			}
		}

		unsigned int cur_line = 1;

		while(std::getline(input, line))
		{
			boost::trim(line);

			if(line.empty())
			{
				continue;
			}

			boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

			if(fields.size() != num_fields)
			{
				throw IOError(
					"Expected " + boost::lexical_cast<std::string>(num_fields) +
					" columns in line " + boost::lexical_cast<std::string>(cur_line) +
					", got " + boost::lexical_cast<std::string>(fields.size())
				);
			}

			if(start)
			{
				row_names.push_back(fields[0]);
			}

			for(size_t i = start; i < fields.size(); ++i)
			{
				try
				{
					data.push_back(boost::lexical_cast<double>(fields[i]));
				}
				catch(boost::bad_lexical_cast& e)
				{
					if(nan_like_symbols.find(fields[i]) == nan_like_symbols.end()) {
						throw IOError(e.what());
					} else {
						data.push_back(std::numeric_limits<double>::quiet_NaN());
					}
				}
			}

			++cur_line;
		}

		if(opts & TRANSPOSE)
		{
			DenseMatrix result(num_fields - start, data.size() / (num_fields - start));

			int i = 0;

			for(auto it = data.begin(); it != data.end(); ++it, ++i)
			{
				result(i % (num_fields - start), i / (num_fields - start)) = *it;
			}

			if(col_names.size() == result.rows()) result.setRowNames(col_names);
			if(row_names.size() == result.cols()) result.setColNames(row_names);

			return result;
		}
		else
		{
			DenseMatrix result(data.size() / (num_fields - start), num_fields - start);

			int i = 0;

			for(auto it = data.begin(); it != data.end(); ++it, ++i)
			{
				result(i / (num_fields - start), i % (num_fields - start)) = *it;
			}

			if(row_names.size() == result.rows()) result.setRowNames(row_names);
			if(col_names.size() == result.cols()) result.setColNames(col_names);

			return result;
		}
	}
}

