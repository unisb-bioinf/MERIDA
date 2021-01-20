#include "matrixTools.h"

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>

#include <iostream>
#include <fstream>

namespace GeneTrail
{

	std::vector<unsigned int>
	getIndices(const DenseMatrix& matrix,
	           const std::vector<std::string>& colnames,
	           const std::string& groupname)
	{
		std::vector<unsigned int> indices;
		indices.reserve(colnames.size());
		for(const auto& s : colnames) {
			if(matrix.hasCol(s)) {
				indices.emplace_back(matrix.colIndex(s));
			} else {
				std::cerr << "WARNING: Could not find column \"" + s + "\".\n";
			}
		}

		if(indices.empty()) {
			throw EmptyGroup(groupname);
		}

		return indices;
	}

	std::tuple<DenseMatrixSubset, DenseMatrixSubset>
	splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference,
	            const std::vector<std::string>& test)
	{
		return std::make_tuple(
		    DenseMatrixSubset::createColSubset(
		        &matrix, getIndices(matrix, reference, "reference")),
		    DenseMatrixSubset::createColSubset(
		        &matrix, getIndices(matrix, test, "test")));
	}

	std::tuple<DenseMatrixSubset, DenseMatrixSubset>
	splitMatrixRows(DenseMatrix& matrix, const std::vector<std::string> & reference, const std::vector<std::string> & test)
	{
		return std::make_tuple(
				DenseMatrixSubset::createRowSubset(&matrix, reference),
				DenseMatrixSubset::createRowSubset(&matrix, test)
		);
	}

	DenseMatrix buildDenseMatrix(const std::string& expr1,
	                             const std::string& expr2,
	                             const MatrixReaderOptions& options)
	{
		auto m1 = readDenseMatrix(expr1, options);
		if(expr2 != "") {
			auto m2 = readDenseMatrix(expr2, options);
			m1.cbind(m2);
		}
		return m1;
	}

	DenseMatrix readDenseMatrix(const std::string& matrix,
	                            const MatrixReaderOptions& options)
	{
		unsigned int opts = DenseMatrixReader::NO_OPTIONS;

		if(!options.no_rownames) {
			opts |= DenseMatrixReader::READ_ROW_NAMES;
		}

		if(!options.no_colnames) {
			opts |= DenseMatrixReader::READ_COL_NAMES;
		}

		if(options.additional_colname) {
			opts |= DenseMatrixReader::ADDITIONAL_COL_NAME;
		}

		DenseMatrixReader reader;

		std::ifstream strm(matrix, std::ios::binary);
		return reader.read(strm, opts);
	}
}
