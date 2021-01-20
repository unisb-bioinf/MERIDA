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

#ifndef GT2_MATRIX_H
#define GT2_MATRIX_H

#include "macros.h"

#include <map>
#include <vector>
#include <string>

#include "Exception.h"

namespace GeneTrail
{
	class GT2_EXPORT Matrix
	{
		public:
			/// The precision used in the matrix
			typedef double       value_type;

			/// The index type used in the matrix
			typedef unsigned int index_type;

			/**
			 * Virtual Destructor
			 */
			virtual ~Matrix() {};

			///\ingroup Accessors
			///@{

			/**
			 * Sets the row names of the matrix to the values specified in row_names
			 */
			virtual void setRowNames(const std::vector<std::string>& row_names) = 0;

			/**
			 * Sets the column names of the matrix to the values specified in col_names
			 */
			virtual void setColNames(const std::vector<std::string>& col_names) = 0;

			/**
			 * Return the column names
			 */
			virtual const std::vector<std::string>& colNames() const = 0;

			/**
			 * Return the row names
			 */
			virtual const std::vector<std::string>& rowNames() const = 0;

			/**
			 * Rename row old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			virtual void setRowName(const std::string& old_name, const std::string& new_name) = 0;

			/**
			 * Rename row i to new_name
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			virtual void setRowName(index_type i, const std::string& new_name) = 0;

			/**
			 * Get the name of row i
			 */
			virtual const std::string& rowName(index_type i) const = 0;

			/**
			 * Rename column old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			virtual void setColName(const std::string& old_name, const std::string& new_name) = 0;

			/**
			 * Rename column j to new_name
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			virtual void setColName(index_type j, const std::string& new_name) = 0;

			/**
			 * Get the name of column j
			 */
			virtual const std::string& colName(index_type j) const = 0;

			/**
			 * Return the row index of the column identified by
			 * row.
			 *
			 * \warning the result of this method is undefined if
			 *          row is not a valid row name
			 */
			virtual index_type rowIndex(const std::string& row) const = 0;

			/**
			 * Return the column index of the column identified by
			 * column.
			 *
			 * \warning the result of this method is undefined if
			 *          column is not a valid column name
			 */
			virtual index_type colIndex(const std::string& col) const = 0;

			/**
			 * Return true if the matrix contains a row with row name
			 * "name". False otherwise.
			 */
			virtual bool hasRow(const std::string& name) const = 0;

			/**
			 * Return true if the matrix contains a column with column name
			 * "name". False otherwise.
			 */
			virtual bool hasCol(const std::string& name) const = 0;

			/**
			 * The number of columns in the matrix
			 */
			virtual index_type cols() const = 0;

			/**
			 * The number of rows in the matrix
			 */
			virtual index_type rows() const = 0;

			virtual value_type& operator()(index_type i, index_type j) = 0;
			virtual value_type  operator()(index_type i, index_type j) const = 0;

			value_type set(index_type i, index_type j, value_type v) {
				if(i >= rows()) {
					throw InvalidIndex(i, rows());
				}

				if(j >= cols()) {
					throw InvalidIndex(i, cols());
				}

				return (*this)(i,j) = v;
			}

			value_type  get(index_type i, index_type j) const {
				if(i >= rows()) {
					throw InvalidIndex(i, rows());
				}

				if(j >= cols()) {
					throw InvalidIndex(i, cols());
				}

				return (*this)(i,j);
			}

			///@}
			///\ingroup Matrix Operations
			///@{

			/**
			 * Reorder the rows of the matrix according to a given
			 * permutation.
			 *
			 * \param perm A vector of matrix indices. The new matrix
			 *             will consist of the rows at the indices
			 *             contained in perm in the order they were
			 *             specified.
			 */
			virtual void shuffleRows(const std::vector<index_type>& perm) = 0;

			/**
			 * Reorder the columns of the matrix according to a given
			 * permutation.
			 *
			 * \see Matrix::shuffleRows
			 */
			virtual void shuffleCols(const std::vector<index_type>& perm) = 0;

			/**
			 * Remove the rows identified by the passed indices
			 *
			 * \param indices a vector of row indices
			 * \warning the passed indices must be sorted in ascending order
			 */
			virtual void removeRows(const std::vector<index_type>& indices) = 0;

			/**
			 * Remove the columns identified by the passed indices
			 *
			 * \see Matrix::removeRows
			 */
			virtual void removeCols(const std::vector<index_type>& indices) = 0;

			/**
			 * Transpose the matrix.
			 *
			 * Row names and col names will also be switched.
			 */
			virtual void transpose() = 0;

			///@}
	};
}

#endif //GT2_MATRIX_H

