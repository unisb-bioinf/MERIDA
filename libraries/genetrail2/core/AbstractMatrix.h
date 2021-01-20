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
#ifndef GT2_ABSTRACT_MATRIX_H
#define GT2_ABSTRACT_MATRIX_H

#include "macros.h"

#include "Matrix.h"

#include <map>
#include <vector>
#include <string>

namespace GeneTrail
{
	class GT2_EXPORT AbstractMatrix : public Matrix
	{
		public:
			/// The precision used in the matrix
			typedef double       value_type;

			/// The index type used in the matrix
			typedef unsigned int index_type;

			///\ingroup Constructors
			///@{
			/**
			 * Default constructor
			 *
			 * Reserves the storage for a rows * cols matrix.
			 * The matrix will not be initialised.
			 */
			AbstractMatrix(index_type rows, index_type cols);

			/**
			 * Row-name constructor
			 *
			 * This constructs a rows.size() x cols.size() matrix and sets the row and column names
			 * to rows and cols respectively
			 */
			AbstractMatrix(std::vector<std::string>   rows, std::vector<std::string>   cols);

			/**
			 * Default copy constructor
			 */
			AbstractMatrix(const AbstractMatrix&) = default;

			/**
			 * Move Constructor
			 */
			AbstractMatrix(AbstractMatrix&& matrix);

			/**
			 * Virtual Destructor
			 */
			virtual ~AbstractMatrix() {};

			///@}
			///\ingroup Operators
			///@{

			/**
			 * Move assignment operator
			 */
			AbstractMatrix& operator=(AbstractMatrix&& matrix);

			///@}
			///\ingroup Accessors
			///@{

			/**
			 * Sets the row names of the matrix to the values specified in row_names
			 */
			void setRowNames(const std::vector<std::string>& row_names) override;

			/**
			 * Sets the column names of the matrix to the values specified in col_names
			 */
			void setColNames(const std::vector<std::string>& col_names) override;

			/**
			 * Return the column names
			 */
			const std::vector<std::string>& colNames() const override;

			/**
			 * Return the row names
			 */
			const std::vector<std::string>& rowNames() const override;

			/**
			 * Rename row old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			void setRowName(const std::string& old_name, const std::string& new_name) override;

			/**
			 * Rename row i to new_name
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			void setRowName(index_type i, const std::string& new_name) override;

			/**
			 * Get the name of row i
			 */
			const std::string& rowName(index_type i) const override;

			/**
			 * Rename column old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			void setColName(const std::string& old_name, const std::string& new_name) override;

			/**
			 * Rename column j to new_name
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			void setColName(index_type j, const std::string& new_name) override;

			/**
			 * Get the name of column j
			 */
			const std::string& colName(index_type j) const override;

			/**
			 * Return the row index of the column identified by
			 * row.
			 *
			 * \warning the result of this method is undefined if
			 *          row is not a valid row name
			 */
			index_type rowIndex(const std::string& row) const override;

			/**
			 * Return the column index of the column identified by
			 * column.
			 *
			 * \warning the result of this method is undefined if
			 *          column is not a valid column name
			 */
			index_type colIndex(const std::string& col) const override;

			/**
			 * Return true if the matrix contains a row with row name
			 * "name". False otherwise.
			 */
			bool hasRow(const std::string& name) const override;

			/**
			 * Return true if the matrix contains a column with column name
			 * "name". False otherwise.
			 */
			bool hasCol(const std::string& name) const override;

			/**
			 * The number of columns in the matrix
			 */
			index_type cols() const override;

			/**
			 * The number of rows in the matrix
			 */
			index_type rows() const override;

			///@}
			///\ingroup Matrix Operations
			///@{

			/**
			 * \copydoc Matrix::transpose
			 */
			virtual void transpose() override;

			///@}

		protected:
			// Containers for mapping row names
			std::vector<std::string> index_to_rowname_;
			std::map<std::string, index_type> rowname_to_index_;

			// Containers for mapping column names
			std::vector<std::string> index_to_colname_;
			std::map<std::string, index_type> colname_to_index_;

			void updateRowAndColNames_();

			// Helper for renaming rows or columns
			void setName_(index_type j,
			              const std::string& new_name,
			              std::map<std::string, index_type>& name_to_index,
			              std::vector<std::string>& index_to_name);

			void setName_(const std::string& old_name,
			              const std::string& new_name,
			              std::map<std::string, index_type>& name_to_index,
			              std::vector<std::string>& index_to_name);
	};
}

#endif //GT2_MATRIX_H

