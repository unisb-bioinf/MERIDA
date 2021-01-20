/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include "macros.h"

#include "AbstractMatrix.h"

#include <Eigen/Core>

#include <vector>
#include <map>

namespace GeneTrail
{
	/**
	 * A wrapper around MatrixXd that attaches row and column names.
	 *
	 * Additionally some convenience matrix manipulation and reshaping
	 * functions have been added. Serialize and deserialize this class
	 * using DenseMatrixWriter/Reader.
	 *
	 * \note DenseMatrix implements move constructor and assignment operators
	 * this means it is efficient to return temporary objects from a
	 * function.
	 */
	class GT2_EXPORT DenseMatrix : public AbstractMatrix
	{
		public:
			/// The Eigen class used for representing the internal matrix
			typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> DMatrix;

			/// The Eigen class used for representing rows and columns
			typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> Vector;

			/**
			 * Default constructor
			 *
			 * Reserves the storage for a rows * cols matrix.
			 * The matrix will not be initialised.
			 */
			DenseMatrix(index_type rows, index_type cols);

			/**
			 * Row-name constructor
			 *
			 * This constructs a rows.size() x cols.size() matrix and sets the row and column names
			 * to rows and cols respectively
			 */
			DenseMatrix(std::vector<std::string> rows, std::vector<std::string> cols);

			/**
			 * Default copy constructor
			 */
			DenseMatrix(const DenseMatrix&) = default;

			/**
			 * Move Constructor
			 */
			DenseMatrix(DenseMatrix&& matrix);

			/**
			 * Return the i-th row.
			 */
			DMatrix::RowXpr row(index_type i);

			/**
			 * Return the i-th row. Const version.
			 */
			DMatrix::ConstRowXpr row(index_type i) const;

			/**
			 * Set the row identified by "name" to v
			 *
			 * This method is a noop if "name" is not present in the matrix
			 */
			void setRow(const std::string& name, const Vector& v);

			/**
			 * Set the i-th row to v
			 */
			void setRow(index_type i, const Vector& v);

			/**
			 *
			 */
			void rbind(const DenseMatrix& m);

			/**
			 *
			 */
			void cbind(const DenseMatrix& m);

			/**
			 * Return the j-th column
			 */
			DMatrix::ColXpr col(index_type j);

			/**
			 * Return the j-th column
			 */
			DMatrix::ConstColXpr col(index_type j) const;

			/**
			 * Set the column identified by "name" to v
			 *
			 * This method is a noop if "name" is not present in the matrix
			 */
			void setCol(const std::string& name, const Vector& v);

			/**
			 * Set the i-th column to v
			 */
			void setCol(index_type j, const Vector& v);

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			DMatrix& matrix();

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			const DMatrix& matrix() const;

			///@}
			///\ingroup Operators
			///@{

			/**
			 * Returns a reference to the matrix coefficient at position (i,j)
			 */
			value_type& operator()(index_type i, index_type j) override
			{
				return m_(i, j);
			}

			/**
			 * Returns the matrix coefficient at position (i,j)
			 */
			value_type operator()(index_type i, index_type j) const override
			{
				return m_(i, j);
			}

			/**
			 * Move assignment operator
			 */
			DenseMatrix& operator=(DenseMatrix&& matrix);
			///@}

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
			void shuffleRows(const std::vector<index_type>& perm) override;

			/**
			 * Reorder the columns of the matrix according to a given
			 * permutation.
			 *
			 * \see DenseMatrix::shuffleRows
			 */
			void shuffleCols(const std::vector<index_type>& perm) override;

			/**
			 * Remove the rows identified by the passed indices
			 *
			 * \param indices a vector of row indices
			 * \warning the passed indices must be sorted in ascending order
			 */
			void removeRows(const std::vector<index_type>& indices) override;

			/**
			 * Remove the columns identified by the passed indices
			 *
			 * \see DenseMatrix::removeRows
			 */
			void removeCols(const std::vector<index_type>& indices) override;

			/**
			 * \copydoc AbstractMatrix::transpose
			 *
			 * \note If you just need the transpose of the matrix in a computation use
			 *       .matrix().transpose()
			 */
			void transpose() override;

			///@}
		private:
			// Actual matrix payload
			DMatrix m_;


			// Helper for removing rows or columns
			void remove_(const std::vector<index_type>& indices,
			                   std::map<std::string, index_type>& name_to_index,
			                   std::vector<std::string>& index_to_name,
			                   std::function<void(index_type, index_type)> copy);

			void shuffle_(std::vector<index_type> perm,
			              std::map<std::string, index_type>& name_to_index,
			              std::vector<std::string>& index_to_name,
			              std::function<void(index_type, index_type)> swap);
	};
}

#endif // DENSEMATRIX_H

