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

#include "DenseMatrixSubset.h"

#include "Exception.h"

namespace GeneTrail
{
	/*
	 * Static members
	 */
	DenseMatrixSubset DenseMatrixSubset::createColSubset(DenseMatrix* mat, ISubset cols)
	{
		ISubset rows(mat->rows());
		for(size_t i = 0; i < rows.size(); ++i) {
			rows[i] = i;
		}

		return DenseMatrixSubset(mat, std::move(rows), std::move(cols));
	}

	DenseMatrixSubset DenseMatrixSubset::createColSubset(DenseMatrix* mat, const SSubset& cols)
	{
		ISubset subs(cols.size());

		int i = 0;
		for(const auto& s : cols) {
			subs[i] = mat->colIndex(s);
			++i;
		}

		return createColSubset(mat, subs);
	}

	DenseMatrixSubset DenseMatrixSubset::createRowSubset(DenseMatrix* mat, ISubset rows)
	{
		ISubset cols(mat->cols());
		for(size_t i = 0; i < cols.size(); ++i) {
			cols[i] = i;
		}

		return DenseMatrixSubset(mat, std::move(rows), std::move(cols));
	}

	DenseMatrixSubset DenseMatrixSubset::createRowSubset(DenseMatrix* mat, const SSubset& rows)
	{
		ISubset subs(rows.size());

		int i = 0;
		for(const auto& s : rows) {
			subs[i] = mat->rowIndex(s);
			++i;
		}

		return createRowSubset(mat, subs);
	}

	DenseMatrixSubset DenseMatrixSubset::createSubset(DenseMatrix* mat, const SSubset& rows, const SSubset& cols)
	{
		ISubset rsubs(rows.size());
		ISubset csubs(cols.size());

		int i = 0;
		for(const auto& s : rows) {
			rsubs[i] = mat->rowIndex(s);
			++i;
		}

		i = 0;
		for(const auto& s : cols) {
			csubs[i] = mat->colIndex(s);
			++i;
		}

		return DenseMatrixSubset(mat, rsubs, csubs);
	}

	/*
	 * Constructors
	 */
	DenseMatrixSubset::DenseMatrixSubset(DenseMatrix* mat, ISubset rows, ISubset cols)
		: mat_(mat),
		  row_subset_(std::move(rows)),
		  col_subset_(std::move(cols))
	{
	}

	DenseMatrixSubset::DenseMatrixSubset(DenseMatrixSubset&& subs)
		: mat_(subs.mat_),
		  row_subset_(std::move(subs.row_subset_)),
		  col_subset_(std::move(subs.col_subset_))
	{
	}

	/*
	 * Operators
	 */
	DenseMatrixSubset& DenseMatrixSubset::operator=(DenseMatrixSubset&& subs)
	{
		assert(this != &subs);

		mat_ = subs.mat_;
		row_subset_ = std::move(subs.row_subset_);
		col_subset_ = std::move(subs.col_subset_);

		return *this;
	}

	Matrix::index_type DenseMatrixSubset::cols() const
	{
		return col_subset_.size();
	}

	Matrix::index_type DenseMatrixSubset::rows() const
	{
		return row_subset_.size();
	}

	bool DenseMatrixSubset::hasRow(const std::string& name) const
	{
		for(auto i : row_subset_) {
			if(mat_->rowName(i) == name) {
				return true;
			}
		}

		return false;
	}

	bool DenseMatrixSubset::hasCol(const std::string& name) const
	{
		for(auto i : col_subset_) {
			if(mat_->colName(i) == name) {
				return true;
			}
		}

		return false;
	}

	Matrix::index_type DenseMatrixSubset::rowIndex(const std::string& row) const
	{
		int idx = 0;

		for(auto i : row_subset_) {

			if(mat_->rowName(i) == row) {
				return idx;
			}
			++idx;
		}

		return -1;
	}

	const std::string& DenseMatrixSubset::rowName(Matrix::index_type i) const
	{
		return mat_->rowName(row_subset_[i]);
	}

	const std::vector< std::string >& DenseMatrixSubset::rowNames() const
	{
		row_names_cache_.resize(row_subset_.size());

		for(size_t i = 0; i < row_subset_.size(); ++i) {
			row_names_cache_[i] = mat_->rowName(row_subset_[i]);
		}

		return row_names_cache_;
	}

	Matrix::index_type DenseMatrixSubset::colIndex(const std::string& col) const
	{
		int idx = 0;

		for(auto i : col_subset_) {

			if(mat_->colName(i) == col) {
				return idx;
			}
			++idx;
		}

		return -1;
	}

	const std::string& DenseMatrixSubset::colName(Matrix::index_type j) const
	{
		return mat_->colName(col_subset_[j]);
	}

	const std::vector< std::string >& DenseMatrixSubset::colNames() const
	{
		col_names_cache_.resize(col_subset_.size());

		for(size_t i = 0; i < col_subset_.size(); ++i) {
			col_names_cache_[i] = mat_->colName(col_subset_[i]);
		}

		return col_names_cache_;
	}

	void DenseMatrixSubset::setColName(Matrix::index_type j, const std::string& new_name)
	{
		mat_->setColName(col_subset_[j], new_name);
	}

	void DenseMatrixSubset::setColName(const std::string& old_name, const std::string& new_name)
	{
		setColName(colIndex(old_name), new_name);
	}

	void DenseMatrixSubset::setColNames(const std::vector< std::string >& col_names)
	{
		int i = 0;
		for(const auto& s : col_names) {
			mat_->setColName(col_subset_[i], s);
			++i;
		}
	}

	void DenseMatrixSubset::setRowName(Matrix::index_type i, const std::string& new_name)
	{
		mat_->setRowName(row_subset_[i], new_name);
	}

	void DenseMatrixSubset::setRowName(const std::string& old_name, const std::string& new_name)
	{
		setRowName(rowIndex(old_name), new_name);
	}

	void DenseMatrixSubset::setRowNames(const std::vector< std::string >& row_names)
	{
		int i = 0;
		for(const auto& s : row_names) {
			mat_->setRowName(row_subset_[i], s);
			++i;
		}
	}

	/*
	 * Matrix Operations
	 */
	void DenseMatrixSubset::remove_(const std::vector< Matrix::index_type >& indices, ISubset& subset)
	{
		size_t next_idx = 1;
		size_t write_idx = indices[0];

		for(size_t read_idx = indices[0] + 1; read_idx < subset.size(); ++read_idx) {
			if(next_idx < indices.size() && read_idx == indices[next_idx]) {
				assert(indices[next_idx -1] < indices[next_idx]);
				++next_idx;
			}
			else
			{
				subset[write_idx] = subset[read_idx];
				++write_idx;
			}
		}

		subset.resize(subset.size() - indices.size());
	}

	void DenseMatrixSubset::removeCols(const std::vector< Matrix::index_type >& indices)
	{
		remove_(indices, col_subset_);
	}

	void DenseMatrixSubset::removeRows(const std::vector< Matrix::index_type >& indices)
	{
		remove_(indices, row_subset_);
	}

	void DenseMatrixSubset::shuffleCols(const std::vector< Matrix::index_type >& perm)
	{
		ISubset tmp(col_subset_.size());

		for(size_t i = 0; i < perm.size(); ++i) {
			tmp[i] = col_subset_[perm[i]];
		}

		std::swap(col_subset_, tmp);
	}

	void DenseMatrixSubset::shuffleRows(const std::vector< Matrix::index_type >& perm)
	{
		ISubset tmp(row_subset_.size());

		for(size_t i = 0; i < perm.size(); ++i) {
			tmp[i] = row_subset_[perm[i]];
		}

		std::swap(row_subset_, tmp);
	}

	void DenseMatrixSubset::transpose()
	{
		throw NotImplemented(__FILE__, __LINE__, "DenseMatrixSubset::transpose()");
	}

}

