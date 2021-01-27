#include "matrix.hpp"
#include <iostream>
#include <stdexcept>
#include <stdio.h>

template <typename T>
Matrix<T>::Matrix() {}

template <typename T>
Matrix<T>::Matrix(int M, int N): m_rows(M), m_cols(N), m_data(new T[M * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  T();
}

template <typename T>
Matrix<T>::Matrix(int N): m_rows(N), m_cols(N), m_data(new T[N * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  T(i==j);

}

template <typename T>
Matrix<T>::Matrix(int M, int N, T* dat): m_rows(M), m_cols(N), m_data(new T[M * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  dat[i * m_cols + j];
}

template <typename T>
void Matrix<T>::display() const {
    std::string res="";
    for(int i = 0; i < m_rows; i++) {
        res+="[";
        for(int j = 0; j < m_cols; j++) {
            res+=std::to_string(m_data[i * m_cols + j]);
            res+=", ";
        }
        res.erase(res.end()-2, res.end());
        res+="]\n";
    }
    std::cout << res << std::endl;
}

template <typename T>
T Matrix<T>::elem_at(int i, int j) const {
	if (i >= m_rows) {
		throw std::invalid_argument("Either the index of row is larger than the number of rows");
	} else if (j >= m_cols) {
		throw std::invalid_argument("Either the index of column is larger than the number of columns");
	}
    return m_data[i * m_cols + j];
}

template <typename T>
Matrix<T>* Matrix<T>::set_elem_at(int i, int j, T val) {
    m_data[i * m_cols + j] = val;
    return this;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> res = Matrix<T>(m_cols, m_rows);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++) {
            res.set_elem_at(j, i, elem_at(i, j));
        }
    return res;
}

template <typename T>
void Matrix<T>::switch_row(int row1, int row2) {
	for (int j = 0; j < m_cols; j++)
	{
		double temp = elem_at(row1, j);
		set_elem_at(row1, j, elem_at(row2, j));
		set_elem_at(row2, j, temp);
	}
}


template <typename T>
// in order the pivot is equal to 1 
void Matrix<T>::reduction(int row, double alpha) {
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row, j, elem_at(row, j) * alpha);
	}
}

template <typename T>
void Matrix<T>::transvection(int row1, int row2, double alpha) {
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row1, j, elem_at(row1, j) + elem_at(row2, j) * alpha);
	}
}

template <typename T>
// We create 0 abose and below the pivot
void Matrix<T>::zeros(int row, int col) {
	for (int k = 0; k < m_rows; k++) {
		double c = elem_at(k, col);
		if (k != row) {
			for (int l = 0; l < m_cols; l++)
				set_elem_at(k, l, elem_at(k, l) - c * elem_at(row, l));
		}
	}
}

template <typename T>
void Matrix<T>::inverse_row(Matrix<T> &M, int row, int col) {
	for (int k = 0; k < m_rows; k++) {
		double c = elem_at(k, col);
		if (k != row) {
			for (int l = 0; l < m_cols; l++)
				M.set_elem_at(k, l, M.elem_at(k, l) - c * M.elem_at(row, l));
		}
	}
}

template <typename T>
// we make sure the column is not full of 0, otherwise, it means that the matrix is not inversible
bool Matrix<T>::is_still_inversible(int row, int col) {
	for (int l = row; l < m_rows; l++) {
		if (elem_at(l, col) != 0)
			return true;
	}
	return false;
}

template <typename T>
Matrix<T> Matrix<T>::inverse() {
	Matrix<T> iden = Matrix<T>(m_cols);
	Matrix<T> rhs = Matrix<T>(m_rows, m_cols, m_data);

	int l = 0;
	int i = 0;
	for (int j = 0; j < m_cols; j++) {
		if (!rhs.is_still_inversible(i, j)) {
			throw std::invalid_argument("The coef matrix is not inversible");
		}

		l = i;
		while (rhs.elem_at(l, j) == 0) {
			l++;
		}

		iden.switch_row(l, i);
		rhs.switch_row(l, i);

		iden.reduction(i, 1 / elem_at(i, j));
		rhs.reduction(i, 1 / elem_at(i, j));

		rhs.inverse_row(iden, i, j);

		rhs.zeros(i, j);
		i++;
	}
	return iden;
}


template <typename T>
Matrix<T> Matrix<T>::dot(const Matrix &B) const {
    Matrix<T> res = Matrix<T>(m_rows,B.m_cols);
	for (int i = 0; i < m_rows; i++) {
		for (int j = 0; j < B.m_cols; j++) {
			T tot = 0.;
			for (int k = 0; k < m_cols; k++)
				tot += (elem_at(i, k) * B.elem_at(k, j));
			res.set_elem_at(i, j, tot);
		}
	}
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::column(int j) {
	Matrix<T> res = Matrix<T>(m_rows, 1);
	for (int i = 0; i < m_rows; i++)
		res.set_elem_at(i, 0, elem_at(i, j));
	return res;
}

template <typename T>
void Matrix<T>::fill_column(int j, const Matrix<T> &data) {
	for (int i = 0; i < m_rows; i++)
		set_elem_at(i, j, data.elem_at(i, 0));
}

//Hadamard product
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &B) const {
    Matrix<T> res = Matrix<T>(m_rows, m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            res.set_elem_at(i, j, B.elem_at(i,j) * elem_at(i, j));
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &scal) const{
    Matrix<T> res = Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            res.set_elem_at(i,j,scal * elem_at(i,j));
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &B) const {
    Matrix<T> res = Matrix<T>(m_rows, m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            res.set_elem_at(i, j, B.elem_at(i, j) + elem_at(i, j));
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T &scal) const {
    Matrix<T> res = Matrix<T>(m_rows, m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            res.set_elem_at(i, j, scal + elem_at(i, j));
    return res;
}


template <typename T>
Matrix<T>::~Matrix() {
    delete[] m_data;
}

template class Matrix<int>;
template class Matrix<long>;
template class Matrix<float>;
template class Matrix<double>;
