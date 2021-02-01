#include "Matrix.h"
#include <stdexcept>
#include <math.h>

template <typename T>
Matrix<T>::Matrix(){
}

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
std::string Matrix<T>::to_string() const{
    std::string res="";
    for(int i = 0; i < m_rows; i++){
        res+="[";
        for(int j = 0; j < m_cols; j++){
            res+=std::to_string(m_data[i * m_cols + j]);
            res+=", ";
        }
        res.erase(res.end()-2, res.end());
        res+="]\n";
    }
    return res;
}

template <typename T>
T Matrix<T>::elem_at(int i, int j) const{
    return m_data[i * m_cols + j];
}

template <typename T>
Matrix<T>* Matrix<T>::set_elem_at(int i, int j, T val){
    m_data[i * m_cols + j]=val;
    return this;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T>* res = new Matrix<T>(m_cols,m_rows);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(j,i,elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::extract_column(int j) const{
	Matrix<T> res(m_rows, 1);
	for (int i = 0; i < m_rows; i++)
		res.set_elem_at(i, 0, elem_at(i, j));
	return res;
}

template <typename T>
void Matrix<T>::fill_column(int j, const Matrix &col){
	for (int i = 0; i < m_rows; i++)
		set_elem_at(i, j, col.elem_at(i, 0));
}

template <typename T>
Matrix<T> Matrix<T>::dot(const Matrix &B) const{
    Matrix<T> res(m_rows,B.m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < B.m_cols; j++){
            T tot{};
            for(int k = 0; k < m_cols; k++){
                tot+=(elem_at(i,k)*B.elem_at(k,j));
            }
            res.set_elem_at(i,j,tot);
        }
    return res;
}

template <typename T>
void Matrix<T>::switch_row(int row1, int row2){
	for (int j = 0; j < m_cols; j++)
	{
		T temp = elem_at(row1, j);
		set_elem_at(row1, j, elem_at(row2, j));
		set_elem_at(row2, j, temp);
	}
}

// in order to set the pivot to 1
template <typename T>
void Matrix<T>::reduction(int row, T alpha){
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row, j, elem_at(row, j) * alpha);
	}
}

template <typename T>
void Matrix<T>::transvection(int row1, int row2, T alpha){
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row1, j, elem_at(row1, j) + elem_at(row2, j) * alpha);
	}
}

template <typename T>
Matrix<T> Matrix<T>::inverse() const{
    if (m_cols!=m_rows)
			throw std::domain_error("The matrix is not square");

	Matrix<T> Id(m_cols);
	Matrix<T> cop(m_rows, m_cols, m_data);

	int r = -1;
	for(int j = 0; j < m_cols; j++){
        int k=r+1;
        T cand_val = fabs(cop.elem_at(k,j));
        for(int i = r+2; i<m_rows; i++){
            if(fabs(cop.elem_at(i,j))>cand_val){
                k=i;
                cand_val=fabs(cop.elem_at(k,j));
            }
        }
		if (cand_val==T())
            throw std::invalid_argument("The matrix is not invertible");

        r++;
        Id.reduction(k, 1 / cop.elem_at(k, j));
		cop.reduction(k, 1 / cop.elem_at(k, j));
        if(k!=r){
            Id.switch_row(k, r);
            cop.switch_row(k, r);
        }

        for(int i = 0; i<m_rows; i++){
            if(i!=r){
                Id.transvection(i,r,-cop.elem_at(i,j));
                cop.transvection(i,r,-cop.elem_at(i,j));
            }
        }
	}
	return Id;
}

//Hadamard product
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &B) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,B.elem_at(i,j) * elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &scal) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,scal * elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &B) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,B.elem_at(i,j) + elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T &scal) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,scal + elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix &cop){
    m_rows=cop.m_rows;
    m_cols=cop.m_cols;
    //delete m_data ????????
    m_data = new T[m_rows*m_cols];
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  cop.m_data[i * m_cols + j];
    return *this;
}

template <typename T>
Matrix<T>::~Matrix(){
    delete[] m_data;
}

template class Matrix<int>;
template class Matrix<long>;
template class Matrix<float>;
template class Matrix<double>;
