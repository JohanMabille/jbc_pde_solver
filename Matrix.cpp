#include "Matrix.h"
#include <iostream>

template <typename T>
Matrix<T>::Matrix(){
    //ctor
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
Matrix<T> Matrix<T>::dot(const Matrix &B) const{
    Matrix<T>* res = new Matrix<T>(m_rows,B.m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < B.m_cols; j++){
            T tot{};
            for(int k = 0; k < m_cols; k++)
                tot+=(elem_at(i,k)*B.elem_at(k,j));
            res->set_elem_at(i,j,tot);
        }
    return *res;
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
Matrix<T>::~Matrix(){
    delete[] m_data;
}

template class Matrix<int>;
template class Matrix<long>;
template class Matrix<float>;
template class Matrix<double>;
