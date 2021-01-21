#ifndef MATRIX_H
#define MATRIX_H
#include <string>


template <typename T=double>
class Matrix
{
    public:
        Matrix();
        explicit Matrix(int N);
        explicit Matrix(int M, int N);
        explicit Matrix(int M, int N, T* dat);

        std::string to_string() const;

        T elem_at(int i, int j) const;
        Matrix* set_elem_at(int i, int j, T val);

        Matrix transpose() const;
        Matrix dot(const Matrix &B) const;
        Matrix operator*(const Matrix &B) const;
        Matrix operator*(const T &scal) const;
        Matrix operator+(const Matrix &B) const;
        Matrix operator+(const T &scal) const;

        virtual ~Matrix();

    protected:

    private:
        int m_rows;
        int m_cols;
        T* m_data;

};

#endif // MATRIX_H
