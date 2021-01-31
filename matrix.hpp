#ifndef MATRIX_H
#define MATRIX_H

#include <string>

typedef unsigned long ul;

template <typename T=double>
class Matrix {
    public:
        Matrix();
        explicit Matrix(ul N);
        explicit Matrix(ul M, ul N);
        explicit Matrix(ul M, ul N, T* dat);

        void display() const;

        T elem_at(ul i, ul j) const;
        Matrix* set_elem_at(ul i, ul j, T val);

		Matrix inverse();
		void switch_row(ul row1, ul row2);
		void reduction(ul row, double alpha);
		void transvection(ul row1, ul row2, double alpha);  
		void zeros(ul row, ul col);
		void inverse_row(Matrix<T> &M, ul row, ul col);
		bool is_still_inversible(ul row, ul col);

        Matrix transpose() const;
        Matrix dot(const Matrix &B) const;
		Matrix column(ul j);
		void fill_column(ul j, const Matrix &data);
        Matrix operator*(const Matrix &B) const;
        Matrix operator*(const T &scal) const;
        Matrix operator+(const Matrix &B) const;
        Matrix operator+(const T &scal) const;

        virtual ~Matrix();

    private:
        ul m_rows;
        ul m_cols;
        T* m_data;

};

#endif
