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
        // Why returning a pointer to the object?
        // You could return a reference if you want to use
        // chained methods, but that's unlikely for classes
        // like this
        Matrix* set_elem_at(int i, int j, T val);

        Matrix extract_column(int j) const;
        void fill_column(int j, const Matrix &data);

        Matrix inverse() const;

        Matrix transpose() const;
        Matrix dot(const Matrix &B) const;
        Matrix operator*(const Matrix &B) const;
        Matrix operator*(const T &scal) const;
        Matrix operator+(const Matrix &B) const;
        Matrix operator+(const T &scal) const;
        Matrix& operator=(const Matrix &cop);

        // Why a virtual destructor? This class
        // has value semantics and should not be
        // inherited
        virtual ~Matrix();

    protected:
        int m_rows;
        int m_cols;
        // Why not std::vector?
        // You didn't implement copy constructor,
        // not move constructor and move assign operators
        // Having a pointer here without the right value
        // semantics implementation epxoses you to a lot
        // of potential crashes.
        T* m_data;

        void inverse_row(Matrix* M, int row, int col);
        void switch_row(int row1, int row2);
        void reduction(int row, T alpha);
        void transvection(int row1, int row2, T alpha);

};

#endif // MATRIX_H
