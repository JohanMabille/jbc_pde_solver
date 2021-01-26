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
		void display() const;

        T elem_at(int i, int j) const;
        Matrix* set_elem_at(int i, int j, T val);

		Matrix inverse();
		void switch_row(int row1, int row2);
		void reduction(int row, double alpha);
		void transvection(int row1, int row2, double alpha);  
		void zeros(int row, int col);
		void inverse_row(Matrix<T>* M, int row, int col);
		bool is_still_inversible(int row, int col);

        Matrix transpose() const;
        Matrix dot(const Matrix &B) const;
		Matrix column(int j);
		void fill_column(int j, const Matrix data);
        Matrix operator*(const Matrix &B) const;
        Matrix operator*(const T &scal) const;
        Matrix operator+(const Matrix &B) const;
        Matrix operator+(const T &scal) const;
		T operator()(int i, int j) const; 

		//void test(); 
        virtual ~Matrix();

    protected:

    private:
        int m_rows;
        int m_cols;
        T* m_data;

};

#endif // MATRIX_H
