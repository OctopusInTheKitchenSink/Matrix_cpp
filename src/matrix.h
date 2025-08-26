#ifndef H_MATRIX
#define H_MATRIX

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int rows_{0}, cols_{0};
  double **matrix_{nullptr};

 public:
  S21Matrix(int r, int c);
  S21Matrix() noexcept;
  S21Matrix(const S21Matrix &other) noexcept;
  ~S21Matrix();
  S21Matrix(S21Matrix &&other) noexcept;
  int get_rows() const;
  int get_cols() const;
  double **get_matrix() const;
  void set_rows(int r);
  void set_cols(int c);
  void set_matrix(int rows, int cols, double **matrix);
  void set_matrix(int rows, int cols, double *matrix);
  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix &Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
  void round_for_inverse();
  void get_det(double &res);
  int zeros_counter(char &param);
  void minor_matrix(int row, int col, S21Matrix &matrix);
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double &num);
  bool operator==(const S21Matrix &other);
  S21Matrix &operator=(const S21Matrix &other) noexcept;
  S21Matrix &operator=(const S21Matrix &&other) noexcept;
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double &num);
  double operator()(int rows, int cols) const;
};

#endif