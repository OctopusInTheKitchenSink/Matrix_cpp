#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../matrix.h"

TEST(Constructor, Correct_matrix1) {
  S21Matrix m;
  EXPECT_EQ(m.get_matrix(), nullptr);
  EXPECT_EQ(m.get_rows(), 0);
  EXPECT_EQ(m.get_cols(), 0);
};

TEST(Constructor, Correct_matrix2) {
  S21Matrix m{3, 3};
  EXPECT_EQ(m.get_rows(), 3);
  EXPECT_EQ(m.get_cols(), 3);
};

TEST(Constructor, Correct_matrix4) {
  S21Matrix m{3, 4};
  EXPECT_EQ(m.get_rows(), 3);
  EXPECT_EQ(m.get_cols(), 4);
};

TEST(Constructor, Correct_matrix3) {
  S21Matrix m(4, 5);
  S21Matrix m1(m);
  EXPECT_EQ(m1.get_rows(), m.get_rows());
  EXPECT_EQ(m1.get_cols(), m.get_cols());
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 5; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], m.get_matrix()[i][j]);
    }
  }
};

TEST(Constructor, Incorrect_matrix1) {
  EXPECT_THROW(S21Matrix m(-3, 4), std::invalid_argument);
};

TEST(Constructor, Incorrect_matrix2) {
  EXPECT_THROW(S21Matrix m(3, -4), std::invalid_argument);
};

TEST(Set_matrix, setter1) {
  S21Matrix m(4, 5);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  int counter = 0;
  EXPECT_EQ(2, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[counter++]);
    }
  }
};

TEST(Set_matrix, setter2) {
  S21Matrix m(4, 5);
  double **buff = new double *[3]{0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    buff[i] = new double[2]{0};
    for (int j = 0; j < 2; ++j) {
      buff[i][j] = buff1[counter++];
    }
  }
  m.set_matrix(3, 2, buff);
  EXPECT_EQ(3, m.get_rows());
  EXPECT_EQ(2, m.get_cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[i][j]);
    }
  }
};

TEST(Set_matrix, setter3) {
  S21Matrix m(2, 2);
  double **buff = new double *[3]{0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    buff[i] = new double[3]{0};
    for (int j = 0; j < 3; ++j) {
      buff[i][j] = buff1[counter++];
    }
  }
  m.set_matrix(3, 3, buff);
  EXPECT_EQ(3, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[i][j]);
    }
  }
};

TEST(set_matrix, setter4) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  int counter = 0;
  EXPECT_EQ(2, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[counter++]);
    }
  }
};

TEST(set_matrix, setter5) {
  S21Matrix m(1, 1);
  double **buff = new double *[3]{0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    buff[i] = new double[3]{0};
    for (int j = 0; j < 3; ++j) {
      buff[i][j] = buff1[counter++];
    }
  }
  EXPECT_THROW(m.set_matrix(-2, 3, buff), std::invalid_argument);
};

TEST(set_matrix, setter6) {
  S21Matrix m(1, 1);
  double **buff = new double *[3]{0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    buff[i] = new double[3]{0};
    for (int j = 0; j < 3; ++j) {
      buff[i][j] = buff1[counter++];
    }
  }
  EXPECT_THROW(m.set_matrix(2, -3, buff), std::invalid_argument);
};

TEST(set_matrix, setter7) {
  S21Matrix m(1, 1);
  double **buff = nullptr;
  EXPECT_THROW(m.set_matrix(2, 3, buff), std::invalid_argument);
};

TEST(set_rows, setter5) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  m.set_rows(3);
  int counter = 0;
  EXPECT_EQ(3, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i != 2) {
        EXPECT_EQ(m.get_matrix()[i][j], buff[counter++]);
      } else {
        EXPECT_EQ(m.get_matrix()[i][j], 0);
      }
    }
  }
};

TEST(set_rows, setter7) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  m.set_rows(1);
  int counter = 0;
  EXPECT_EQ(1, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int j = 0; j < 3; ++j) {
    EXPECT_EQ(m.get_matrix()[0][j], buff[j]);
  }
};

TEST(set_rows, setter8) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  EXPECT_THROW(m.set_rows(-1), std::out_of_range);
};

TEST(set_cols, setter6) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(3, 2, buff);
  m.set_cols(3);
  int counter = 0;
  EXPECT_EQ(3, m.get_rows());
  EXPECT_EQ(3, m.get_cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (j != 2) {
        EXPECT_EQ(m.get_matrix()[i][j], buff[counter++]);
      } else {
        EXPECT_EQ(m.get_matrix()[i][j], 0);
      }
    }
  }
};

TEST(set_cols, setter9) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  m.set_cols(1);
  int counter = 0;
  EXPECT_EQ(2, m.get_rows());
  EXPECT_EQ(1, m.get_cols());
  EXPECT_EQ(m.get_matrix()[0][0], 1.0);
  EXPECT_EQ(m.get_matrix()[1][0], 4.0);
};

TEST(set_rows, setter10) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  EXPECT_THROW(m.set_cols(-1), std::out_of_range);
};

TEST(set_rows, setter11) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  EXPECT_THROW(m.set_matrix(-1, 2, buff), std::invalid_argument);
};

TEST(set_rows, setter12) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  EXPECT_THROW(m.set_matrix(2, -1, buff), std::invalid_argument);
};

TEST(set_cols, setter13) {
  S21Matrix m(1, 1);
  double *buff = {nullptr};
  EXPECT_THROW(m.set_matrix(1, 2, buff), std::invalid_argument);
};

TEST(EqMatrix, eq1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_EQ(m1.EqMatrix(m2), false);
};

TEST(EqMatrix, eq2) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff);
  EXPECT_EQ(m1.EqMatrix(m2), true);
};

TEST(EqMatrix, eq3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  EXPECT_EQ(m1.EqMatrix(m2), false);
};

TEST(SumMatrix, sum1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1.SumMatrix(m2), std::invalid_argument);
};

TEST(SumMatrix, sum2) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1.SumMatrix(m2), std::invalid_argument);
};

TEST(SumMatrix, sum3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  m1.SumMatrix(m2);
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], buff[counter] + buff1[counter]);
      counter++;
    }
  }
};

TEST(SubMatrix, sub1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1.SubMatrix(m2), std::invalid_argument);
};

TEST(SubMatrix, sub2) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1.SubMatrix(m2), std::invalid_argument);
};

TEST(SubMatrix, sub3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  m1.SubMatrix(m2);
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], buff[counter] - buff1[counter]);
      counter++;
    }
  }
};

TEST(MulNumber, mul_n1) {
  S21Matrix m(2, 3);
  m.MulNumber(4);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], 0);
    }
  }
};

TEST(MulNumber, mul_n2) {
  S21Matrix m(2, 3);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  m.MulNumber(4);
  int counter = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[counter++] * 4);
    }
  }
};

TEST(MulMatrix, mul_matr1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1.MulMatrix(m2), std::invalid_argument);
};

TEST(MulMatrix, mul_matr2) {
  S21Matrix m1(2, 1), m2(2, 1);
  EXPECT_THROW(m1.MulMatrix(m2), std::invalid_argument);
};

TEST(MulMatrix, mul_matr3) {
  S21Matrix m1(3, 3), m2(3, 3);
  double buff[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9},
         result[9] = {30, 36, 42, 66, 81, 96, 102, 126, 150};
  int counter = 0;
  m1.set_matrix(3, 3, buff), m2.set_matrix(3, 3, buff);
  m1.MulMatrix(m2);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(MulMatrix, mul_matr4) {
  S21Matrix m1(1, 1), m2(1, 1);
  double buff[] = {100}, buff1[] = {0.001};
  m1.set_matrix(1, 1, buff), m2.set_matrix(1, 1, buff1);
  m1.MulMatrix(m2);
  EXPECT_EQ(m1.get_matrix()[0][0], 0.1);
};

TEST(Transpose, tr1) {
  S21Matrix m1(3, 3);
  double buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9},
         result[] = {1, 4, 7, 2, 5, 8, 3, 6, 9};
  int counter = 0;
  m1.set_matrix(3, 3, buf);
  m1.Transpose();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(Transpose, tr2) {
  S21Matrix m1(1, 1);
  double buff[] = {0.001};
  m1.set_matrix(1, 1, buff);
  m1.Transpose();
  EXPECT_EQ(m1.get_matrix()[0][0], buff[0]);
};

TEST(Transopse, tr3) {
  S21Matrix m1(2, 2);
  double buff[] = {0.01, 0.2, 1, 2.1}, result[] = {0.01, 1, 0.2, 2.1};
  int counter = 0;
  m1.set_matrix(2, 2, buff);
  m1.Transpose();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(CalcCompl, compl1) {
  S21Matrix m1(2, 3);
  EXPECT_THROW(m1.CalcComplements(), std::domain_error);
};

TEST(CalcCompl, compl2) {
  S21Matrix m(3, 3);
  double buff[9] = {5, -1, 1, 2, 3, 4, 1, 0, 3},
         result[9] = {9, -2, -3, 3, 14, -1, -7, -18, 17};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  S21Matrix res = m.CalcComplements();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(CalcCompl, compl3) {
  S21Matrix m(5, 5);
  int counter = 0;
  double buff[25] = {1,  2,  3,  4, 5,  -1, 0, 2, 3, 1, -2, 1, 1,
                     -1, -3, -1, 2, -2, -3, 0, 5, 4, 3, 2,  1},
         result[25] = {-15,  -70, 197, -173, 110,  -48, 168, -330, 348,
                       -138, -30, -42, 198,  -150, 24,  -36, 126,  -174,
                       114,  -30, 33,  56,   -61,  67,  -46};
  m.set_matrix(5, 5, buff);
  S21Matrix res = m.CalcComplements();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      EXPECT_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(CalcCompl, compl4) {
  S21Matrix m(4, 4);
  int counter = 0;
  double buff[16] = {12, 11, 10, 9, 1, 2, 3, 1, -1, -2, 0, 4, 4, 3, 2, 1},
         result[16] = {22,  -41, 25,  -15, -16, 56,  -64, 24,
                       -24, 48,  -24, 0,   -86, 121, -65, 39};
  m.set_matrix(4, 4, buff);
  S21Matrix res = m.CalcComplements();
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(CalcCompl, compl5) {
  S21Matrix m(1, 1);
  double buff[1] = {12};
  m.set_matrix(1, 1, buff);
  S21Matrix res = m.CalcComplements();
  EXPECT_EQ(res.get_matrix()[0][0], 1);
};

TEST(Determinant, det1) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m.Determinant(), std::domain_error);
};

TEST(Determinant, det2) {
  S21Matrix m(2, 2);
  double buff[4] = {1, 2, 3, 4};
  int counter = 0;
  m.set_matrix(2, 2, buff);
  double res = m.Determinant();
  EXPECT_EQ(res, -2);
};

TEST(Determinant, det3) {
  S21Matrix m(3, 3);
  double buff[9] = {1, 2, 3, 4, 10, 6, 7, 8, 9};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  double res = m.Determinant();
  EXPECT_EQ(res, -60);
};

TEST(Determinant, det4) {
  S21Matrix m(3, 3);
  double buff[9] = {1, 2, 3, 4, 0, 6, 7, 8, 9};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  double res = m.Determinant();
  EXPECT_EQ(res, 60);
};

TEST(Determinant, det5) {
  S21Matrix m(3, 3);
  double buff[9] = {1, 0, 3, 4, 0, 6, 7, 8, 9};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  double res = m.Determinant();
  EXPECT_EQ(res, 48);
};

TEST(Determinant, det6) {
  S21Matrix m(7, 7);
  double buff[49] = {1, 0, 1, -1, 2, 3, 4, 1,  2,  3, 4, 5,  6, 7, 7, 6, 5,
                     4, 3, 2, 1,  0, 1, 2, -1, -2, 2, 0, 1,  0, 3, 4, 0, 6,
                     7, 1, 2, 3,  4, 0, 6, 7,  0,  1, 0, -1, 0, 1, 2};
  int counter = 0;
  m.set_matrix(7, 7, buff);
  double res = m.Determinant();
  EXPECT_EQ(res, -1120);
}

TEST(InverseMatrix, inv1) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m.InverseMatrix(), std::domain_error);
};

TEST(InverseMatrix, inv2) {
  S21Matrix m(2, 2);
  EXPECT_THROW(m.InverseMatrix(), std::runtime_error);
};

TEST(InverseMatrix, inv3) {
  S21Matrix m(3, 3);
  double buff[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3},
         result[9] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  S21Matrix res = m.InverseMatrix();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_DOUBLE_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
}

TEST(InverseMatrix, inv4) {
  S21Matrix m(2, 2);
  double buff[4] = {1, 2, 3, 4}, result[4] = {-2, 1, 1.5, -0.5};
  int counter = 0;
  m.set_matrix(2, 2, buff);
  S21Matrix res = m.InverseMatrix();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(InverseMatrix, inv5) {
  S21Matrix m(3, 3);
  double buff[9] = {2, 0, 1, 0, -3, -1, -2, 4, 0},
         result[9] = {2, 2, 1.5, 1, 1, 1, -3, -4, -3};
  int counter = 0;
  m.set_matrix(3, 3, buff);
  S21Matrix res = m.InverseMatrix();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(res.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(InverseMatrix, inv6) {
  S21Matrix m(3, 3);
  double buff[9] = {0, 0, 0, 1, 0, 0, 0, 0, 0};
  m.set_matrix(3, 3, buff);
  EXPECT_THROW(m.InverseMatrix(), std::runtime_error);
};

TEST(InverseMatrix, inv7) {
  S21Matrix m(3, 2);
  EXPECT_THROW(m.InverseMatrix(), std::domain_error);
};

TEST(oper_plus, op1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 + m2, std::invalid_argument);
};

TEST(oper_plus, op2) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1 + m2, std::invalid_argument);
};

TEST(oper_plus, op3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  S21Matrix m3 = m1 + m2;
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], buff[counter] + buff1[counter]);
      counter++;
    }
  }
};

TEST(op_minus, op4) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 - m2, std::invalid_argument);
};

TEST(op_minus, op5) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1 - m2, std::invalid_argument);
};

TEST(op_minus, op6) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  S21Matrix m3 = m1 - m2;
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], buff[counter] - buff1[counter]);
      counter++;
    }
  }
};

TEST(op_mult, op7) {
  S21Matrix m(2, 3);
  S21Matrix m3 = m * 4;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], 0);
    }
  }
}

TEST(op_mult, op8) {
  S21Matrix m(2, 3);
  S21Matrix m3 = m * 4;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], 0);
    }
  }
};

TEST(op_mult, op9) {
  S21Matrix m(2, 3);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  S21Matrix m3 = m * 4;
  int counter = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], buff[counter++] * 4);
    }
  }
};

TEST(op_mult, op10) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 * m2, std::invalid_argument);
};

TEST(op_mult, op11) {
  S21Matrix m1(2, 1), m2(2, 1);
  EXPECT_THROW(m1 * m2, std::invalid_argument);
};

TEST(op_mult, op12) {
  S21Matrix m1(3, 3), m2(3, 3);
  double buff[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9},
         result[9] = {30, 36, 42, 66, 81, 96, 102, 126, 150};
  int counter = 0;
  m1.set_matrix(3, 3, buff), m2.set_matrix(3, 3, buff);
  S21Matrix m3 = m1 * m2;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m3.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(op_mult, op13) {
  S21Matrix m1(1, 1), m2(1, 1);
  double buff[] = {100}, buff1[] = {0.001};
  m1.set_matrix(1, 1, buff), m2.set_matrix(1, 1, buff1);
  S21Matrix m3 = m1 * m2;
  EXPECT_EQ(m3.get_matrix()[0][0], 0.1);
};

TEST(eq_op, eq_op1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_EQ(m1 == m2, false);
}

TEST(eq_op, eq_op2) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_EQ(m1 == m2, false);
};

TEST(eq_op, eq_op3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff);
  EXPECT_EQ(m1 == m2, true);
};

TEST(eq_op, eq_op4) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  EXPECT_EQ(m1 == m2, false);
};

TEST(sum_eq_op, sum_eq_op1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 += m2, std::invalid_argument);
};

TEST(sum_eq_op, sum_eq_op2) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1 += m2, std::invalid_argument);
};

TEST(sum_eq_op, sum_eq_op3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  m1 += m2;
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], buff[counter] + buff1[counter]);
      counter++;
    }
  }
};

TEST(sub_eq_op, sub_eq_op1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 -= m2, std::invalid_argument);
};

TEST(sub_eq_op, sub_eq_op2) {
  S21Matrix m1(2, 1), m2(2, 2);
  EXPECT_THROW(m1 -= m2, std::invalid_argument);
};

TEST(sub_eq_op, sub_eq_op3) {
  S21Matrix m1(3, 2), m2(3, 2);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double buff1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 7.0};
  m1.set_matrix(3, 2, buff), m2.set_matrix(3, 2, buff1);
  m1 -= m2;
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], buff[counter] - buff1[counter]);
      counter++;
    }
  }
};

TEST(mul_n_eq_op, mul_n_eq_op1) {
  S21Matrix m(2, 3);
  m *= 4;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], 0);
    }
  }
};

TEST(mul_n_eq_op, mul_n_eq_op2) {
  S21Matrix m(2, 3);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  m *= 4;
  int counter = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m.get_matrix()[i][j], buff[counter++] * 4);
    }
  }
};

TEST(mul_matr_eq_op, mul_matr_eq_op1) {
  S21Matrix m1(1, 1), m2(2, 2);
  EXPECT_THROW(m1 *= m2, std::invalid_argument);
};

TEST(mul_matr_eq_op, mul_matr_eq_op2) {
  S21Matrix m1(2, 1), m2(2, 1);
  EXPECT_THROW(m1 *= m2, std::invalid_argument);
};

TEST(mul_matr_eq_op, mul_matr_eq_op3) {
  S21Matrix m1(3, 3), m2(3, 3);
  double buff[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9},
         result[9] = {30, 36, 42, 66, 81, 96, 102, 126, 150};
  int counter = 0;
  m1.set_matrix(3, 3, buff), m2.set_matrix(3, 3, buff);
  m1 *= m2;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m1.get_matrix()[i][j], result[counter++]);
    }
  }
};

TEST(mul_matr_eq_op, mul_matr_eq_op4) {
  S21Matrix m1(1, 1), m2(1, 1);
  double buff[] = {100}, buff1[] = {0.001};
  m1.set_matrix(1, 1, buff), m2.set_matrix(1, 1, buff1);
  m1 *= m2;
  EXPECT_EQ(m1.get_matrix()[0][0], 0.1);
};

TEST(assignment_op, test1) {
  S21Matrix m1(1, 1);
  S21Matrix m2 = m1;
  EXPECT_EQ(m1.get_cols(), m2.get_cols());
  EXPECT_EQ(m1.get_rows(), m2.get_rows());
}

TEST(parenthesis, par1) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m(5, 0), std::out_of_range);
};

TEST(parenthesis, par4) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m(0, 5), std::out_of_range);
};

TEST(parenthesis, par2) {
  S21Matrix m(2, 3);
  EXPECT_EQ(m(1, 0), 0);
};

TEST(parenthesis, par3) {
  S21Matrix m(1, 1);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  EXPECT_EQ(m(0, 1), 2.0);
};

TEST(parenthesis, par5) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m(0, -5), std::out_of_range);
};

TEST(parenthesis, par6) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m(-1, 0), std::out_of_range);
};

TEST(parenthesis, par7) {
  S21Matrix m(2, 3);
  EXPECT_THROW(m(2, 3), std::out_of_range);
  EXPECT_THROW(m(0, 3), std::out_of_range);
  EXPECT_THROW(m(2, 0), std::out_of_range);
};

TEST(parenthesis, par8) {
  S21Matrix m(1, 1);
  EXPECT_THROW(m(0, 1), std::out_of_range);
};

TEST(parenthesis, par9) {
  S21Matrix m(1, 1);
  EXPECT_THROW(m(1, 0), std::out_of_range);
};

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
};

TEST(MoveConstructor, move1) {
  S21Matrix m(2, 3);
  double buff[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  m.set_matrix(2, 3, buff);
  double **original_m = m.get_matrix();
  S21Matrix moved(std::move(m));
  EXPECT_EQ(moved.get_rows(), 2);
  EXPECT_EQ(moved.get_cols(), 3);
  EXPECT_EQ(moved.get_matrix(), original_m);
  int counter = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_DOUBLE_EQ(moved.get_matrix()[i][j], buff[counter++]);
    }
  }
};

TEST(MoveConstructor, move2) {
  S21Matrix m(3, 2);
  double buff[] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
  m.set_matrix(3, 2, buff);
  S21Matrix moved(std::move(m));
  EXPECT_EQ(moved.get_rows(), 3);
  EXPECT_EQ(moved.get_cols(), 2);
  EXPECT_EQ(m.get_matrix(), nullptr);
};

TEST(AssignmentOperator, ass1) {
  S21Matrix m(2, 2);
  double buff[] = {1, 2, 3, 4};
  m.set_matrix(2, 2, buff);
  m = m;
  int counter = 0;
  EXPECT_EQ(m.get_rows(), 2);
  EXPECT_EQ(m.get_cols(), 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_DOUBLE_EQ(m.get_matrix()[i][j], buff[counter++]);
    }
  }
};

TEST(AssignmentOperator, ass2) {
  S21Matrix m1(3, 2);
  double buff[] = {1, 2, 3, 4, 5, 6};
  m1.set_matrix(3, 2, buff);
  S21Matrix m2(1, 1);
  m2 = m1;
  EXPECT_EQ(m2.get_rows(), 3);
  EXPECT_EQ(m2.get_cols(), 2);
  EXPECT_DOUBLE_EQ(m2(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m2(2, 1), 6.0);
};

TEST(SetMatrixTest, CorrectDataCopy) {
  double **src = new double *[2];
  src[0] = new double[2]{1.0, 2.0};
  src[1] = new double[2]{3.0, 4.0};

  S21Matrix mat;
  mat.set_matrix(2, 2, src);

  EXPECT_EQ(mat(0, 0), 1.0);
  EXPECT_EQ(mat(0, 1), 2.0);
  EXPECT_EQ(mat(1, 0), 3.0);
  EXPECT_EQ(mat(1, 1), 4.0);

  // Очистка
  for (int i = 0; i < 2; ++i) delete[] src[i];
  delete[] src;
};

TEST(SetMatrixTest, ZeroRows) {
  S21Matrix mat;
  double **src = new double *[1];
  EXPECT_THROW(mat.set_matrix(0, 1, src), std::invalid_argument);
  delete[] src;
}

TEST(SetMatrixTest, NegativeCols) {
  S21Matrix mat;
  double **src = new double *[1];
  EXPECT_THROW(mat.set_matrix(1, -1, src), std::invalid_argument);
  delete[] src;
}

TEST(SetMatrixTest, InvalidDimensions) {
  S21Matrix mat;
  double data[4] = {1, 2, 3, 4};

  EXPECT_THROW(mat.set_matrix(0, 2, data), std::invalid_argument);
  EXPECT_THROW(mat.set_matrix(2, -1, data), std::invalid_argument);
}

TEST(SetMatrixTest, CorrectDataCopy1) {
  S21Matrix mat;
  double data[4] = {1.1, 2.2, 3.3, 4.4};

  mat.set_matrix(2, 2, data);

  EXPECT_DOUBLE_EQ(mat(0, 0), 1.1);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.2);
  EXPECT_DOUBLE_EQ(mat(1, 0), 3.3);
  EXPECT_DOUBLE_EQ(mat(1, 1), 4.4);
}

TEST(SetMatrixTest, MemoryDeallocation) {
  S21Matrix mat;
  double data1[4] = {1, 2, 3, 4};
  double data2[4] = {5, 6, 7, 8};

  mat.set_matrix(2, 2, data1);
  mat.set_matrix(2, 2, data2);

  EXPECT_DOUBLE_EQ(mat(1, 1), 8.0);
}

TEST(SetMatrixTest, LinearToMatrixConversion) {
  S21Matrix mat;
  double data[6] = {1, 2, 3, 4, 5, 6};

  mat.set_matrix(2, 3, data);

  EXPECT_DOUBLE_EQ(mat(0, 0), 1);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2);
  EXPECT_DOUBLE_EQ(mat(0, 2), 3);
  EXPECT_DOUBLE_EQ(mat(1, 0), 4);
  EXPECT_DOUBLE_EQ(mat(1, 1), 5);
  EXPECT_DOUBLE_EQ(mat(1, 2), 6);
}

TEST(SetRowsTest, InvalidRowsLessThanOne) {
  S21Matrix mat(2, 2);
  EXPECT_THROW(mat.set_rows(0), std::out_of_range);
  EXPECT_THROW(mat.set_rows(-1), std::out_of_range);
}

TEST(SetRowsTest, IncreaseRows) {
  S21Matrix mat(2, 2);
  double data[4] = {1, 2, 3, 4};
  mat.set_matrix(2, 2, data);

  mat.set_rows(3);

  EXPECT_EQ(mat.get_rows(), 3);
  EXPECT_EQ(mat.get_cols(), 2);

  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);

  EXPECT_DOUBLE_EQ(mat(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 0.0);
}

TEST(SetRowsTest, DecreaseRows) {
  S21Matrix mat(3, 2);
  double data[6] = {1, 2, 3, 4, 5, 6};
  mat.set_matrix(3, 2, data);

  mat.set_rows(2);

  EXPECT_EQ(mat.get_rows(), 2);
  EXPECT_EQ(mat.get_cols(), 2);

  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);
}

TEST(SetRowsTest, SameRowsCount) {
  S21Matrix mat(2, 2);
  double data[4] = {1, 2, 3, 4};
  mat.set_matrix(2, 2, data);

  mat.set_rows(2);

  EXPECT_EQ(mat.get_rows(), 2);
  EXPECT_EQ(mat.get_cols(), 2);
  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);
}

TEST(SetRowsTest, MemoryReallocation) {
  S21Matrix mat(2, 2);
  double *old_matrix_data = mat.get_matrix()[0];

  mat.set_rows(3);
  double *new_matrix_data = mat.get_matrix()[0];

  EXPECT_NE(old_matrix_data, new_matrix_data);

  EXPECT_EQ(mat.get_rows(), 3);
}

TEST(SetRowsTest, FromEmptyMatrix) {
  S21Matrix mat;

  EXPECT_THROW(mat.set_rows(2), std::out_of_range);
}

TEST(SetRowsTest, MultipleResizes) {
  S21Matrix mat(2, 2);
  double data[4] = {1, 2, 3, 4};
  mat.set_matrix(2, 2, data);
  mat.set_rows(3);
  EXPECT_EQ(mat.get_rows(), 3);
  mat.set_rows(1);
  EXPECT_EQ(mat.get_rows(), 1);
  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
  mat.set_rows(2);
  EXPECT_EQ(mat.get_rows(), 2);
  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 0.0);
};

TEST(S2MatrixTest, SetRowsLarge) {
  S21Matrix m(1, 3);
  m.set_rows(1000);
  EXPECT_EQ(m.get_rows(), 1000);
  EXPECT_EQ(m.get_cols(), 3);
  EXPECT_EQ(m(0, 0), 0.0);
  EXPECT_EQ(m(0, 1), 0.0);
  EXPECT_EQ(m(0, 2), 0.0);
  EXPECT_EQ(m(1, 0), 0.0);
  EXPECT_EQ(m(1, 1), 0.0);
  EXPECT_EQ(m(1, 2), 0.0);
  EXPECT_EQ(m(999, 0), 0.0);
  EXPECT_EQ(m(999, 1), 0.0);
  EXPECT_EQ(m(999, 2), 0.0);
}