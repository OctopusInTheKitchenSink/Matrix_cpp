FLAGS = -lgtest -lgmock -pthread -g -O0 -fprofile-arcs -ftest-coverage

all:
	g++ -Wall -Werror -Wextra -std=c++20 matrix.cpp

clean:
	rm a.out test *.gcov *.css *.html *.gcda *.gcno tests/*html tests/test tests/*.gcda tests/*.gcno tests/*css *.o *.a

test:
	g++ tests/s21_test.cpp -o tests/test matrix.cpp $(FLAGS)
	./tests/test

gcov: test
	gcov tests/*.gcno tests/*.gcda 
	find ./tests/ ! -name "s21_test.cpp" ! -name "*.html" ! -name "*matrix*" -type f -delete
	gcovr -r . --html --html-details -o tests/coverage_report.html
	open tests/coverage_report.html

s21_matrix_oop.a:
	g++ -c matrix.cpp -o matrix.o
	ar -crs s21_matrix_oop.a matrix.o