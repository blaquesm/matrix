.PHONY: all s21_matrix.a test gcov_report check rebuild clean unit_tests.o
CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11
GCOVFLAGS =  -fprofile-arcs -ftest-coverage 
LC = lcov -c -d .
C_FILES = s21*.c
O_FILES = s21*.o
LIBS=-lcheck -lm
LINUX=-lcheck  -lm -lpthread -lrt -lsubunit
OS=$(shell uname -s)

all : s21_matrix.a	
clean : 
	rm -rf *.o *.gcno *.gcda *.info test.out s21_matrix.a
	rm -rf ./report
	rm -rf *.dSYM
rebuild : clean all

$(O_FILES) :
	$(CC) $(CFLAGS) -c $(C_FILES) s21_matrix.h
unit_tests.o :
	$(CC) $(CFLAGS) -c -lcheck unit_tests.c
s21_matrix.a :
	$(CC) $(CFLAGS) -c $(C_FILES)
	ar r s21_matrix.a $(O_FILES)

test: s21_matrix.a unit_tests.c
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) unit_tests.c s21_matrix.a -o test.out $(LIBS)
else
	$(CC) $(CFLAGS) unit_tests.c s21_matrix.a -o test.out $(LIBS) $(LINUX)
endif
	./test.out	
gcov_report : clean
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) $(GCOVFLAGS) -o gcov_test $(C_FILES) unit_tests.c $(LIBS)
else
	$(CC) $(CFLAGS) $(GCOVFLAGS) -o gcov_test  $(C_FILES) unit_tests.c $(LIBS) $(LINUX)
endif
	./gcov_test
	lcov -t "gcov_report" -o report.info --no-external -c -d .
	genhtml -o report/ report.info
	rm -f *.gcno *.gcda *.info report.out *.gcov
	open ./report/index.html


check:
	cp ../materials/linters/.clang-format ./
	clang-format -n *.c
	clang-format -n *.h
	rm -rf .clang-format
