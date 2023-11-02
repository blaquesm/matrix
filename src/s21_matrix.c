#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int ret = OK;
  if (rows > 0 && columns > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix != NULL) {
      for (int c = 0; c < rows && ret == 0; c++) {
        result->matrix[c] = (double *)calloc(columns, sizeof(double));
        if (result->matrix[c] == NULL) ret = 1;
      }
    }
  } else {
    result->matrix = NULL;
    result->columns = 0;
    result->rows = 0;
    ret = ERROR_NULL;
  }

  return ret;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int c = 0; c < A->rows; c++) {
      free(A->matrix[c]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int flag = 1;
  int check = 0;
  if (error_matrix(A) == ERROR_NULL && error_matrix(B) == ERROR_NULL) flag = 0;
  if (A->columns == B->columns && A->rows == B->rows && flag == 1) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) < 1e-7) check++;
      }
    }
    if (check != A->columns * A->rows) flag = 0;
  } else
    flag = 0;
  return flag;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  if (error_matrix(A) == ERROR_NULL || error_matrix(B) == ERROR_NULL)
    flag = ERROR_NULL;
  else if (A->columns != B->columns || A->rows != B->rows)
    flag = ERROR_CALC;
  if (flag == 0) {
    s21_create_matrix(A->rows, A->columns, result);
    if (result->matrix != NULL) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      result->matrix = NULL;
      flag = ERROR_NULL;
    }
  } else
    result->matrix = NULL;

  return flag;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  if (error_matrix(A) == ERROR_NULL || error_matrix(B) == ERROR_NULL)
    flag = ERROR_NULL;
  else if ((A->columns != B->columns || A->rows != B->rows) && flag == 0)
    flag = ERROR_CALC;
  if (flag == 0) {
    s21_create_matrix(A->rows, A->columns, result);
    if (result != NULL) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else
      flag = ERROR_NULL;
  } else
    result->matrix = NULL;

  return flag;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int flag = 0;
  if (error_matrix(A) == ERROR_NULL) flag = ERROR_NULL;
  if (flag == 0 && !isnan(number)) {
    s21_create_matrix(A->rows, A->columns, result);
    if (result != NULL) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else
      flag = ERROR_NULL;
  } else
    result->matrix = NULL;

  return flag;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  if (error_matrix(A) == ERROR_NULL || error_matrix(B) == ERROR_NULL)
    flag = ERROR_NULL;
  if (A->columns == B->rows && flag == OK) {
    s21_create_matrix(A->rows, B->columns, result);
    if (result != NULL) {
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          for (int n = 0; n < A->columns; n++) {
            result->matrix[i][j] += A->matrix[i][n] * B->matrix[n][j];
          }
        }
      }
    } else
      flag = ERROR_NULL;
  } else {
    result->matrix = NULL;
    flag = ERROR_CALC;
  }

  return flag;
}
int s21_transpose(matrix_t *A, matrix_t *result) {
  int flag = 0;
  if (error_matrix(A) == ERROR_NULL) {
    flag = ERROR_NULL;
    result->matrix = NULL;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    if (result->matrix != NULL) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    }
  }
  return flag;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int flag = 0;
  matrix_t calc = {0};
  if (error_matrix(A) == ERROR_NULL)
    flag = ERROR_NULL;
  else if (A->rows == 1 && A->columns == 1 && flag == 0) {
    if (s21_create_matrix(A->rows, A->columns, result) != ERROR_NULL)
      result->matrix[0][0] = A->matrix[0][0];
  } else if (A->rows == A->columns && A->rows > 1 && flag == 0) {
    if (s21_create_matrix(A->rows, A->columns, result) != ERROR_NULL) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          coml(i, j, A, &calc);
          result->matrix[i][j] = pow(-1, i + j) * minor__(&calc);
          s21_remove_matrix(&calc);
        }
      }
    } else {
      result->matrix = NULL;
      flag = ERROR_NULL;
    }
  } else {
    result->matrix = NULL;
    flag = ERROR_CALC;
  }
  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int flag = 0;
  double res = 0.0;
  matrix_t det;
  if (error_matrix(A) == ERROR_NULL)
    flag = ERROR_NULL;
  else if (A->rows == 1 && A->columns == 1 && flag == 0)
    *result = A->matrix[0][0];
  else if (A->rows > 1 && A->columns > 1 && A->rows == A->columns) {
    s21_calc_complements(A, &det);
    for (int i = 0; i < A->rows; i++) res += A->matrix[0][i] * det.matrix[0][i];
    *result = res;
    s21_remove_matrix(&det);
  } else
    flag = ERROR_CALC;
  return flag;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int flag = 0;
  double res = 0.0;
  matrix_t comp = {0};
  matrix_t tran = {0};
  flag = s21_determinant(A, &res);
  if (flag == OK) {
    if (res != 0.0) {
      s21_calc_complements(A, &comp);
      s21_transpose(&comp, &tran);
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = tran.matrix[i][j] / res;
        }
      }
    } else {
      result->matrix = NULL;
      flag = ERROR_CALC;
    }
  } else
    result->matrix = NULL;
  s21_remove_matrix(&comp);
  s21_remove_matrix(&tran);
  return flag;
}

void coml(int row, int columns, matrix_t *A, matrix_t *result) {
  if (s21_create_matrix(A->rows - 1, A->columns - 1, result) != ERROR_NULL) {
    int s = 0;
    for (int i = 0; i < A->rows; i++) {
      if (i != row) {
        int t = 0;
        for (int j = 0; j < A->columns; j++) {
          if (j != columns) {
            result->matrix[s][t] = A->matrix[i][j];
            t++;
          }
        }
        s++;
      }
    }
  } else
    result->matrix = NULL;
}

double minor__(matrix_t *A) {
  double res = 0;
  matrix_t mat = {0};
  if (A->columns == A->rows && A->columns == 1)
    res = A->matrix[0][0];
  else if (A->columns == A->rows && A->columns == 2)
    res = determ_2(*A);
  else if (A->columns > 2 && A->rows > 2 && A->columns == A->rows) {
    for (int k = 0; k < A->columns; k++) {
      if (s21_create_matrix(A->rows - 1, A->columns - 1, &mat) == 0) {
        for (int i = 1; i < A->rows; i++) {
          int t = 0;
          for (int j = 0; j < A->rows; j++) {
            if (j != k) {
              mat.matrix[i - 1][t] = A->matrix[i][j];
              t++;
            }
          }
        }
        res += pow(-1, k + 2) * A->matrix[0][k] * minor__(&mat);
        s21_remove_matrix(&mat);
      }
    }
  }
  return res;
}

double determ_2(matrix_t det) {
  double res = 0;
  if (det.columns == 2 && det.rows == 2)
    res = det.matrix[0][0] * det.matrix[1][1] -
          det.matrix[1][0] * det.matrix[0][1];
  return res;
}

int error_matrix(matrix_t *c) {
  int ret = 0;
  if (c == NULL || c->matrix == NULL || c->columns < 1 || c->rows < 1)
    ret = ERROR_NULL;
  return ret;
}
