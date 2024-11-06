#include "Matrix.h"

//#include <vector>
//#include <cmath>
//#include <stdexcept>
//
//class Matrix {
//public:
//    Matrix(int rows, int columns);
//    Matrix(const Matrix &m);
//
//    int getRows() const;
//    int getColumns() const;
//    double getElement(int row, int column) const;
//    void setElement(int row, int column, double value);
//
//    Matrix transpose() const;
//    Matrix subtract(const Matrix &other) const;
//    Matrix add(const Matrix &other) const;
//    static Matrix multiply(const Matrix &a, const Matrix &b);
//    Matrix inverse() const;
//    Matrix eliminateRow(int row) const;
//    Matrix eliminateColumn(int column) const;
//
//private:
//    int rows;
//    int columns;
//    std::vector<std::vector<double> > data;
//
//    std::vector<double> getRow(int idx) const;
//    void setRow(int idx, const std::vector<double> &row);
//};

Matrix::Matrix(int rows, int columns) : rows(rows), columns(columns), data(rows, std::vector<double>(columns, 0.0)) {}

Matrix::Matrix(const Matrix &m) : rows(m.rows), columns(m.columns), data(m.data) {}

int Matrix::getRows() const {
    return rows;
}

int Matrix::getColumns() const {
    return columns;
}

double Matrix::getElement(int row, int column) const {
    return data[row][column];
}

void Matrix::setElement(int row, int column, double value) {
    data[row][column] = value;
}

std::vector<double> Matrix::getRow(int idx) const {
    return data[idx];
}

void Matrix::setRow(int idx, const std::vector<double> &row) {
    data[idx] = row;
}

Matrix Matrix::transpose() const {
    Matrix result(columns, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            result.setElement(j, i, getElement(i, j));
        }
    }
    return result;
}

Matrix Matrix::subtract(const Matrix &other) const {
    Matrix result(rows, columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            result.setElement(i, j, getElement(i, j) - other.getElement(i, j));
        }
    }
    return result;
}

Matrix Matrix::add(const Matrix &other) const {
    Matrix result(rows, columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            result.setElement(i, j, getElement(i, j) + other.getElement(i, j));
        }
    }
    return result;
}

Matrix Matrix::multiply(const Matrix &a, const Matrix &b) {
    if (a.columns != b.rows) throw std::invalid_argument("Matrix dimensions must match for multiplication.");

    Matrix result(a.rows, b.columns);
    for (int i = 0; i < a.rows; ++i) {
        for (int j = 0; j < b.columns; ++j) {
            double sum = 0;
            for (int k = 0; k < a.columns; ++k) {
                sum += a.getElement(i, k) * b.getElement(k, j);
            }
            result.setElement(i, j, sum);
        }
    }
    return result;
}

Matrix Matrix::inverse() const {
    Matrix identity(columns, rows);
    for (int i = 0; i < rows; ++i) identity.setElement(i, i, 1);

    Matrix output(identity);
    Matrix input(*this);

    for (int j = 0; j < columns; ++j) {
        if (std::fabs(input.getElement(j, j)) < 1e-10) input.setElement(j, j, 0);

        Matrix p(identity);
        double maxElem = input.getElement(j, j);
        int idx = j;

        for (int k = j; k < rows; ++k) {
            if (std::fabs(input.getElement(k, j)) > maxElem) {
                maxElem = std::fabs(input.getElement(k, j));
                idx = k;
            }
        }

        if (idx != j) {
            std::vector<double> temp = p.getRow(idx);
            p.setRow(idx, p.getRow(j));
            p.setRow(j, temp);
            input = multiply(p, input);
            output = multiply(p, output);
        }

        for (int i = j + 1; i < rows; ++i) {
            if (std::fabs(input.getElement(i, j)) < 1e-10) input.setElement(i, j, 0);
            if (input.getElement(i, j) == 0) continue;

            Matrix e(identity);
            e.setElement(i, j, -input.getElement(i, j) / input.getElement(j, j));
            input = multiply(e, input);
            output = multiply(e, output);
        }
    }

    for (int j = columns - 1; j >= 0; j--) {
        for (int i = rows - 1; i >= 0; i--) {
            if (std::fabs(input.getElement(i, j)) < 1e-10) input.setElement(i, j, 0);
            if (input.getElement(i, j) == 0) continue;

            Matrix e(identity);
            e.setElement(i, j, -input.getElement(i, j) / input.getElement(j, j));
            input = multiply(e, input);
            output = multiply(e, output);
        }
    }

    Matrix scale(identity);
    for (int i = 0; i < rows; ++i) scale.setElement(i, i, 1 / input.getElement(i, i));

    return multiply(scale, output);
}

Matrix Matrix::eliminateRow(int row) const {
    Matrix result(rows - 1, columns);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < columns; ++j) {
            result.setElement(i, j, getElement(i, j));
        }
    }
    for (int i = row + 1; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            result.setElement(i - 1, j, getElement(i, j));
        }
    }
    return result;
}

Matrix Matrix::eliminateColumn(int column) const {
    Matrix result(rows, columns - 1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < column; ++j) {
            result.setElement(i, j, getElement(i, j));
        }
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = column + 1; j < columns; ++j) {
            result.setElement(i, j - 1, getElement(i, j));
        }
    }
    return result;
}
