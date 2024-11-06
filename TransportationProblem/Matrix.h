#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cmath>
#include <stdexcept>

class Matrix {
public:
    Matrix(int rows, int columns);
    Matrix(const Matrix &m);

    int getRows() const;
    int getColumns() const;
    double getElement(int row, int column) const;
    void setElement(int row, int column, double value);

    Matrix transpose() const;
    Matrix subtract(const Matrix &other) const;
    Matrix add(const Matrix &other) const;
    static Matrix multiply(const Matrix &a, const Matrix &b);
    Matrix inverse() const;
    Matrix eliminateRow(int row) const;
    Matrix eliminateColumn(int column) const;

private:
    int rows;
    int columns;
    std::vector<std::vector<double> > data;

    std::vector<double> getRow(int idx) const;
    void setRow(int idx, const std::vector<double> &row);
};

#endif // MATRIX_H
