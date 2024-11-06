#include <iostream>
#include "Matrix.h"

class InputReader {
public:
    InputReader();

    int readInt();
    void readMatrixNx1(Matrix &column);
    void readMatrixNxM(Matrix &matrix);

private:
    std::istream &input;
};

InputReader::InputReader() : input(std::cin) {}

/**
 * Reads integer value from input.
 * @return integer value read from input.
 */
int InputReader::readInt() {
    int value;
    input >> value;
    return value;
}

/**
 * Reads column vector (Nx1 matrix) from input and puts inside provided buffer 'column'.
 * @param column buffer for the column vector (Nx1 matrix) read from input.
 */
void InputReader::readMatrixNx1(Matrix &column) {
    int rows = column.getRows();
    for (int i = 0; i < rows; ++i) {
        double value;
        input >> value;
        column.setElement(i, 0, value);
    }
}

/**
 * Reads NxM matrix from input and puts inside provided buffer 'matrix'.
 * @param matrix buffer for the NxM matrix read from input.
 */
void InputReader::readMatrixNxM(Matrix &matrix) {
    int rows = matrix.getRows();
    int cols = matrix.getColumns();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value;
            input >> value;
            matrix.setElement(i, j, value);
        }
    }
}
