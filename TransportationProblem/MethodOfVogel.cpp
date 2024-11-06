#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include "Matrix.h"

class MethodOfVogel {
private:
    Matrix supplyCoefficients;
    Matrix costsCoefficients;
    Matrix demandCoefficients;

public:
    MethodOfVogel(const Matrix &coefficientsOfSupply, const Matrix &coefficientsOfCosts, const Matrix &coefficientsOfDemand)
            : supplyCoefficients(coefficientsOfSupply), costsCoefficients(coefficientsOfCosts), demandCoefficients(coefficientsOfDemand) {}

private:
    int firstMin(int column, int row) {
        if (row == -1 && column == -1) return -1;

        double firstMin = std::numeric_limits<double>::max();
        int index = -1;
        if (column == -1) {
            for (int i = 0; i < costsCoefficients.getColumns(); ++i) {
                double element = costsCoefficients.getElement(row, i);
                if (element < firstMin) {
                    firstMin = element;
                    index = i;
                }
            }
        } else {
            for (int i = 0; i < costsCoefficients.getRows(); ++i) {
                double element = costsCoefficients.getElement(i, column);
                if (element < firstMin) {
                    firstMin = element;
                    index = i;
                }
            }
        }
        return index;
    }

    int secondMin(int column, int row, int firstMinIdx) {
        if (row == -1 && column == -1) return -1;

        double secondMin = std::numeric_limits<double>::max();
        int index = -1;
        if (column == -1) {
            for (int i = 0; i < costsCoefficients.getColumns(); ++i) {
                if (i == firstMinIdx) continue;

                double element = costsCoefficients.getElement(row, i);
                if (element < secondMin) {
                    secondMin = element;
                    index = i;
                }
            }
        } else {
            for (int i = 0; i < costsCoefficients.getRows(); ++i) {
                if (i == firstMinIdx) continue;

                double element = costsCoefficients.getElement(i, column);
                if (element < secondMin) {
                    secondMin = element;
                    index = i;
                }
            }
        }
        return index;
    }

    std::vector<double> getRowDifferences() {
        int size = costsCoefficients.getRows();
        std::vector<double> differences(size, 0.0);

        for (int row = 0; row < size; ++row) {
            if (supplyCoefficients.getElement(row, 0) == 0) {
                differences[row] = std::numeric_limits<double>::max();
                continue;
            }

            int first = firstMin(-1, row);
            int second = secondMin(-1, row, first);

            if (first == -1)
                differences[row] = std::numeric_limits<double>::max();
            else if (second == -1)
                differences[row] = 0.0;
            else
                differences[row] = costsCoefficients.getElement(row, second) - costsCoefficients.getElement(row, first);
        }
        return differences;
    }

    std::vector<double> getColumnDifferences() {
        int size = costsCoefficients.getColumns();
        std::vector<double> differences(size, 0.0);

        for (int column = 0; column < size; ++column) {
            if (demandCoefficients.getElement(column, 0) == 0) {
                differences[column] = std::numeric_limits<double>::max();
                continue;
            }

            int first = firstMin(column, -1);
            int second = secondMin(column, -1, first);

            if (first == -1)
                differences[column] = std::numeric_limits<double>::max();
            else if (second == -1)
                differences[column] = 0.0;
            else
                differences[column] = costsCoefficients.getElement(second, column) - costsCoefficients.getElement(first, column);
        }
        return differences;
    }

    int getIndexOfMaxElement(const std::vector<double> &array) {
        double maxElement = -1.0;
        int index = -1;

        for (int i = 0; i < array.size(); ++i) {
            double element = array[i];
            if (element == std::numeric_limits<double>::max()) continue;

            if (element > maxElement) {
                maxElement = element;
                index = i;
            }
        }
        return index;
    }

    std::vector<int> maxDifference() {
        std::vector<double> rowDifferences = getRowDifferences();
        std::vector<double> columnDifferences = getColumnDifferences();

        int rowIndex = getIndexOfMaxElement(rowDifferences);
        int columnIndex = getIndexOfMaxElement(columnDifferences);

        if (rowIndex == -1 && columnIndex == -1) {
            return std::vector<int>(2, -1);  // Creates a vector of size 2 with both elements set to -1
        } else if (rowIndex == -1) {
            std::vector<int> result(2);  // Create a vector of size 2
            result[0] = 1;                 // Set first element to 1
            result[1] = columnIndex;       // Set second element to columnIndex
            return result;
        } else if (columnIndex == -1) {
            std::vector<int> result(2);  // Create a vector of size 2
            result[0] = 0;                 // Set first element to 0
            result[1] = rowIndex;          // Set second element to rowIndex
            return result;
        }

        // Compare rowDifferences and columnDifferences to determine the result
        std::vector<int> result(2);  // Create a vector of size 2
        if (rowDifferences[rowIndex] >= columnDifferences[columnIndex]) {
            result[0] = 0;               // Set first element to 0
            result[1] = rowIndex;        // Set second element to rowIndex
        } else {
            result[0] = 1;               // Set first element to 1
            result[1] = columnIndex;     // Set second element to columnIndex
        }
        return result;  // Return the result vector
    }




    int minCoefficientIdx(int column, int row) {
        return firstMin(column, row);
    }

    double min(double a, double b) {
        return std::min(a, b);
    }

    void blockEmptyRowAndColumnByPoint(int row, int column) {
        if (supplyCoefficients.getElement(row, 0) == 0) {
            for (int i = 0; i < costsCoefficients.getColumns(); ++i)
                costsCoefficients.setElement(row, i, std::numeric_limits<double>::max());
        }
        if (demandCoefficients.getElement(column, 0) == 0) {
            for (int i = 0; i < costsCoefficients.getRows(); ++i)
                costsCoefficients.setElement(i, column, std::numeric_limits<double>::max());
        }
    }

public:
    Matrix solve() {
        Matrix initialSolution(costsCoefficients.getRows(), costsCoefficients.getColumns());
        for (int i = 0; i < costsCoefficients.getRows(); ++i) {
            for (int j = 0; j < costsCoefficients.getColumns(); ++j) {
                initialSolution.setElement(i, j, 0);
            }
        }

        while (true) {
            std::vector<int> place = maxDifference();
            int type = place[0];
            if (type == -1) break;

            int row, column;
            if (type == 0) {
                row = place[1];
                column = minCoefficientIdx(-1, row);
            } else {
                column = place[1];
                row = minCoefficientIdx(column, -1);
            }

            double supply = supplyCoefficients.getElement(row, 0);
            double demand = demandCoefficients.getElement(column, 0);
            double minValue = min(supply, demand);

            initialSolution.setElement(row, column, minValue);

            costsCoefficients.setElement(row, column, std::numeric_limits<double>::max());
            supplyCoefficients.setElement(row, 0, supply - minValue);
            demandCoefficients.setElement(column, 0, demand - minValue);
            blockEmptyRowAndColumnByPoint(row, column);
        }

        return initialSolution;
    }
};
