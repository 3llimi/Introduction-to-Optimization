#include <limits>
#include "Matrix.h"
#include <vector>

class MethodOfRussell {
private:
    Matrix supplyCoefficients;
    Matrix costsCoefficients;
    Matrix demandCoefficients;

    double maxColumnCoefficient(int column, const std::vector<int> &ignoredRows);
    double maxRowCoefficient(int row, const std::vector<int> &ignoredColumns);
    Matrix deltasMatrix(const std::vector<double> &rowCoefficients, const std::vector<double> &columnCoefficients,
                        const std::vector<int> &ignoredRows, const std::vector<int> &ignoredColumns);
    std::vector<int> largestNegativeDeltaCoordinates(const Matrix &deltas);
    double min(double a, double b);

public:
    MethodOfRussell(const Matrix &coefficientsOfSupply, const Matrix &coefficientsOfCosts, const Matrix &coefficientsOfDemand);
    Matrix solve();
};


MethodOfRussell::MethodOfRussell(const Matrix &coefficientsOfSupply, const Matrix &coefficientsOfCosts, const Matrix &coefficientsOfDemand)
        : supplyCoefficients(coefficientsOfSupply), costsCoefficients(coefficientsOfCosts), demandCoefficients(coefficientsOfDemand) {}

double MethodOfRussell::maxColumnCoefficient(int column, const std::vector<int> &ignoredRows) {
    double max = 0;
    for (int i = 0; i < costsCoefficients.getRows(); ++i) {
        if (std::find(ignoredRows.begin(), ignoredRows.end(), i) == ignoredRows.end() && costsCoefficients.getElement(i, column) > max) {
            max = costsCoefficients.getElement(i, column);
        }
    }
    return max;
}

double MethodOfRussell::maxRowCoefficient(int row, const std::vector<int> &ignoredColumns) {
    double max = 0;
    for (int i = 0; i < costsCoefficients.getColumns(); ++i) {
        if (std::find(ignoredColumns.begin(), ignoredColumns.end(), i) == ignoredColumns.end() && costsCoefficients.getElement(row, i) > max) {
            max = costsCoefficients.getElement(row, i);
        }
    }
    return max;
}

Matrix MethodOfRussell::deltasMatrix(const std::vector<double> &rowCoefficients, const std::vector<double> &columnCoefficients,
                                     const std::vector<int> &ignoredRows, const std::vector<int> &ignoredColumns) {
    Matrix deltas(costsCoefficients.getRows(), costsCoefficients.getColumns());
    for (size_t i = 0; i < rowCoefficients.size(); ++i) {
        for (size_t j = 0; j < columnCoefficients.size(); ++j) {
            if (std::find(ignoredRows.begin(), ignoredRows.end(), i) == ignoredRows.end() &&
                std::find(ignoredColumns.begin(), ignoredColumns.end(), j) == ignoredColumns.end()) {
                deltas.setElement(i, j, costsCoefficients.getElement(i, j) - rowCoefficients[i] - columnCoefficients[j]);
            }
        }
    }
    return deltas;
}

std::vector<int> MethodOfRussell::largestNegativeDeltaCoordinates(const Matrix &deltas) {
    std::vector<int> coordinates(2);  // Initialize with size 2
    coordinates[0] = 0;
    coordinates[1] = 0;

    double min = 0;
    for (int i = 0; i < deltas.getRows(); ++i) {
        for (int j = 0; j < deltas.getColumns(); ++j) {
            if (deltas.getElement(i, j) < min) {
                min = deltas.getElement(i, j);
                coordinates[0] = i;
                coordinates[1] = j;
            }
        }
    }
    return coordinates;
}


double MethodOfRussell::min(double a, double b) {
    return a < b ? a : b;
}

Matrix MethodOfRussell::solve() {
    Matrix initialSolution(costsCoefficients.getRows(), costsCoefficients.getColumns());
    std::vector<int> ignoredColumns;
    std::vector<int> ignoredRows;

    while (true) {
        std::vector<double> rowCoefficients(costsCoefficients.getRows(), 0.0);
        std::vector<double> columnCoefficients(costsCoefficients.getColumns(), 0.0);

        for (int i = 0; i < costsCoefficients.getRows(); ++i) {
            if (std::find(ignoredRows.begin(), ignoredRows.end(), i) == ignoredRows.end()) {
                rowCoefficients[i] = maxRowCoefficient(i, ignoredColumns);
            }
        }
        for (int i = 0; i < costsCoefficients.getColumns(); ++i) {
            if (std::find(ignoredColumns.begin(), ignoredColumns.end(), i) == ignoredColumns.end()) {
                columnCoefficients[i] = maxColumnCoefficient(i, ignoredRows);
            }
        }

        Matrix deltas = deltasMatrix(rowCoefficients, columnCoefficients, ignoredRows, ignoredColumns);
        std::vector<int> minDeltaCoordinates = largestNegativeDeltaCoordinates(deltas);
        int x = minDeltaCoordinates[0];
        int y = minDeltaCoordinates[1];
        double supply = supplyCoefficients.getElement(x, 0);
        double demand = demandCoefficients.getElement(y, 0);

        if (min(supply, demand) == supply) {
            initialSolution.setElement(x, y, supply);
            demandCoefficients.setElement(y, 0, demand - supply);
            ignoredRows.push_back(x);
        } else {
            initialSolution.setElement(x, y, demand);
            supplyCoefficients.setElement(x, 0, supply - demand);
            ignoredColumns.push_back(y);
        }

        if (ignoredRows.size() == static_cast<size_t>(costsCoefficients.getRows()) ||
            ignoredColumns.size() == static_cast<size_t>(costsCoefficients.getColumns())) {
            break;
        }
    }
    return initialSolution;
}
