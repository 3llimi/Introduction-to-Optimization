#include <algorithm>
#include "Matrix.h"

class NorthWestCornerMethod {
private:
    Matrix supplyCoefficients;
    Matrix costsCoefficients;
    Matrix demandCoefficients;

    /**
     * Returns the minimum of two numbers.
     *
     * @param a 1st double number
     * @param b 2nd double number
     * @return minimum from numbers a and b
     */
    static double min(double a, double b) {
        return std::min(a, b);
    }

public:
    /**
     * Constructs a NorthWestCornerMethod solver with specified supply, demand,
     * and costs matrices.
     * The sum of supply coefficients is assumed to be equal to the sum of demand
     * coefficients. If this condition is not satisfied, the behavior of the solver
     * is undefined.
     *
     * @param coefficientsOfSupply Supply coefficients
     * @param coefficientsOfCosts  Cost coefficients
     * @param coefficientsOfDemand Demand coefficients
     */
    NorthWestCornerMethod(const Matrix& coefficientsOfSupply,
                          const Matrix& coefficientsOfCosts,
                          const Matrix& coefficientsOfDemand)
            : supplyCoefficients(coefficientsOfSupply),
              costsCoefficients(coefficientsOfCosts),
              demandCoefficients(coefficientsOfDemand) {}

    /**
     * Solves using the North-West corner method to find the initial solution X0.
     * This method does not change the solver's state, so repeated invocations yield
     * consistent results.
     *
     * @return initial solution obtained using the North-West corner method
     */
    Matrix solve() {
        Matrix solution(costsCoefficients.getRows(), costsCoefficients.getColumns());
        Matrix demand(demandCoefficients);
        Matrix supply(supplyCoefficients);

        int i = 0;
        int j = 0;
        while (i < solution.getRows() && j < solution.getColumns()) {
            double allocation = min(demand.getElement(j, 0), supply.getElement(i, 0));
            demand.setElement(j, 0, demand.getElement(j, 0) - allocation);
            supply.setElement(i, 0, supply.getElement(i, 0) - allocation);

            solution.setElement(i, j, allocation);

            if (demand.getElement(j, 0) == 0) {
                j += 1;
            }

            if (supply.getElement(i, 0) == 0) {
                i += 1;
            }
        }

        return solution;
    }
};
