#include <iostream>
#include <string>
#include "Matrix.h"

class OutputPrinter {
public:
    /**
     * Prints initial solution in a beautiful way.
     * @param solution Initial basic feasible solution - x0
     */
    void printSolution(const Matrix &solution) const {
        for (int i = 0; i < solution.getRows(); i++) {
            for (int j = 0; j < solution.getColumns(); j++) {
                if (solution.getElement(i, j) == 0) {
                    continue;
                }
                std::string s = "X";
                s += std::to_string(i);
                s += std::to_string(j);
                s += " = ";
                s += std::to_string(solution.getElement(i, j));
                std::cout << s << std::endl;
            }
        }
    }

    /**
     * Prints input parameter table (a table constructed using matrix C, vectors S
     * and D) in a beautiful way.
     * @param supply A vector of coefficients of supply - S
     * @param costs A matrix of coefficients of costs - C
     * @param demand A vector of coefficients of demand - D
     */
    void printTable(const Matrix &supply, const Matrix &costs, const Matrix &demand) const {
        Matrix table(costs.getRows() + 1, costs.getColumns() + 1);

        for (int i = 0; i < supply.getRows(); i++) {
            table.setElement(i, costs.getColumns(), supply.getElement(i, 0));
        }
        for (int i = 0; i < costs.getRows(); i++) {
            for (int j = 0; j < costs.getColumns(); j++) {
                table.setElement(i, j, costs.getElement(i, j));
            }
        }
        for (int i = 0; i < demand.getRows(); i++) {
            table.setElement(costs.getRows(), i, demand.getElement(i, 0));
        }

        std::cout << "    COSTS               SUPPLY  " << std::endl;

        for (int i = 0; i < table.getRows(); i++) {
            std::string row = "";
            for (int j = 0; j < table.getColumns(); j++) {
                if (j != costs.getColumns()) {
                    row += std::to_string(table.getElement(i, j)) + " ";
                } else if (i != costs.getRows()) {
                    row += "   |" + std::to_string(table.getElement(i, j));
                }
            }
            std::cout << row << std::endl;

            if (i == costs.getRows() - 1) {
                std::string s = "";
                for (int k = 0; k < costs.getColumns(); k++) {
                    s += "_____";
                }
                std::cout << s << std::endl;
            }
        }
        std::cout << "      DEMAND      " << std::endl;
    }
};
