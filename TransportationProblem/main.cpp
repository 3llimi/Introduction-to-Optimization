#include <iostream>
#include "InputReader.cpp"
#include "InputValidator.cpp"
#include "OutputPrinter.cpp"
#include "NorthWestCornerMethod.cpp"
#include "MethodOfVogel.cpp"
#include "MethodOfRussell.cpp"
#include "Matrix.cpp"
#include "Matrix.h"

int main() {
    InputReader reader;
    std::cout << "Enter the number of sources:" << std::endl;
    int numOfSupply = reader.readInt();
    std::cout << "Enter the number of destinations:" << std::endl;
    int numOfDemand = reader.readInt();

    Matrix supplyCoefficients(numOfSupply, 1);
    std::cout << "Enter a vector of coefficients of supply - S:" << std::endl;
    reader.readMatrixNx1(supplyCoefficients);

    Matrix costsCoefficients(numOfSupply, numOfDemand);
    std::cout << "Enter a matrix of coefficients of costs - C:" << std::endl;
    reader.readMatrixNxM(costsCoefficients);

    Matrix demandCoefficients(numOfDemand, 1);
    std::cout << "Enter a vector of coefficients of demand - D:" << std::endl;
    reader.readMatrixNx1(demandCoefficients);

    InputValidator validator;
    bool isApplicable = validator.isMethodApplicable(supplyCoefficients, costsCoefficients, demandCoefficients);
    bool isBalanced = validator.isProblemBalanced(supplyCoefficients, demandCoefficients);

    if (!isApplicable) {
        std::cout << "The method is not applicable!" << std::endl;
        return 1;
    }
    if (!isBalanced) {
        std::cout << "The problem is not balanced!" << std::endl;
        return 1;
    }

    NorthWestCornerMethod method1(supplyCoefficients, costsCoefficients, demandCoefficients);
    MethodOfVogel method2(supplyCoefficients, costsCoefficients, demandCoefficients);
    MethodOfRussell method3(supplyCoefficients, costsCoefficients, demandCoefficients);

    Matrix initialSolution1 = method1.solve();
    Matrix initialSolution2 = method2.solve();
    Matrix initialSolution3 = method3.solve();

    OutputPrinter printer;
    printer.printTable(supplyCoefficients, costsCoefficients, demandCoefficients);

    std::cout << "Initial basic feasible solution - x0 using North-West corner method:" << std::endl;
    printer.printSolution(initialSolution1);
    std::cout << "Initial basic feasible solution - x0 using Vogel’s approximation method:" << std::endl;
    printer.printSolution(initialSolution2);
    std::cout << "Initial basic feasible solution - x0 using Russell’s approximation method:" << std::endl;
    printer.printSolution(initialSolution3);

    return 0;
}
