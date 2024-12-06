#include <bits/stdc++.h>

using namespace std;


// Class for implementing the Bisection Method to find the root of a function
class BisectionMethod {
public:
    BisectionMethod(double tolerance); // Constructor to initialize the tolerance
    double findRoot(double a, double b); // Method to find the root within an interval [a, b]
    double getRoot() const; // Getter for the computed root

private:
    double tolerance; // Tolerance value for stopping criterion
    double root; // Variable to store the computed root

    double f(double x); // Function for which the root is to be found
};

// Constructor initializes tolerance and sets root to 0
BisectionMethod::BisectionMethod(double tolerance) : tolerance(tolerance), root(0) {}

// Define the function for which the root will be calculated
double BisectionMethod::f(double x) {
    return x * x * x - 6 * x * x + 11 * x - 6; // A cubic polynomial
}

// Method to find the root using the Bisection Metho
double BisectionMethod::findRoot(double a, double b) {
    // Check if the endpoints of the interval are roots
    if (f(a) == 0) {
        root = a;
        return root;
    }
    if (f(b) == 0) {
        root = b;
        return root;
    }

    // Ensure the function values at a and b have opposite signs
    if (f(a) * f(b) > 0) {
        cerr << "Error: f(a) = " << f(a) << ", f(b) = " << f(b) << ". "
                  << "They must have opposite signs. Please try again with a different interval.\n";
        throw invalid_argument("f(a) and f(b) must have opposite signs.");
    }

    double c; // Midpoint of the interval
    while ((b - a) / 2.0 > tolerance) { // Continue until interval size is within tolerance
        c = (a + b) / 2.0; // Compute the midpoint
        if (abs(f(c)) < tolerance) { // If the function value at c is close to 0, stop
            root = c;
            return root;
        }
        // Determine which subinterval contains the root
        if (f(a) * f(c) < 0) {
            b = c; // Root is in [a, c]
        } else {
            a = c; // Root is in [c, b]
        }
    }
    root = (a + b) / 2.0; // Set the root as the midpoint of the final interval
    return root;
}

// Getter method to return the computed root
double BisectionMethod::getRoot() const {
    return root;
}


// Class for reading user input
class InputReader {
public:
    InputReader();
    double readDouble();

private:
    istream &input;
};

// Constructor initializes the input stream to standard input
InputReader::InputReader() : input(cin) {}

double InputReader::readDouble() {
    double value;
    input >> value;
    return value;
}


int main() {
    InputReader reader; // Create an InputReader object

    double a, b, epsilon;  // Variables for the interval and tolerance

    while (true) { // Loop until valid input is provided
        try {
            cout << "Enter interval [a, b]:\n";
            cout << "a: ";
            a = reader.readDouble(); // Read the left endpoint of the interval
            cout << "b: ";
            b = reader.readDouble(); // Read the right endpoint of the interval

            cout << "Enter tolerance (ϵ): ";
            epsilon = reader.readDouble(); // Read the tolerance

            BisectionMethod bisection(epsilon); // Create a BisectionMethod object with the given tolerance
            double root = bisection.findRoot(a, b); // Find the root within the interval

            cout << fixed << setprecision(6);
            cout << "Approximate root: " << root << endl;
            break;
        } catch (const invalid_argument &e) { // Handle invalid intervals
            cerr << e.what() << "\nPlease enter a valid interval.\n";
        }
    }

    return 0;
}
