#include <bits/stdc++.h>

using namespace std;


// Class for implementing the Golden Section Search Method for finding the minimum of a function
class GoldenSectionMethod {
public:
    GoldenSectionMethod(double tolerance);  // Constructor to initialize the method with a tolerance value
    void findMinimum(double a, double b); // Finds the minimum of the function within the interval [a, b]
    double getXMin() const; // Returns the x-coordinate of the minimum point
    double getFMin() const; // Returns the function value at the minimum point

private:
    double tolerance; // Tolerance value for stopping
    double x_min; // Stores the x-coordinate of the minimum point
    double f_min; // Stores the function value at the minimum point

    double f(double x); // The function to be minimized
};

// Constructor initializes tolerance and default values for x_min and f_min
GoldenSectionMethod::GoldenSectionMethod(double tolerance)
    : tolerance(tolerance), x_min(0), f_min(0) {}

// Function to be minimized
double GoldenSectionMethod::f(double x) {
    return (x - 2) * (x - 2) + 3;  // Minimum at x = 2
}

// Implements the Golden Section Search algorithm
void GoldenSectionMethod::findMinimum(double a, double b) {
    const double phi = (1 + sqrt(5)) / 2; // Golden ratio
    double c = b - (b - a) / phi; // First intermediate point
    double d = a + (b - a) / phi; // Second intermediate point

    // Continue the search until the interval width is within the tolerance
    while ((b - a) > tolerance) {
        if (f(c) < f(d)) {
            b = d; // Narrow the interval to [a, d]
        } else {
            a = c; // Narrow the interval to [c, b]
        }
        // Recalculate intermediate points
        c = b - (b - a) / phi;
        d = a + (b - a) / phi;
    }

    x_min = (a + b) / 2; // Approximate x_min as the midpoint of the final interval
    f_min = f(x_min); // Evaluate the function at x_min
}

// Getter for x_min
double GoldenSectionMethod::getXMin() const {
    return x_min;
}

// Getter for f_min
double GoldenSectionMethod::getFMin() const {
    return f_min;
}

// Class for reading input values
class InputReader {
public:
    InputReader();
    double readDouble();

private:
    istream &input;
};


InputReader::InputReader() : input(cin) {}


double InputReader::readDouble() {
    double value;
    input >> value;
    return value;
}


int main() {
    InputReader reader; // Object for handling input

    double a, b, epsilon; // Variables for interval and tolerance

    while (true) { // Loop until valid input is provided
        try {
            cout << "Enter interval [a, b]:\n";
            cout << "a: ";
            a = reader.readDouble(); // Read value for a
            cout << "b: ";
            b = reader.readDouble(); // Read value for b

            cout << "Enter tolerance (ϵ): ";
            epsilon = reader.readDouble(); // Read value for tolerance

            // Check if the interval is valid
            if (a >= b) {
                throw invalid_argument("Error: a must be less than b.");
            }

            // Create an instance of GoldenSectionMethod with the given tolerance
            GoldenSectionMethod goldenSection(epsilon);
            goldenSection.findMinimum(a, b); // Perform the search

            cout << fixed << setprecision(6);
            cout << "Approximate x_min: " << goldenSection.getXMin() << endl;
            cout << "f(x_min): " << goldenSection.getFMin() << endl;
            break;
        } catch (const invalid_argument &e) {
            // Handle invalid input and prompt the user to try again
            cerr << e.what() << "\nPlease enter a valid interval.\n";
        }
    }

    return 0;
}
