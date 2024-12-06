#include <bits/stdc++.h>

using namespace std;
// Class to handle Gradient Ascent
class GradientAscent {
public:
    GradientAscent(double alpha, int iterations);
    void maximize(double x0);
    double getXMax() const;
    double getFMax() const;

private:
    double alpha;  // Learning rate
    int iterations; // Number of iterations
    double x_max;   // x value that maximizes f(x)
    double f_max;   // Maximum value of f(x)

    double f(double x);        // Function f(x)
    double fPrime(double x);   // Derivative f'(x)
};

// Constructor for GradientAscent
GradientAscent::GradientAscent(double alpha, int iterations)
    : alpha(alpha), iterations(iterations), x_max(0), f_max(0) {}

// Function f(x) = -x^2 + 4x + 1
double GradientAscent::f(double x) {
    return -x * x + 4 * x + 1;
}

// Derivative f'(x) = -2x + 4
double GradientAscent::fPrime(double x) {
    return -2 * x + 4;
}

// Perform Gradient Ascent
void GradientAscent::maximize(double x0) {
    double x = x0;
    for (int i = 0; i < iterations; ++i) {
        double gradient = fPrime(x);
        x += alpha * gradient; // Update rule
    }
    x_max = x;
    f_max = f(x_max);
}

// Get the x value that maximizes f(x)
double GradientAscent::getXMax() const {
    return x_max;
}

// Get the maximum value of f(x)
double GradientAscent::getFMax() const {
    return f_max;
}

// InputReader class for reading inputs
class InputReader {
public:
    InputReader();
    double readDouble();
    int readInt();

private:
    istream &input;
};

// Constructor for InputReader
InputReader::InputReader() : input(cin) {}

// Read a double value
double InputReader::readDouble() {
    double value;
    input >> value;
    return value;
}

// Read an integer value
int InputReader::readInt() {
    int value;
    input >> value;
    return value;
}

// Main function
int main() {
    InputReader reader;

    cout << "Enter initial guess (x0): ";
    double x0 = reader.readDouble();

    cout << "Enter learning rate (alpha): ";
    double alpha = reader.readDouble();

    cout << "Enter number of iterations (N): ";
    int iterations = reader.readInt();

    GradientAscent ascent(alpha, iterations);
    ascent.maximize(x0);

    cout << fixed << setprecision(4);
    cout << "Approximate x_max: " << ascent.getXMax() << endl;
    cout << "f(x_max): " << ascent.getFMax() << endl;

    return 0;
}
