#include <bits/stdc++.h>

using namespace std;
double EPSILON;

class ColumnVector {
public:
    vector<double> columnValues;
    unsigned long long columnSize;
    ColumnVector() {
        columnValues = vector<double>();
        columnSize = 0;
    }
    ColumnVector(ColumnVector const &C) {
        columnValues = C.columnValues;
        columnSize = C.columnSize;
    }
    ColumnVector(ColumnVector &C) {
        columnSize = C.columnSize;
        columnValues = C.columnValues;
    }
    ColumnVector(vector<double>& columnValues) {
        this->columnValues = columnValues;
        columnSize = columnValues.size();
    }
    ColumnVector(unsigned long long size) {
        columnSize = size;
        columnValues = vector<double>(size);
    }
    ColumnVector operator+(ColumnVector& b) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] + b.columnValues[i];
        }
        return result;
    }
    ColumnVector operator-(ColumnVector& b) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] - b.columnValues[i];
        }
        return result;
    }
    ColumnVector operator*(double scalar) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] * scalar;
        }
        return result;
    }
    ColumnVector operator/(double scalar) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] / scalar;
        }
        return result;
    }
    double operator*(const ColumnVector& b) const {
        if (columnSize != b.columnSize) {
            throw std::invalid_argument("Column sizes must match for dot product.");
        }
        double result = 0.0;
        for (int i = 0; i < columnSize; i++) {
            result += columnValues[i] * b.columnValues[i];
        }
        return result;
    }
    friend istream& operator>>(istream& stream, ColumnVector& columnVector)
    {
        for (int i = 0; i < columnVector.columnSize; i++) {
            stream >> columnVector.columnValues[i];
        }
        return stream;
    }

    friend ostream& operator<<(ostream& stream, ColumnVector& columnVector)
    {
        for (int i = 0; i < columnVector.columnSize; i++) {
            stream << columnVector.columnValues[i] << endl;
        }
        return stream;
    }
    double norm() {
        double result = 0;
        for (int i = 0; i < columnSize; i++) {
            result += (columnValues[i]*columnValues[i]);
        }
        return sqrt(result);
    }
};


class Matrix {
protected:
    int nbRows, nbColumns;
    vector<ColumnVector> elements;
public:
    int getNbRows() const {
        return nbRows;
    }
    int getNbColumns() const {
        return nbColumns;
    }
    void setNbRows(int newNbRows) {
        nbRows = newNbRows;
    }
    void setNbColumns(int newNbColumns) {
        nbColumns = newNbColumns;
    }
    Matrix(int nbRows, int nbColumns) {
        this->nbRows = nbRows;
        this->nbColumns = nbColumns;
        for (int i = 0; i < nbColumns; i++) {
            elements.emplace_back(nbRows);
            for (int j = 0; j < nbColumns; j++){
                elements[i].columnValues.push_back(0);
            }
        }
    }
    Matrix(const ColumnVector& c){
        this->nbRows = c.columnSize;
        this->nbColumns = c.columnSize;
        for (int i = 0; i < nbColumns; i++) {
            elements.emplace_back(nbRows);
            for (int j = 0; j < nbColumns; j++){
                elements[i].columnValues.push_back(0);
            }
        }
        for (int i = 0; i < c.columnSize; i++) {
            this->setElement(i, i, c.columnValues[i]);
        }
    }
    void setElement(int i, int j, double value) {
        this->elements[j].columnValues[i] = value;
    }
    double& getElement(int i, int j) {
        return this->elements[j].columnValues[i];
    }
    Matrix operator+(Matrix& m1) {

        if (getNbRows() != m1.getNbRows() || getNbColumns() != m1.getNbColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix D(nbRows,nbColumns);
        for (int i = 0; i < nbRows; i++) {
            for (int j = 0; j < nbColumns; j++) {
                D.setElement(i,j, getElement(i,j) + m1.getElement(i,(j)));
            }
        }
        return D;
    }
    Matrix operator-(Matrix& m1) {
        if (getNbRows() != m1.getNbRows() || getNbColumns() != m1.getNbColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix D(nbRows,nbColumns);
        for (int i = 0; i < nbRows; i++) {
            for (int j = 0; j < nbColumns; j++) {
                D.setElement(i,j, getElement(i,j) - m1.getElement(i,(j)));
            }
        }
        return D;
    }
    Matrix operator*(Matrix& m1) {
        if (this->getNbColumns() != m1.getNbRows()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix F(this->getNbRows(),m1.getNbColumns());
        for (int i = 0; i < this->getNbRows(); i++) {
            for (int j = 0; j < m1.getNbColumns(); j++) {
                double result = 0;
                for (int k = 0; k < nbColumns; k++){
                    result += this->getElement(i, k)*m1.getElement(k,j);
                }
                F.setElement(i,j,result);
            }
        }
        return F;
    }

    ColumnVector operator*(ColumnVector& C) {
        if (this->getNbColumns() != C.columnSize) {
            cout << "Error: the dimensional problem occurred\n";
            ColumnVector Empty(0);
            return Empty;
        }
        ColumnVector F(getNbRows());
        for (int i = 0; i < getNbRows(); i++) {
            double result = 0;
            for (int j = 0; j < C.columnSize; j++) {
                result += getElement(i,j)*C.columnValues[j];
            }
            F.columnValues[i] = result;
        }
        return F;
    }

    Matrix& operator=(Matrix m1) {
        this->setNbRows(m1.getNbRows());
        this->setNbColumns(m1.getNbColumns());
        for (int i = 0; i < nbRows; i++) for (int j = 0; j < nbColumns; j++) {
                this->setElement(i,j, m1.getElement(i,j));
            }
        return *this;
    }

    friend istream& operator>>(istream& stream, Matrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream >> matrix.getElement(i,j);
            }
        return stream;
    }

    friend ostream& operator<<(ostream& stream, Matrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) {
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << matrix.getElement(i,j) << " ";
            }
            stream << endl;
        }
        return stream;
    }
    Matrix transpose() {
        Matrix G(nbColumns,nbRows);
        for (int i = 0; i < nbColumns; i++) {
            for (int j = 0; j < nbRows; j++) {
                G.setElement(i,j, getElement(j,i));
            }
        }
        return G;
    }
};

class SquareMatrix: public Matrix {
public:
    explicit SquareMatrix(int size) : Matrix(size, size) {};
    SquareMatrix(Matrix matrix) : Matrix(matrix.getNbRows(), matrix.getNbColumns()) {
        for (int i = 0; i < nbRows; i++)
            for (int j = 0; j < nbColumns; j++)
                setElement(i, j, matrix.getElement(i, j));
    }
    SquareMatrix operator+(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m + *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    SquareMatrix operator-(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m - *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    SquareMatrix operator*(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m * *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    ColumnVector operator*(ColumnVector& C) {
        auto *m = (Matrix*) this;
        ColumnVector *c = &C;
        ColumnVector result = (Matrix)(*m) * *c;
        auto *sResult = (ColumnVector *) &result;
        return *sResult;
    }
    SquareMatrix& operator=(SquareMatrix sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        *m = *m1;
        return *this;
    }
    SquareMatrix transpose() {
        Matrix *m = this;
        Matrix result = m->transpose();
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
};

class IdentityMatrix: public SquareMatrix {
public:
    explicit IdentityMatrix(int size): SquareMatrix(size) {
        for (int i = 0; i < getNbRows(); i++) for (int j = 0; j < getNbColumns(); j++) {
                if (i == j) {
                    this->setElement(i,j,1);
                } else {
                    this->setElement(i,j,0);
                }
            }
    };
};

class EliminationMatrix: public SquareMatrix {
private:

    int iNullified;
    int jNullified;
public:
    EliminationMatrix(SquareMatrix &matrix, int i, int j) : SquareMatrix(matrix.getNbRows()) {
        IdentityMatrix identityMatrix(matrix.getNbRows());
        *this = *((EliminationMatrix *)(SquareMatrix *) (&identityMatrix));
        iNullified = i;
        jNullified = j;
        this->setElement(i,j,-(matrix.getElement(i,j)/matrix.getElement(j,j)));
    }
    int getINullified() const {
        return iNullified;
    }

    int getJNullified() const {
        return jNullified;
    }

};

class PermutationMatrix: public SquareMatrix {
private:
    int firstSwapped;
    int secondSwapped;
    void swapRows(int i1, int i2) {
        for (int i = 0; i < nbRows; i++) {
            swap(getElement(i1,i), getElement(i2,i));
        }
    }
public:
    PermutationMatrix(SquareMatrix& matrix, int i1, int i2) : SquareMatrix(matrix.getNbRows()) {
        IdentityMatrix identityMatrix(matrix.getNbRows());
        *this = *((PermutationMatrix *)(SquareMatrix *) (&identityMatrix));
        firstSwapped = i1;
        secondSwapped = i2;
        swapRows(firstSwapped,secondSwapped);
    }
};

class AugmentedMatrix: public Matrix {
public:
    SquareMatrix m1;
    SquareMatrix m2;
    explicit AugmentedMatrix(SquareMatrix &M1) : Matrix(M1.getNbRows(), M1.getNbColumns()), m1(M1), m2(IdentityMatrix(M1.getNbRows())) {};

    friend ostream& operator<<(ostream& stream, AugmentedMatrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) {
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << fixed << setprecision(2) << matrix.m1.getElement(i,j) << " ";
            }
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << fixed << setprecision(2) << matrix.m2.getElement(i,j) << " ";
            }
            stream << endl;
        }
        return stream;
    }
};

class FindInverseMatrix {
private:
    int step;

    static int findPermutation(AugmentedMatrix& matrix, int column) {
        int pivot = column;
        for (int j = column+1; j < matrix.getNbColumns(); j++) {
            if (abs(matrix.m1.getElement(j,column)) > abs(matrix.m1.getElement(pivot,column))) {
                pivot = j;
            }
        }
        return pivot;
    }
    void makePermutation(AugmentedMatrix& matrix, int row1, int row2) {
        if (row1 != row2) {
            PermutationMatrix P(matrix.m1, row1, row2);
            matrix.m1 = P*matrix.m1;
            matrix.m2 = P*matrix.m2;
            step++;
        }
    }
    void eliminatePosition(AugmentedMatrix& matrix, int i, int j) {
        EliminationMatrix E(matrix.m1,i,j);
        matrix.m1 = E*matrix.m1;
        matrix.m2 = E*matrix.m2;
        step++;
    }

    void directWay(AugmentedMatrix& augmentedMatrix) {
        for (int i = 0; i < augmentedMatrix.getNbColumns()-1; i++) {
            int pivot = findPermutation(augmentedMatrix,i);
            makePermutation(augmentedMatrix,i,pivot);
            for (int j = i+1; j < augmentedMatrix.getNbColumns(); j++) {
                if (augmentedMatrix.m1.getElement(j,i) != 0) {
                    eliminatePosition(augmentedMatrix,j,i);
                }
            }
        }
    }

    void wayBack(AugmentedMatrix& augmentedMatrix) {
        for (int i = augmentedMatrix.getNbColumns()-1; i > 0; i--) {
            for (int j = i-1; j >= 0; j--) {
                if (augmentedMatrix.m1.getElement(j,i) != 0) {
                    eliminatePosition(augmentedMatrix,j,i);
                }
            }
        }
    }

    static void normalizeDiagonal(AugmentedMatrix& augmentedMatrix) {
        for (int i = 0; i < augmentedMatrix.getNbRows(); i++) {
            double pivot = augmentedMatrix.m1.getElement(i,i);
            augmentedMatrix.m1.setElement(i,i,augmentedMatrix.m1.getElement(i,i)/pivot);
            for (int j = 0; j < augmentedMatrix.getNbColumns(); j++) {
                augmentedMatrix.m2.setElement(i,j,augmentedMatrix.m2.getElement(i,j)/pivot);
            }
        }
    }
public:
    SquareMatrix findInverse(SquareMatrix matrix) {
        step = 0;
        AugmentedMatrix augmentedMatrix(matrix);
        step++;
        directWay(augmentedMatrix);
        wayBack(augmentedMatrix);
        normalizeDiagonal(augmentedMatrix);
        matrix = augmentedMatrix.m2;
        return matrix;
    }
};


void InteriorPointsolve(ColumnVector c, Matrix A, ColumnVector rightHand, ColumnVector initial, double alpha, int nb_variables, int nb_equations, int itr) {
    Matrix D(initial);
    Matrix A_ = A * D;
    ColumnVector c_ = D * c;
    Matrix A_T = A_.transpose();
    SquareMatrix A_A_T(A_ * A_T);
    FindInverseMatrix inverseMatrix;
    Matrix A_A_T_1 = inverseMatrix.findInverse(A_A_T);
    IdentityMatrix I(nb_variables);
    Matrix P = A_A_T_1 * A_;
    P = A_T * P;
    P = Matrix(I) - P;
    ColumnVector c_p = P * c_;
    double mini = 0;
    for (int i = 0; i < nb_variables; i++) {
        mini = min(c_p.columnValues[i], mini);
    }
    double v = abs(mini);
    ColumnVector ones(nb_variables);
    for (int i = 0; i < nb_variables; i++) {
        ones.columnValues[i] = 1;
    }
    ColumnVector x_ = c_p * (alpha / v);
    x_ = x_ + ones;
    ColumnVector x = D * x_;
    bool found = true;
    for (int i = 0; i < nb_variables; i++) {
        if (abs(x.columnValues[i] - initial.columnValues[i]) > 1e-4)
            found = false;
    }
    if (found) {
        cout << "Optimal Solution (alpha = " << alpha << "):\n" << x;
        cout << "Objective Function Value: " << (c * x) << endl;
        return;
    } else {
        InteriorPointsolve(c, A, rightHand, x, alpha, nb_variables, nb_equations, itr + 1);
    }
}

bool checkValid(Matrix A, ColumnVector rightHand, ColumnVector initial, int nb_variables, int nb_equations) {
    for (int i = 0; i < nb_variables; i++) {
        if (initial.columnValues[i] <= 0) return false;
    }
    for (int i = 0; i < nb_equations; i++) {
        double sum = 0;
        for (int j = 0; j < nb_variables; j++) {
            sum += A.getElement(i, j) * initial.columnValues[j];
        }
        if (sum != rightHand.columnValues[i]) return false;
    }
    return true;
}



struct Simplex {
    int m, n;
    vector<int> B, N;
    vector<vector<double> > D;

    Simplex(const vector<vector<double> > &A, const vector<double> &b, const vector<double> &c) :
        m(b.size()), n(c.size()), B(m), N(n + 1), D(m + 2, vector<double>(n + 2)) {
        for (int i = 0; i < m; i++) B[i] = n + i;
        for (int j = 0; j < n; j++) N[j] = j;
        N[n] = -1;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) D[i][j] = A[i][j];
            D[i][n] = -1;
            D[i][n + 1] = b[i];
        }
        for (int j = 0; j < n; j++) D[m][j] = -c[j];
        D[m + 1][n] = 1;
    }

    void pivot(int r, int s) {
        double inv = 1.0 / D[r][s];
        for (int i = 0; i < m + 2; i++) if (i != r)
                for (int j = 0; j < n + 2; j++) if (j != s)
                        D[i][j] -= D[r][j] * D[i][s] * inv;
        for (int j = 0; j < n + 2; j++) if (j != s) D[r][j] *= inv;
        for (int i = 0; i < m + 2; i++) if (i != r) D[i][s] *= -inv;
        D[r][s] = inv;
        swap(B[r], N[s]);
    }

    bool simplex(int phase) {
        int x = phase == 1 ? m + 1 : m;
        while (true) {
            int s = -1;
            for (int j = 0; j <= n; j++) {
                if (N[j] == -1) continue;
                if (s == -1 || D[x][j] < D[x][s] || (D[x][j] == D[x][s] && N[j] < N[s])) s = j;
            }
            if (D[x][s] > -EPSILON) return true;
            int r = -1;
            for (int i = 0; i < m; i++) {
                if (D[i][s] < EPSILON) continue;
                if (r == -1 || D[i][n + 1] / D[i][s] < D[r][n + 1] / D[r][s] ||
                    (D[i][n + 1] / D[i][s] == D[r][n + 1] / D[r][s] && B[i] < B[r])) r = i;
            }
            if (r == -1) return false;
            pivot(r, s);
        }
    }

    double solve(vector<double> &x) {
        int r = 0;
        for (int i = 1; i < m; i++) if (D[i][n + 1] < D[r][n + 1]) r = i;
        if (D[r][n + 1] < -EPSILON) {
            pivot(r, n);
            if (!simplex(1) || D[m + 1][n + 1] < -EPSILON) return -numeric_limits<double>::infinity();
            for (int i = 0; i < m; i++) if (B[i] == -1) {
                    int s = -1;
                    for (int j = 0; j <= n; j++) if (s == -1 || D[i][j] < D[i][s]) s = j;
                    pivot(i, s);
                }
        }
        if (!simplex(2)) return numeric_limits<double>::infinity();
        x = vector<double>(n);
        for (int i = 0; i < m; i++) if (B[i] < n) x[B[i]] = D[i][n + 1];
        return D[m][n + 1];
    }
};

int main() {
    int nb_variables, nb_equations;
    double epsilon;

    cout << "Enter number of variables and constraints: ";
    cin >> nb_variables >> nb_equations;

    ColumnVector c(nb_variables);
    cout << "Enter coefficients of objective function: ";
    cin >> c;

    Matrix A(nb_equations, nb_variables);
    cout << "Enter constraint coefficient matrix:\n";
    cin >> A;

    ColumnVector b(nb_equations);
    cout << "Enter right-hand side vector b: ";
    cin >> b;

    ColumnVector initial(nb_variables);
    cout << "Enter initial point vector: ";
    cin >> initial;

    cout << "Enter approximation accuracy ε: ";
    cin >> epsilon;


    cout << "Solution with alpha = 0.5:\n";
    InteriorPointsolve(c, A, b, initial, 0.5, nb_variables, nb_equations, 1);

    cout << "Solution with alpha = 0.9:\n";
    InteriorPointsolve(c, A, b, initial, 0.9, nb_variables, nb_equations, 1);
    vector<vector<double> > B(nb_variables, vector<double>(nb_equations));
    vector<double> b1(nb_equations);
    vector<double> c1(nb_variables);

    Simplex simplex(B, b1, c1);
    vector<double> x;
    double result = simplex.solve(x);

    if (result == numeric_limits<double>::infinity()) {
        cout << "The method is not applicable!" << endl;
    } else if (result == -numeric_limits<double>::infinity()) {
        cout << "The method is not applicable!" << endl;
    } else {
        cout << "Optimal solution: ";
        for (double xi : x) cout << fixed << setprecision(6) << xi << " ";
        cout << endl;
        cout << "Optimal value: " << fixed << setprecision(6) << result << endl;
    }

    return 0;
}
