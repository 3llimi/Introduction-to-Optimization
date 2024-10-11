#include <bits/stdc++.h>


using namespace std;

double EPSILON;

struct Simplex {
    int m, n;
    vector<int> B, N;
    vector<vector<double>> D;

    Simplex(const vector<vector<double>> &A, const vector<double> &b, const vector<double> &c) :
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
    int n, m;
    cin >> n;
    cin >> m;

    vector<double> c(n);
    vector<vector<double>> A(m, vector<double>(n));
    vector<double> b(m);

    for (int i = 0; i < n; i++) cin >> c[i];

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) cin >> A[i][j];

    for (int i = 0; i < m; i++) cin >> b[i];
    double accuracy;
    cin >> accuracy;
    EPSILON = accuracy;
    Simplex simplex(A, b, c);
    vector<double> x;
    double result = simplex.solve(x);

    if (result == numeric_limits<double>::infinity()) {
        cout << "The method is not applicable!" << endl;
    } else {
        cout << "Optimal solution: ";
        for (double xi : x) cout << fixed << setprecision(6) << xi << " ";
        cout << endl;
        cout << "Optimal value: " << fixed << setprecision(6) << result << endl;
    }

    return 0;
}