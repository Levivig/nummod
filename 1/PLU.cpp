#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>

using namespace std;

int main(void) {
    //  Field to hold size of matrix (nxn)
    int n = 0;

    // Field to hold max element of current column
    double max = 0;

    //  Field to hold the row of current max element
    int maxRow = 0;

    // Getting size of matrix from input
    cin >> n;

    //  Matrix singularity
    bool singular = false;

    //  If the size of the matrix is 0, we finish processing (End of input)
    while (n != 0) {

        //  Field to hold matrix
        vector<vector<double> > matrixA(n, vector<double>(n));

        // Getting matrix from input
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
              cin >> matrixA[i][j];
          }
        }

        //  Permutation matrix, simplified to vector form
        vector<int> P(n);
        for (int i = 0; i < n; i++) {
            P[i] = i;
        }

        //  Field to hold number of b vectors
        int m;

        //  Getting number of b vectors
        cin >> m;

        //  Left vector of equation system
        vector<double> B(n);

        //  Solution vector of Ly=Bp
        vector<double> Y(n);

        //  Final solution vector
        vector<double> X(n);

        //  PLU factorization
        for (int k = 0; k < n; k++) {
            max = abs(matrixA[k][k]);
            maxRow = k;

            //  Getting max element of current column
            for (int s = k; s < n; s++)
                if (max < abs(matrixA[s][k])) {
                    max = abs(matrixA[s][k]);
                    maxRow = s;
                }

            //  If max element is 0, we can't proceed with the factorization
            if (max < 10e-15) {
                cout << "szingularis" << endl;
                singular = true;
                break;
            }

            //  If not zero, we swap the s-th row of A with the k-th row of A
            else {
                swap(matrixA[maxRow], matrixA[k]);

                //  and also the s-th element of P with the k-th element of P
                swap(P[maxRow], P[k]);
            }

            for (int i = k+1; i < n; i++) {
                matrixA[i][k] = matrixA[i][k] / matrixA[k][k];

                for (int j = k+1; j < n; j++) {
                    matrixA[i][j] = matrixA[i][j] - (matrixA[i][k] * matrixA[k][j]);
                }
            }
        }

        //  Solving equation system with one b vector at a time
        for (int ii = 0; ii < m; ii++) {

            //  Getting a B vector from input
            for (int i = 0; i < n; i++) {
                cin >> B[i];
            }

            //  New vector to hold the reorderd B vector
            vector<double> Bp(n);

            //  Reordering B according to P
            for (int i = 0; i < n; i++) {
                Bp[i] = B[P[i]];
            }

            //  Ly = b
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += matrixA[i][j] * Y[j];
                }

                Y[i] = Bp[i] - sum;
            }

            // Ux = y
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0;
                for (int j = i + 1; j < n; j++) {
                    sum += matrixA[i][j] * X[j];
                }

                X[i] = (Y[i] - sum) / matrixA[i][i];
            }

            //  If the matrix is NOT singular we print the solution
            if (!singular) {
                for (int i = 0; i < n; i++) {
                    cout << fixed << setprecision(8) << X[i] << " ";
                }

                cout << endl;
            }
        }

        //  Getting size of next matrix
        cin >> n;

        //  Resetting
        singular = false;
    }

    return 0;
}
