#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>

using namespace std;

double belso_szorzat(const vector<double> x, const vector<double> y){
    /*  <x, y>
    *   A paraméterként kapott vektorok belső szorzatával tér vissza
    */

    double res = 0;

    for (int i = 0; i < x.size(); i++) {
        res += x[i] * y[i];
    }

    return res;
}

double norm(const vector<double> x){
    /*  ||x||2
    *   A paraméterként kapott vektor 2 normájával tér vissza
    */

    return sqrt(belso_szorzat(x,x));
}

vector<double> mxv(const vector<vector<double> > A, const vector<double> x){
    /*
    *   A paraméterként kapot mátrix és vektor szorzatával tér vissza
    */

    vector<double> res(x.size(), 0);

    for (int i = 0; i < x.size(); i++) {
      for (int j = 0; j < x.size(); j++) {
        res[i] += A[i][j] * x[j];
      }
    }

    return res;
}

vector<double> cxv(const vector<double> v, const double c){
    /*
    *   A paraméterként kapott vektor c-szeresével tér vissza, ahol
    *   c a másik paraméter
    */


    vector<double> res(v.size(), 0);

    for (int i = 0; i < v.size(); i++) {
        res[i] = c * v[i];
    }

    return res;
}

vector<double> vmv(const vector<double> x, const vector<double> y){
    /*
    *   A paraméterként kapott két vektor különbségével tér vissza
    */

    vector<double> res(x.size(), 0);

    for (int i = 0; i < x.size(); i++) {
        res[i] = x[i] - y[i];
    }

    return res;
}


int main(){

    //  Kimenet formázása, minden szám 8 tizedesjegy pontossággal
    cout << fixed << setprecision(8);

    bool singular = false;

    //  Mátrixok száma
    int N = 0;
    cin >> N;

    for (int ni = 0; ni < N; ni++) {

        //  Mátrix mérete
        int n = 0;
        cin >> n;

        //  Mátrix
        vector<vector<double> > A(n, vector<double>(n));

        vector<double> y0(n, 0);
        vector<double> y(n, 0);

        //  Getting matrix from input
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
              cin >> A[i][j];
          }
        }

        vector<vector<double> > ogA(n, vector<double>(n,0));
        ogA = A;

        int m = 0;
        cin >> m;

        for (int mi = 0; mi < m; mi++) {

            double c = 0.0;
            double epsilon = 0.0;
            int maxit = 0;

            singular = false;

            cin >> c;
            cin >> maxit;
            cin >> epsilon;

            for (int i = 0; i < n; i++) {
                cin >> y0[i];
            }

            vector<double> lambda(maxit+1);

            A = ogA;

            //  A-cE
            for (int i = 0; i < n; i++) {
                A[i][i] -= c;
            }

            double max = 0.0;
            int maxRow = 0;

            //  Permutation matrix, simplified to vector form
            vector<int> P(n, 0);
            for (int i = 0; i < n; i++) {
                P[i] = i;
            }

            //  PLU factorization
            for (int k = 0; k < n; k++) {
                max = abs(A[k][k]);
                maxRow = k;

                //  Getting max element of current column
                for (int s = k; s < n; s++)
                    if (max < abs(A[s][k])) {
                        max = abs(A[s][k]);
                        maxRow = s;
                    }

                //  If max element is 0, we can't proceed with the factorization
                if (max < 10e-15) {
                    singular = true;
                    break;
                }

                //  If not zero, we swap the s-th row of A with the k-th row of A
                else {
                    swap(A[maxRow], A[k]);

                    //  and also the s-th element of P with the k-th element of P
                    swap(P[maxRow], P[k]);
                }

                for (int i = k+1; i < n; i++) {
                    A[i][k] = A[i][k] / A[k][k];

                    for (int j = k+1; j < n; j++) {
                        A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
                    }
                }
            }

            if (singular) {
                cout << c << " " << endl;
            }
            else {

                if (norm(y0) < 10e-15) {
                    cout << "kezdovektor" << endl;
                }
                else {

                    //  Normálás
                    double y0norm = norm(y0);

                    for (int i = 0; i < n; i++) {
                        y[i] = y0[i] / y0norm;
                    }

                    //lambda.push_back(belso_szorzat(mxv(ogA, y), y));
                    lambda[0] = belso_szorzat(mxv(ogA, y), y);

                    int k = 0;
                    for (k = 1; k <= maxit; k++) {

                        vector<double> yp(n, 0);

                        for (int i = 0; i < n; i++) {
                            yp[i] = y[P[i]];
                        }

                        //  Solution vector of Lz=yp
                        vector<double> Z(n, 0);

                        //  Lz = yp
                        for (int i = 0; i < n; i++) {
                            double sum = 0;
                            for (int j = 0; j < i; j++) {
                                sum += A[i][j] * Z[j];
                            }

                            Z[i] = yp[i] - sum;
                        }

                        //  Final solution vector
                        vector<double> X(n, 0);

                        // Ux = Z
                        for (int i = n - 1; i >= 0; i--) {
                            double sum = 0;
                            for (int j = i + 1; j < n; j++) {
                                sum += A[i][j] * X[j];
                            }

                            X[i] = (Z[i] - sum) / A[i][i];
                        }

                        double xnorm = norm(X);
                        for (int i = 0; i < n; i++) {
                            y[i] = X[i] / xnorm;
                        }

                        //lambda.push_back(belso_szorzat(mxv(ogA, y), y));
                        lambda[k] = belso_szorzat(mxv(ogA, y), y);

                        // Leállási feltétel
                        if ( abs(lambda[k] - lambda[k-1]) <= epsilon * (1 + abs(lambda[k]))) {
                            break;
                        }
                    }

                    double sikerTesztFeltetel = belso_szorzat(
                        vmv(mxv(ogA, y), cxv(y, lambda[k])),
                        vmv(mxv(ogA, y), cxv(y, lambda[k]))
                    );

                    if (k == maxit+1) {
                        cout << "maxit" << endl;
                    }
                    // Sikerteszt
                    else if (sikerTesztFeltetel <= epsilon) {
                        cout << "siker " << lambda[k] << " ";
                        for (auto yi : y)
                            cout << yi << " ";
                        cout << sikerTesztFeltetel << " " << k-1 << endl;
                    }
                    else {
                        cout << "sikertelen " << lambda[k] << " ";
                        for (auto yi : y)
                            cout << yi << " ";
                        cout << sikerTesztFeltetel << " " << k-1 << endl;
                    }
                }
            }
        }
    }

    return 0;
}
