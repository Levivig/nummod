#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


using namespace std;

int fakt(int n) {
    if(n == 0) return 1;

    int eredm = 1;

    for (int i = 1; i <= n; i++) {
        eredm *= i;
    }

    return eredm;
}

int get_index(vector<double> vect, double itemToFind) {
    int index = 0;
    for (int i = 0; i < vect.size(); i++) {
        if( abs(vect[i] - itemToFind) < 10e-15)
            index = i;
    }
    return index;
}

int main() {
    bool debug = true;

    cout << fixed << setprecision(0);

    //  Feladatok száma
    int N;
    cin >> N;


    for (int in = 0; in < N; in++) {
        //  adatok száma
        int n;
        cin >> n;

        //  Illesztési feltételek száma
        int M;
        cin >> M;

        //  pontok vektora
        vector<double> x(n + 1);

        //  függvényértékek mátrixa
        vector<vector<double> > f(n+1, vector<double>());

        vector<int> m(n+1);

        // Getting x from input
        for (int i = 0; i < n + 1; i++) {
              cin >> x[i];
              cin >> m[i];

              f[i].resize(m[i]);
              for (int j = 0; j < m[i]; j++) {
                  cin >> f[i][j];
              }
        }

        vector<double> xNew;

        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m[i]; j++) {
                xNew.push_back(x[i]);
            }
        }

        vector<double> f_new;

        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m[i]; j++) {
                f_new.push_back(f[i][0]);
            }
        }

        int ym;
        cin >> ym;

        vector<double> y(ym);
        for (int i = 0; i < ym; i++) {
            cin >> y[i];
        }

        if (debug) {
            cout << "X new: " << endl;
            for (int i = 0; i < M; i++) {
                cout << xNew[i] << " ";
            }
            cout << endl;

            cout << "F new: " << endl;
            for (int i = 0; i < M; i++) {
                    cout << f_new[i] << " ";
                }
                cout << endl;

            cout << "m: " << endl;
            for (int i = 0; i < n + 1; i++) {
                cout << m[i] << " ";
            }
            cout << endl;


            cout << "y: " << endl;
            for (int i = 0; i < ym; i++) {
                cout << y[i] << " ";
            }
            cout << endl;

        }

        vector<vector<double> > diffTabla(M, vector<double>(M));

        for (int i = 0; i < M; i++) {
            diffTabla[i][0] = f_new[i];
        }


        int columnNo = 0;
        int rowNo = 0;

        for (int k = 1; k < M; k++) {
            columnNo++;
            for (int i = k; i < M; i++) {
                if ( abs(xNew[i] - xNew[i-k]) < 10e-15) {
                    rowNo = get_index(x, xNew[i]);

                    diffTabla[i][k] = f[rowNo][columnNo] / fakt(columnNo);
                    if(rowNo == m[columnNo] - 1) rowNo=1;
                }
                else{
                    diffTabla[i][k] = (diffTabla[i][k-1] - diffTabla[i-1][k-1]) / (xNew[i] - xNew[i-k]);
                }
            }
        }

        if (debug) {
            cout << "diffTabla" << endl;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j <= i; j++) {
                    cout << setw(10) << diffTabla[i][j] << " ";
                }
                cout << endl;
            }
        }
        cout << "------------------------------" << endl;


        vector<double> b(M);

        for (int i = 0; i < M; i++) {
            b[i] = diffTabla[i][i];
        }

        for (int i = 0; i < M; i++) {
            cout << b[i] << " ";
        }
        cout << endl;

        for (int i = 0; i < ym; i++) {
            double c = b[M - 1];

            //  Horner
            for (int k = 1; k < M; k++) {
                c = c*(y[i] - xNew[M - 1 - k]) + b[M - 1 - k];
            }

            cout << c << " ";
        }
        cout << endl;
    }

    return 0;
}
