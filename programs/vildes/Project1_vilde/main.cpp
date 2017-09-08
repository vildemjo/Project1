#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
    ofstream outfile;
    outfile.open("matrix_values.txt");

    int n = 5;
    double h = 1.0/(1 + n);
    double x[n];
    x[0] = 0;
    b[0] = 0;
    double b[n];
    int Matrix[n][n];

    double a1[n];     // this is the superdiagonal
    double a2[n];     // this is the diagonal
    double a3[n];     // this is the subdiagonal

    for(int i = 0; i<n; i++){
        a1[i] = -1;
        a2[i] = 2;
        a3[i] = -1;

        for(int j=0; j<n ;j++){
        Matrix[i][j] = 0;
        }

        Matrix[i][i] = a2[i];
        Matrix[i][i+1] = a1[i];
        Matrix[i][i-1] = a3[i];

        x[i+1] = x[i] + h;
        b[i+1] = close_function(x[i+1])

    }

    for(j)
       new_a[j][k] = a[j][k] - (a[j][m]*a[m][k])/ (double) a[m][m];
       // Finding coeffientist for the tilde_b
       new_b[j] = b[j] - (a[j][m]*b[m])/ (double) a[m][m];


       a[j][k] = new_a[j][k];

       b[j] = new_b[j];

    outfile.close();
    return 0;
}
