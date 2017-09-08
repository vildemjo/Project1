#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double source_term(double x);
double closed_form_solution(double x);
double RelativeError(double v, double u);

int main()
{
    ofstream outfile;
    outfile.open("matrix_values.txt");

    outfile << fixed;
    outfile << setprecision(4);

    // The variables
    int n = 20;
    double h = 1.0/(double)(1 + n);

    //The vectors // new means after forward substitution
    double x[n], v[n], b[n], new_b[n];

    // The matrix values // new means after backward substitution
    double e[n];     // this is the superdiagonal and the subdiagonal - symmetrical
    double d[n], new_d[n];     // this is the diagonal

    double Error[n];

     x[0] = h;

    for(int i = 0; i<n; i++){
        e[i] = -1.0;
        d[i] = 2.0;

        b[i] = source_term(x[i])*h*h;
        x[i+1] = x[i] + h;
    }

    new_b[0] = b[0];
    new_d[0] = d[0];
    v[0] = 0.0;

    // printing to file - title and initial values
    outfile << "d" << "     " << "d_new" << "   " << "b" << "    " << "b_new" << endl;
    outfile << d[0] << " " << new_d[0] << " " << b[0] << " " << new_b[0] << endl;

    // forward substitution
    for(int i = 1; i<n; i++){
        new_d[i] = d[i] - (e[i-1]*e[i-1])/(double)new_d[i-1];

        new_b[i] = b[i] - (new_b[i-1]*e[i-1])/(double)new_d[i-1];

        outfile << d[i] << " " << new_d[i]<< " " << b[i] << " " <<new_b[i] << endl;
    }

    // backward substitution
    v[n-1] = new_b[n-1]/(double) new_d[n-1];

    for(int j = 2; j<n ; j++){
        int i = n-j;
        v[i] = (new_b[i] -e[i]*v[i+1])/new_d[i];
    }

    double MaxError = 0.0;

    for(int i = 1; i<n; i++){
        Error[i] = RelativeError(v[i], closed_form_solution(x[i]));
        if(MaxError<Error[i]){
            MaxError = Error[i];
        }
        outfile << Error[i] << ", ";
    }


    for(int i = 1; i<n; i++)
        cout << v[i] << "  " << closed_form_solution(x[i]) << endl;

    outfile << " " << endl;

    outfile << "MaxError = " << MaxError;

    outfile.close();

    cout << closed_form_solution(x[1]);

    return 0;
}

double source_term(double x){
    double f = 100*exp(-10*x);
    return f;
}

double closed_form_solution(double x){
    double u = 1- (1-exp(-10.0))*x - exp(-10.0*x);
    return u;
}

double RelativeError(double v, double u){
    double epsilon = (abs((v-u)/u));
    return epsilon;
}
