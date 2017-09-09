#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

double source_term(double x);
double closed_form_solution(double x);
double RelativeError(double v, double u);
double* Calculated_solution(int n, double* v, double* x, double h);
void print_to_file(int n, double *v, double *x);

int main(){

    for(int i= 1; i<5; i++){
        int n = pow(10,i);//10^i;
        double h = 1.0/((double)(1 + n));
        double x[n];
        double v[n];
        x[0] = h;

        for(int i = 1; i<n; i++){
            x[i] = x[i-1] + h;
        }

        Calculated_solution(n, v, x, h);
        print_to_file(n, v, x);

        double MaxError = 0.0;
        double Error;

        for(int i = 1; i<n; i++){
            Error = RelativeError(v[i], closed_form_solution(x[i]));
            if(MaxError<Error){
                MaxError = Error;
            }
        }
        cout << "n=" << n <<" MaxError=" << MaxError << endl;
    }


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
    double epsilon = log10(abs((v-u)/u));
    return epsilon;
}

double* Calculated_solution(int n, double *v, double *x, double h){

    // The variables

    //The vectors // new means after forward substitution
    double b[n], new_b[n];

    // The matrix values // new means after backward substitution
    double* e = new double[n];     // this is the superdiagonal and the subdiagonal - symmetrical
    double d[n], new_d[n];     // this is the diagonal

    for(int i = 0; i<n; i++){
        e[i] = -1.0;
        d[i] = 2.0;

        b[i] = source_term(x[i])*h*h;
    }

    new_b[0] = b[0];
    new_d[0] = d[0];

    // forward substitution
    for(int i = 1; i<n; i++){
        new_d[i] = d[i] - (e[i]*e[i-1])/(double)new_d[i-1];

        new_b[i] = b[i] - (new_b[i-1]*e[i-1])/(double)new_d[i-1];

    }

    // backward substitution
    v[n-1] = new_b[n-1]/(double) new_d[n-1];

    for(int i = n-2; i>=0 ; i--){
        //        int i = n-j;
        v[i] = (new_b[i] -e[i]*v[i+1])/new_d[i];
    }

    return v;

}

void print_to_file(int n, double* v, double* x){

    ofstream outfile;

    string str = to_string(n);
    string filename = string("result_n_is") + str + string(".txt");

    outfile.open(filename);

    outfile << fixed;
    outfile << setprecision(4);

    outfile << "Calculated solution:" << " , " << "Analytical solution:" << endl;

    for(int i = 0; i<n;i++){
        outfile << v[i] << "            " << closed_form_solution(x[i]) << endl;
    }

    outfile.close();
}
