#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
using namespace std;
using namespace arma;


double function_f(double x);
double deriv_u(double x);
double RelativeError(double v, double u);

int main(){
    int n = 10;
    double *a, *b, *c, *b_tilde, *v;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    b_tilde = new double[n];
    v = new double[n];


    double h = 1.0/(double (1 + n));

    //x[0] = h;
//Filling in the diagonals and b_tilde as we have in 1b
for(int i =0; i<n;i++){
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
    b_tilde[i] = function_f(i*h)*h*h;
}

//Forward substitution

for(int i =1; i<n;i++){
    b[i] = b[i]-a[i]*c[i-1]/b[i-1];
    b_tilde[i] = b_tilde[i] - a[i]*b_tilde[i-1]/b[i-1];
}


//Backward substitution
v[n]
for(int i =n-2; i>=0 ; i--){
    cout <<i<<endl;
}







    return 0;
}



double function_f(double x){
    double f = 100*exp(-10*x);
    return f;
}


double deriv_u(double x){
    double u = 1- (1-exp(-10.0))*x - exp(-10.0*x);
    return u;
}

double RelativeError(double v, double u){
    double epsilon = (abs((v-u)/u));
    return epsilon;
}


















/*
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
    double h = 1.0/(double (1 + n));

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

    // forward substitution
    for(int i = 1; i<n; i++){
        new_d[i] = d[i] - (e[i-1]*e[i-1])/(double)new_d[i-1];

        new_b[i] = b[i] - (new_b[i-1]*e[i-1])/(double)new_d[i-1];
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
*/
