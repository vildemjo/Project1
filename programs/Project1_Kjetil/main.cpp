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

ofstream outfile;

int main(){

for(int n = 10;n<11; n=n*10){
    //int n = 10;
    double *a, *b, *c, *b_tilde, *v;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    b_tilde = new double[n];
    v = new double[n];
    double h = 1.0/(double (n+1));
    string outfilename = "../Plotting_kjetil/1b_n_";

    outfilename += to_string(n);
    outfilename += ".txt";

//Filling in the diagonals and b_tilde as we have in 1b
for(int i =0; i<n;i++){
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
    b_tilde[i] = function_f((i+1)*h)*h*h;
}



//Forward substitution

for(int i =1; i<n;i++){
    b[i] = b[i]-a[i-1]*c[i-1]/b[i-1];
    b_tilde[i] = b_tilde[i] - a[i-1]*b_tilde[i-1]/b[i-1];
    cout << b_tilde[i]<<endl;
}
}
cout<< "        " <<endl;

/*
// prepare file for writing
outfile.open(outfilename);
outfile << fixed;
outfile << setprecision(4);
outfile << "Nummerical derivartion in task 1b" <<endl;
outfile << "n="<<"  "<< n << endl;;
outfile <<"Nummerical_solution   Analytical_solution        rel.error"<< endl;


//Backward substitution
v[n-1] = b_tilde[n-1]/b[n-1];

for(int i =n-2; i>=0 ; i--){
    v[i] = b_tilde[i]- c[i]*v[i+1]/b[i];
}
for(int i=0; i<n; i++){
    outfile << v[i] << "               " << deriv_u((i+1)*h) <<"                       "<<RelativeError(v[i],deriv_u((i+1)*h) ) << endl;
}

// ALl arrays are internal, have not x=0 and x=1, but everything inside.



outfile.close();
*/

// -----------------------------------------------------------
// 1c:
// -----------------------------------------------------------


for(int n = 10;n<11; n=n*10){
    //int n = 10;
    double *a, *b, *c, *b_tilde, *v;
    //a = new double[n];
    b = new double[n];
    c = new double[n];
    b_tilde = new double[n];
    v = new double[n];
    double h = 1.0/(double (n+1));
    string outfilename = "../Plotting_kjetil/1b_n_";

    outfilename += to_string(n);
    outfilename += ".txt";

//Filling in the diagonals and b_tilde as we have in 1b
for(int i =0; i<n;i++){
    //a[i] = -1;
    b[i] = 2;
    c[i] = -1;
    b_tilde[i] = function_f((i+1)*h)*h*h;
    }

int flops;
flops = 4*(n-1);

//Forward substitution

for(int i =1; i<n;i++){
    b[i] = (i+1.0)/i;
    b_tilde[i] = b_tilde[i] + ((i-1.0)*b_tilde[i-1])/double(i);
    cout <<b_tilde[i]<<endl;
}




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
