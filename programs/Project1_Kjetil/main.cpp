#include <armadillo>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>

#include <iomanip>
using namespace std;
using namespace arma;


double function_f(double x);
double deriv_u(double x);
double RelativeError(double v, double u);
mat general_solver(int n,mat a, mat c, mat b, mat v, mat f_tilde);
mat spess_solver(int n, mat b, mat f_tilde, mat v);

void writing(int n, string oppg, double h, mat u);
mat LU_solver(mat w, mat A);

ofstream outfile;


int main(){

    //double *a, *b, *c, *f_tilde, *v, **A, *ar;

    int n_maks;

    n_maks = 7;

    mat tid_1b = ones<vec>(n_maks);
    mat tid_1c = ones<vec>(n_maks);
    mat tid_lu = ones<vec>(n_maks);
    mat feil = ones<vec>(n_maks);

    for(int N = 1; N< n_maks+ 1; N++){
        int n = pow(10,N);
        mat a = ones<vec>(n);a[0] = 0;
        mat b = (ones<vec>(n))*(-2);
        mat c = ones<vec>(n);c[n-1] = 0;
        mat f_tilde = zeros<vec>(n);
        mat v = ones<vec>(n);
        double u;
        //mat x = ones<vec>(n);
        //mat u_analytic = ones<vec>(n);
        double h = 1.0/(double (n+1));


    //Filling in the diagonals and f_tilde as we have in 1b
    for(int i =0; i<n;i++){
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        if (i == 0) a[i] = 0;
        if (i==n-1) c[i] = 0;
        f_tilde[i] = function_f((i+1)*h)*h*h;

    }

    clock_t start1, finish1;
    start1 = clock();
    v = general_solver( n, a, c,  b,  v,  f_tilde);
    finish1 = clock();
    tid_1b[N-1] = (double) (finish1 - start1)/(CLOCKS_PER_SEC );


    clock_t start, finish;
    start = clock();
    v = spess_solver( n,  b,  f_tilde,  v);
    finish = clock();
    tid_1c[N-1] = (double) (finish - start)/(CLOCKS_PER_SEC );

    //cout << setprecision(10) << setw(20) << "Time used  for vector addition=" << timeused  << endl;


    //armadillo solves everything
if(N<=4){
    mat A = zeros<mat>(n,n);
      for (int i=0; i<n; i++) {
        A(i,i)=b[i];
        if (i!=n-1) A(i,i+1) = c[i];
        if (i!=n-1) A(i+1,i) = a[i+1];
    }

    //vec ar = solve(A,f_tilde);

    clock_t start2, finish2;
    start2 = clock();

    v = LU_solver(f_tilde, A);
    finish2 = clock();
    tid_lu[N-1] = (double) (finish2 - start2)/(CLOCKS_PER_SEC );

}
else{tid_lu[N-1] =0;}

    string filnavn = "../Plotting_kjetil/figures/1c_n_";
    filnavn += to_string(n);
    double max_error = 0;
    double errorr;
    // RelativeError(double v, double u){
    for (int i=0;i<n;i++){
        u = deriv_u((i+1)*h);
        errorr = RelativeError(v[i], u);
        if (abs(errorr) >max_error)max_error = abs(errorr);
    }
    feil[N-1] = log10(max_error);

    writing( n, filnavn,  h,  v);
    cout << N<< endl;
}

/*
    outfile.open("../../tex_files/figures/feil.tex");
    outfile << fixed;
    outfile << setprecision(3);
    outfile << "n & relative error  \\\\ " << endl;
    outfile << "\\hline  " <<endl;

    for(int i=0; i<n_maks; i++){
        outfile <<  "1e" << i+1 <<"    &         " << feil[i] << "\\\\ " <<       endl;
    }

    outfile.close();
*/



outfile.open("../../tex_files/figures/tid1.tex");
outfile << fixed;
outfile << setprecision(1);
outfile <<  std::scientific;

//outfile << "Time comparisson 1b, 1c" <<endl;
outfile << "n & general($ s$)  & special($ s$) & LU(s)  \\\\ " << endl;
outfile << "\\hline  " <<endl;

for(int i=0; i<n_maks; i++){
    outfile <<  "1e" << i+1 <<"    &         "<<tid_1b[i]<<"      &       "<< tid_1c[i]<<  "  &"<< tid_lu[i]   << "\\\\ " <<       endl;
}

outfile.close();

outfile.open("../../tex_files/figures/feil.tex");
outfile << fixed;
outfile << setprecision(2);
//outfile << "Time comparisson 1b, 1c" <<endl;
outfile << "n  &   max rel. error \\\\ " << endl;
outfile << "\\hline  " <<endl;


for(int i=0; i<n_maks; i++){
    outfile <<  i+1 <<" & "<< "  &"<< feil[i]<< " \\\\ "<<  endl;
}

outfile.close();


return 0;



}


mat spess_solver(int n, mat b, mat f_tilde, mat v){

for(int i =1; i<n+1;i++){
    b[i-1] = (i+1.0)/(i);
    f_tilde[i-1] = f_tilde[i-1] + (i-1)*f_tilde[i-2]/double(i);
}


v[n-1] = f_tilde[n-1]/b[n-1];
for(int i =n; i>=2; i--){
    v[i-2] = double(i-1)/(i)*(f_tilde[i-2]+ v[i-1]);
    //v_[i] = (f_tilde_[i] + v_[i+1])/b_[i];

}

return v;
}

mat general_solver(int n,mat a, mat c, mat b, mat v, mat f_tilde){
//Forward substitution

for (int i = 1; i<n;i++){
    double alpha = a[i]/b[i-1];
    b[i] = b[i] - alpha*c[i-1];
    f_tilde[i] = f_tilde[i] - alpha*f_tilde[i-1];
}

//Backward substitution
v[n-1] = f_tilde[n-1]/b[n-1];
for(int i =n-2; i>=0 ; i--){
    v[i] = (f_tilde[i]- c[i]*v[i+1]) / b[i];
}
return v;
}

void writing(int n, string filnavn, double h, mat u){

        string outfilename = "../Plotting_kjetil/";
        outfilename += filnavn;
        outfilename += ".txt";

        // prepare file for writing
        outfile.open(outfilename);
        outfile << fixed;
        outfile << setprecision(4);
        outfile << "Nummerical derivartion in task" << filnavn <<endl;
        outfile << "n="<<"  "<< n << endl;;
        outfile <<"Nummerical_solution   Analytical_solution        rel.error"<< endl;


        for(int i=0; i<n; i++){
            outfile << u[i] << "               " << deriv_u((i+1)*h) <<"                       "<<RelativeError(u[i],deriv_u((i+1)*h) ) << endl;
        }

        // ALl arrays are internal, have not x=0 and x=1, but everything inside.



        outfile.close();


}

mat LU_solver(mat w, mat A){
    mat L,U;
    lu(L,U,A);
    mat Y = solve(L,w);
    mat X = solve(U,Y);
    return X;}

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

