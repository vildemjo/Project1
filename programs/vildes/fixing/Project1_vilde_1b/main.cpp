#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include <Windows.h>

using namespace std;

double source_term(double x);
double closed_form_solution(double x);
double RelativeError(double v, double u);
double* Calculated_solution(int n, double* v, double* x, double h);
void print_to_file(int n, double *v, double *u);
double get_cpu_time();

int main(){

    ofstream outfile;
    ofstream outfile_latex;

    outfile.open("error.txt");
    outfile_latex.open("error_table.tex");

    outfile_latex << fixed;
    outfile_latex << setprecision(2);

    outfile << "log10(h):" << "    " <<" MaxError[log10(RelativeError)]:" << endl;
    outfile_latex << "log(h):" << " & " << "log(RelativeError):" << " \\\\ " << "\\hline" << endl;

    for(int i= 1; i<=7; i++){
        int n = pow(10,i);//10^i;
        double h = 1.0/((double)(1 + n));
        double* x = new double[n];
        double* v = new double[n];
        double* u = new double[n];

        // Making a position array with n points
        // Will add the two known boundary points when writing to file
        x[0] = h;
        u[0] = closed_form_solution(x[0]);
        for(int i = 1; i<n; i++){
            x[i] = x[i-1] + h;
            u[i] = closed_form_solution(x[i]);
        }

        Calculated_solution(n, v, x, h);

        // print_to_file(n, v, u);
/*
        double MaxError = 100.0;
        double Error;

        for(int i = 1; i<n; i++){
            Error = RelativeError(v[i], u[i]);
            if(abs(MaxError)>Error){
                MaxError = Error;
            }
        }

        outfile << log10(h) << "    " << MaxError << endl;
        outfile_latex << log10(h) << " & " << MaxError << " \\\\ " << " \\hline" <<endl;
*/
    }

    outfile.close();

    cout << "CPU time: " << get_cpu_time() << " s";

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
    double* b = new double[n];
    double* new_b = new double[n];

    // The matrix values // new means after backward substitution
    double* e = new double[n];     // this is the superdiagonal and the subdiagonal - symmetrical
    double* d = new double[n];
    double* new_d = new double[n];

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
        v[i] = (new_b[i] -e[i]*v[i+1])/(double) new_d[i];
    }
    return v;
}

void print_to_file(int n, double* v, double *u){

    ofstream outfile;

    string str = to_string(n);
    string filename = string("result_n_") + str + string(".txt");

    outfile.open(filename);

    outfile << fixed;
    outfile << setprecision(5);

    // setting the boundary values
    double v_0 = 0.0;
    double u_0 = v_0;
    double v_n = v_0;
    double u_n = v_0;

    outfile << "Calculated solution:" << " , " << "Analytical solution:" << endl;
    outfile << v_0 << "                    " << u_0 << endl;
    for(int i = 0; i<n;i++){
        outfile << v[i] << "                    " << u[i] << endl;
    }
    outfile << v_n << "                    " << u_n << endl;
    outfile.close();
}

double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}
