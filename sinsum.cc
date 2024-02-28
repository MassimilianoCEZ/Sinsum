#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

// Struct to hold parameters
struct Parameters {
    string type;
    int nbN;
    double tmin;
    double tmax;
    double ampMin;
    double ampMax;
    int nbL;
};

// constants and piece of code to be used for the project SINUSSUM
const double EPSIL_DICHO(1e-9);
const double EPSIL_T(1e-13);

// Error messages
const string BAD_SIGNAL("Error: signal type must be SAWTOOTH, SQUARE or TRIANGLE");
const string NBN_TOO_SMALL("Error: the number of sine terms must be greater than 0");
const string NBL_TOO_SMALL("Error: the number of lines must be greater than 2");
const string TIME_MIN_MAX("Error: time max is not bigger than min");
const string SIGNAL_MIN_MAX("Error: signal max is not bigger than min");
const string WRONG_TIME_VAL("Error: both time values must belong to [0., 1.]");
const string NBL_NOT_ODD("Error: the number of lines must be odd");

// Custom Prototype Functions
void print_error(string message);
void error(Parameters params);
int test_nbN(int nbN);
int test_time_min_max(double tmin, double tmax);
int test_amp_min_max(double ampMin, double ampMax);
int test_n_rows(int rows);
int test_type(string type);
/////////////////////////////////////////

bool check_interval(double lhs, double epsilon, double rhs);
void matrix_chosen(vector<vector<char>>& matrix, Parameters given, double epsil_t );
vector<vector<char>> empty_matrix(int nrows);
int time_index_i(double ampmin, double ampmax, int nrows );
void assign_time_axe(Parameters given, vector<vector<char>>& matrix);
int row_index(double function_value, double min, double deltaS);
void sawtooth_theory(vector<vector<char>>& matrix, Parameters given, double epsil_t);
void square_theory(vector<vector<char>>& matrix,Parameters given, double epsil_t);
void triangular_theory(vector<vector<char>>& matrix, Parameters given, double epsil_t);
double value_dent_de_scie(double t);
double value_square(double t);
double value_triangular(double t, double epsil_t);
double value_square_approx(double t, int nbN );
void functions_approx(vector<vector<char>>& matrix, Parameters given, double epsil_t);
double value_sawtooth_approx(double t, int nbN );
double value_triangle_approx(double t, int nbN );
void print(vector<vector<char>> matrix_input, int nbL, int nbC);
void print_bars(int nbC);
double approximated_matrix_chosen(string type, double t, int nbN);
////////////////////////////////////////////////////

vector<double> time_dichotomy(int nbN, string type);
double max_dicho_research(double t_start, double t_finish, string type, int nbN);

int main() {
    Parameters params;
    cin >> params.type >> params.nbN >> params.tmin >> params.tmax 
    >> params.ampMin >> params.ampMax >> params.nbL;
    vector<vector<char>> matrix_start;
    int nbC;
    nbC = 2 * params.nbL - 1;
    error(params);
    matrix_start = empty_matrix(params.nbL);
    assign_time_axe(params, matrix_start);
    matrix_chosen(matrix_start, params, EPSIL_T);
    print_bars(nbC);
    print(matrix_start, params.nbL, nbC);
    print_bars(nbC);
    vector<double> time_research = time_dichotomy(params.nbN, params.type);
    cout << setprecision(8) << fixed ;
    double max;
    max=max_dicho_research(time_research[0],time_research[1],params.type, params.nbN);
    cout << max << endl;
    return 0;
}



//print (task 2)

void print(vector<vector<char>> matrix_input, int nbL, int nbC){
    for (int i = 0; i < nbL; ++i) {
        for (int j = 0; j < nbC; ++j) {
            cout << matrix_input[i][j];
        }
        cout << endl;
    }

}

void print_bars(int nbC){
    for(int i(0); i < nbC; ++i){
        cout << "-" ;
    }
    cout << endl;
}


// Function Tache 1
//////////////////////////////////////////////////////////////////////////////////////


void error(Parameters params) {
    int condition_type;
    condition_type = test_type(params.type);
    if (condition_type == 0){
        print_error(BAD_SIGNAL);
    }
    int condition_nbN;
    condition_nbN = test_nbN(params.nbN);
    if (condition_nbN == 1) {
        print_error(NBN_TOO_SMALL);
    }
    int condition_t_min_max;
    condition_t_min_max = test_time_min_max(params.tmin, params.tmax);
    switch (condition_t_min_max) {
        case 1:
            print_error(TIME_MIN_MAX);
            break;
        case 2:
            print_error(WRONG_TIME_VAL);
            break;
    }
    int condition_amp_min_max;
    condition_amp_min_max = test_amp_min_max(params.ampMin, params.ampMax);
    if (condition_amp_min_max == 1) {
        print_error(SIGNAL_MIN_MAX);
    }
    int condition_n_rows;
    condition_n_rows = test_n_rows(params.nbL);
    switch (condition_n_rows) {
        case 1:
            print_error(NBL_TOO_SMALL);
            break;
        case 2:
            print_error(NBL_NOT_ODD);
            break;
    }
}

int test_nbN(int nbN) {
    if (nbN <= 0) {
        return 1; // 1 false, error
    }
    return 0;
}

int test_time_min_max(double tmin, double tmax) {
    if (tmin >= tmax) {
        return 1;
    }
    if ((tmax < 0. || tmin < 0.) || ((tmax > 1.0) || (tmin > 1.0))) {
        return 2;
    }
    return 0;
}

int test_amp_min_max(double ampMin, double ampMax) {
    if (ampMin >= ampMax) {
        return 1;
    }
    return 0;
}

int test_n_rows(int rows) {
    if (rows <= 2) {
        return 1;
    }
    if (rows % 2 == 0) {
        return 2;
    }
    return 0;
}

int test_type(string type){
    if (type == "SAWTOOTH"){
        return 1;
    }

    if (type == "SQUARE"){
        return 2;
    }

    if (type == "TRIANGLE"){
        return 3;
    }
    return 0;
}


void print_error(string message) {
    cout << message;
    cout << endl;
    exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////

// Functions Task 2

//////////////////////////////////////////////////////////////////////////////////////

bool check_interval(double lhs, double epsilon, double rhs){
    if((lhs - epsilon <= rhs)&&(lhs + epsilon >= rhs)){
        return 1;
    }
    return 0;
}




vector<vector<char>> empty_matrix(int nrows){
    int ncols;
    ncols = 2 * nrows - 1;
    vector<vector<char>> matrix_start(nrows, vector<char>(ncols, ' '));
    return matrix_start;
}


void matrix_chosen(vector<vector<char>>& matrix, Parameters given, double epsil_t ){
    if(given.type == "SAWTOOTH"){
        sawtooth_theory(matrix,given,epsil_t);
    }
    if(given.type == "SQUARE"){
        square_theory(matrix,given,epsil_t);
    }
    if(given.type == "TRIANGLE"){
        triangular_theory(matrix,given,epsil_t);
    }
    functions_approx(matrix,given,epsil_t);
}


int time_index_i(double ampmin, double ampmax, int nrows ){
    double deltas;
    deltas = (ampmax - ampmin) / (nrows - 1);
    for(int i(0); i < nrows; ++i){
        if( (ampmax - (deltas / 2 ) <= 0) && (0 < ampmax + (deltas / 2 ) )){
            return i;
        }
        ampmax -= deltas;
    }
    return -1;
}

void assign_time_axe(Parameters given, vector<vector<char>>& matrix){
    int time_row = time_index_i(given.ampMin, given.ampMax, given.nbL);
    int nbC;
    nbC = 2 * given.nbL - 1;
    if (time_row != -1){
        for (int j = 0; j < nbC; ++j) {
            matrix[time_row][j] = '.';
        }
    }
}


int row_index(double function_value, double min, double deltaS){
    double v;
    v = (function_value - min) / deltaS + 0.5;
    int i;
    if (v >= 0.0){
        i = floor(v);
        if (v - i >= -0.5 && v - i < 0.5) {
            return i;   
        }
    } else {
        return -1;
    }

}

void sawtooth_theory(vector<vector<char>>& matrix,Parameters given,double epsil_t){
    int nbC = 2 * given.nbL - 1;
    double delta_t = (given.tmax - given.tmin) / (nbC - 1);
    double delta_s = (given.ampMax - given.ampMin) / (given.nbL - 1);
    double t(given.tmin);
    double s_t;
    int i;
    for(int j(0); t <= given.tmax; ++j){
        if(check_interval(t,epsil_t,0)||check_interval(t,epsil_t,1)){
            s_t = 0 ;
            i = given.nbL - 1 - row_index(s_t,given.ampMin,delta_s) ; 
            if (i >= 0 && i < given.nbL) { // Verify index before acceding matrix
                matrix[i][j] = '+';
            }
        } else {
            s_t = value_dent_de_scie(t);
            i = given.nbL - 1 - row_index(s_t, given.ampMin, delta_s);
            if (i >= 0 && i < given.nbL) { 
                matrix[i][j] = '+';
            }
        }
        t += delta_t;
    }

}


double value_dent_de_scie(double t){
    double y;
    y = 2 * t - 1;
    return y;
}


void square_theory(vector<vector<char>>& matrix,Parameters given,double epsil_t){
    int nbC = 2 * given.nbL - 1;
    double delta_t = (given.tmax - given.tmin) / (nbC - 1);
    double delta_s = (given.ampMax - given.ampMin) / (given.nbL - 1);
    double t(given.tmin);
    double s_t;
    int i;
    for(int j(0); t <= given.tmax; ++j){
        if( check_interval(t,epsil_t,0.) || check_interval(t,epsil_t,0.5) ||
         check_interval(t,epsil_t,1.) ){
            s_t = 0 ;
            i = given.nbL - 1 - row_index(s_t,given.ampMin,delta_s) ; 
            if (i >= 0 && i < given.nbL) {
                matrix[i][j] = '+';
            }
        } else {
            s_t = value_square(t);
            i = given.nbL - 1 - row_index(s_t, given.ampMin, delta_s);
            if (i >= 0 && i < given.nbL) { 
                matrix[i][j] = '+';
            }
        }
        t += delta_t;
    }
}



double value_square(double t){
    double y;
    if ((t > 0 ) && (t < 0.5)){
        y = 1;
    }
    if ((t > 0.5) && (t < 1)){
        y = -1;
    }
    
    return y;
}


void triangular_theory(vector<vector<char>>& matrix, Parameters given, double epsil_t){
    int nbC = 2 * given.nbL - 1;
    double delta_t = (given.tmax - given.tmin) / (nbC - 1);
    double delta_s = (given.ampMax - given.ampMin) / (given.nbL - 1);
    double t(given.tmin);
    double s_t;
    int i;
    for(int j(0); t <= given.tmax; ++j){
        if( check_interval(t,epsil_t,0.5) ){
            s_t = 1 ;
            i = given.nbL - 1 - row_index(s_t,given.ampMin,delta_s) ; 
            if (i >= 0 && i < given.nbL) { 
                matrix[i][j] = '+';
            }
        } else {
            s_t = value_triangular(t,epsil_t);
            i = given.nbL - 1 - row_index(s_t, given.ampMin, delta_s);
            if (i >= 0 && i < given.nbL) { 
                matrix[i][j] = '+';
            }
        }
        t += delta_t;
    }
}

double value_triangular(double t, double epsil_t){
    double y;
    if ((t >= 0) && (t < 0.5)){
        y = 4*t - 1;
    }
    if ((t >= 0.5) && (t <= 1)){
        y = -4*t + 3;
    }
    return y;
}


double approximated_matrix_chosen(string type, double t, int nbN){
        double s_t;
        if(type == "SQUARE"){
            s_t = value_square_approx(t,nbN);  
        }
        if(type == "TRIANGLE"){
            s_t = value_triangle_approx(t,nbN);
        }
        if(type == "SAWTOOTH"){
            s_t = value_sawtooth_approx(t,nbN);
        }
        return s_t;
}




void functions_approx(vector<vector<char>>& matrix, Parameters given, double epsil_t){
    int nbC = 2 * given.nbL - 1;
    double delta_t = (given.tmax - given.tmin) / (nbC - 1);
    double delta_s = (given.ampMax - given.ampMin) / (given.nbL - 1);
    double t(given.tmin);
    double s_t;
    int i;
    for(int j(0); t <= given.tmax; ++j){
        s_t = approximated_matrix_chosen(given.type, t, given.nbN);
        i = given.nbL - 1 - row_index(s_t, given.ampMin, delta_s);
        if (i >= 0 && i < given.nbL) { 
            matrix[i][j] = '*';
        }
        t += delta_t; 
    }
}


double value_square_approx(double t, int nbN ){
    double y = 0.0;
    for(int i(1); i <= nbN; ++i){
        y += sin(2.0 * M_PI * (2.0 * i - 1) * t ) / (2.0 * i - 1);
    }
    y = y * (4.0) / M_PI;
    return y;
}

double value_sawtooth_approx(double t, int nbN ){
    double y = 0.0;
    for(int i(1); i <= nbN; ++i){
        y += pow(-1.0,i)/i * sin(2.0 * M_PI * i * (t - 0.5) );
    }
    y = y * (-2.0) / M_PI;
    return y;
}

double value_triangle_approx(double t, int nbN ){
    double y = 0.0;
    for(int i(1); i <= nbN; ++i){
        y += pow(-1.0,i)/pow((2.0*i-1),2)*sin(2.0* M_PI * (2.0 * i - 1) * (t - 0.25) );
    }
    y = y * (-8.0) / pow(M_PI, 2);
    return y;
}

//////////////////////////////////////////////////////////////////////////////////////

// Functions Task 3

//////////////////////////////////////////////////////////////////////////////////////



vector<double> time_dichotomy(int nbN, string type) {
    vector<double> time_interval;

    if (type == "SAWTOOTH") {
        time_interval = {1.0 - 1.0 / (2 * nbN + 1), 1.0};
    } else if (type == "SQUARE") {
        time_interval = {0.0, 1.0 / (2 * nbN + 1)};
    } else if (type == "TRIANGLE") {
        time_interval = {0.5 - 1.0 / (2*(2 * nbN + 1)), 0.5 + 1.0 / (2*(2 * nbN + 1))};
    }

    return time_interval;
}


double max_dicho_research(double t_start, double t_finish, string type, int nbN){
    double t_avrg, t_avrg_old, f_start, f_finish, f_avrg, f_old_avrg;
    do{
        t_avrg_old = t_avrg;
        f_old_avrg = approximated_matrix_chosen(type, t_avrg_old, nbN);
        t_avrg =(t_start + t_finish) / 2.0; 
        f_start = approximated_matrix_chosen(type, t_start, nbN);
        f_avrg = approximated_matrix_chosen(type, t_avrg, nbN );
        f_finish = approximated_matrix_chosen(type, t_finish, nbN);
        if((f_avrg < f_start) && (f_avrg > f_finish)){
            t_start = t_start;
            t_finish = t_avrg;
        }
        if((f_avrg > f_start) && (f_avrg < f_finish)){
            t_start = t_avrg;
            t_finish = t_finish;
        }
        if((f_avrg>f_start)&&(f_avrg>f_finish)||(f_avrg < f_start)&&(f_avrg<f_finish)){
            if(f_start < f_finish){
                t_start = t_avrg;
            }
            if(f_start > f_finish){
                t_finish = t_avrg;
            }
        }
        if((f_avrg == f_start) || (f_avrg == f_finish)){
            if(f_start < f_finish){
                t_start = t_avrg;
            }
            if(f_start > f_finish){
                t_finish = t_avrg;
            }
        }
    }while(not(((f_avrg-f_old_avrg)-EPSIL_DICHO<0)&&
        ((f_avrg-f_old_avrg)+EPSIL_DICHO > 0)));
    return f_avrg;
}



















