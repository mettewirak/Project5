#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <cstdlib>
#include <string>
#include <sstream>
using namespace std;

ifstream infile;

double calculate(double m, double estimate){

    double analytic = exp(-m);
    return abs(analytic-estimate);
}


int main()
{
    cout << "Hello World!" << endl;

    infile.open("/uio/hume/student-u87/mettewir/Documents/Computational Physics/Project5/Histogram data/Histogram N=500 time=7 runs=1.txt");

    //double a, b, c, d;
    double sum = 0.0;
    int counter = 0;

    string line;
    while(getline(infile, line)){
        istringstream iss(line);
        double a, b, c;
        if(!(iss >> a >> b >> c)) {break;}

        cout << a << "  " << b << endl;
        sum += calculate(a, c);
        counter ++;
    }

    cout << "The sum of errors for all steps is: " << sum << endl;
    cout << "Number of steps: " << counter << endl;
    cout << "The error per step is: " << sum/counter << endl;

    return 0;
}
