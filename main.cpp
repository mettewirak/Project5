#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <cstdlib>
#include <string>
using namespace std;

ofstream ofile;  // To contain histogram data
ofstream ofile2; // To contain the mean values of money (sorted)

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

void Transaction(double *values, int agent_1, int agent_2);                         // From 5a)
void Transaction(double *values, int agent_1, int agent_2, double lambda);          // From 5c)
int P_with_alpha(double M1, double M2, double alpha);                               // From 5d)
int P_with_alpha_gamma(double M1, double M2, double alpha, double gamma, int *c, int *c2);
int cmp(const void *x, const void *y);
void print_histogram(int N, int total_time, int total_runs, int *histogram, int money_steps_per, int money_steps_total,int alpha, int gamma, double lambda);
void print_averages(int N, int total_time, int total_runs, double *average_bank,int alpha, int gamma, double lambda);



int main()
{
    int N = 1000;

    int initial_money = 1; // Expressed in multiples of M_0
    double bank[N];
    double average_bank[N];


    int agent_1 = 0, agent_2 = 1;

    int total_runs = pow(10,4);
    int total_time = pow(10,7);

    double lambda = 0.5;
    double alpha = 1;
    double gamma = 0;
    double mu=0;

    int money_steps_per = 100; // Per unit aka d_money = 0.01
    int maximum_money = initial_money*100; // Guessed value, to inlcude most of the measurements.
    int money_steps_total = money_steps_per*maximum_money;

    // cout << "Number of possible values: " << money_steps_total << endl;

    int histogram[money_steps_total];



    int v = 0;
    int c[N*N]; // Number of interactions between two agents.
    int *C1;
    int *C2;
    for(alpha=2;alpha>0;alpha--){
    for(gamma =0; gamma<5; gamma++){
        for(int index = 0; index < money_steps_total; index++){
            histogram[index] = 0;
        }
        for(int agent = 0; agent < N; agent++){
            bank[agent] = initial_money;
            average_bank[agent] = 0;
        }
        for(int index = 0; index < money_steps_total; index++){
            histogram[index] = 0;
        }
        for(int current_run = 0; current_run < total_runs; current_run++){

            for(int agent = 0; agent < N; agent ++){
                bank[agent] = initial_money;
                for(int agent2 = 0; agent2 < N; agent2++){
                    c[agent+agent2*N] = 0;
                }
            }
            for(int year=0; year<4;year++){
            for(int current_time = 0; current_time < total_time/4; current_time++){
                for(int agent=0;agent<N;agent++){
                    bank[agent]=bank[agent]-(bank[agent]-1)*mu;

                }
                agent_1 = RandomNumberGenerator(gen)*N;
                agent_2 = RandomNumberGenerator(gen)*N;

                while(agent_1 == agent_2){
                    agent_2 = RandomNumberGenerator(gen)*N;
                }

                // 5a and 5b:
                //Transaction(bank, agent_1, agent_2);


                // 5c:
                //Transaction(bank, agent_1, agent_2, lambda);


                // 5d:
                /*if(P_with_alpha(bank[agent_1],bank[agent_2], alpha)==1){
                Transaction(bank, agent_1, agent_2, lambda);
            }*/


                // 5e:
                C1 = &c[agent_1+agent_2*N];
                C2 = &c[agent_2+agent_1*N];
                if(P_with_alpha_gamma(bank[agent_1],bank[agent_2], alpha, gamma,C1,C2) ==1){
                    c[agent_1+agent_2*N] = c[agent_1+agent_2*N]+1;
                    c[agent_2+agent_1*N] = c[agent_2+agent_1*N]+1;
                    Transaction(bank, agent_1, agent_2, lambda);
                }
            }
            }
            cout << "Run " << (current_run+1) << " finished.\n";

            // Sort the array, so that the averages of highest/lowerst(...) moneys can be found
            qsort(bank, sizeof(bank)/sizeof(bank[0]), sizeof(bank[0]), cmp);
            for(int agent = 0; agent < N; agent++){
                average_bank[agent] += (bank[agent]/(double)total_runs);
            }
        }

        // Checking that the amount of money in the system is conserved.
        double sum=0;
        for(int agent = 0; agent < N; agent++){
            sum += bank[agent];
            histogram[(int)(money_steps_per*bank[agent])] ++;
        }
        cout << "Penger i systemet = " << sum << " i snitt = " << sum/N << endl;
        /*
    for(int i = 0; i < N; i++){
        for(int k = 0; k < N; k++){
            cout << c[i+k*N] << " ";
        }
        cout << endl;
    }*/

        print_histogram(N, total_time, total_runs, histogram, money_steps_per, money_steps_total, alpha, gamma, lambda);
        print_averages(N, total_time, total_runs, average_bank, alpha, gamma, lambda);
    }}
    return 0;
}


void Transaction(double *values, int agent_1, int agent_2){

    double epsilon = RandomNumberGenerator(gen);
    double sum = values[agent_1] + values[agent_2];
    values[agent_1] = epsilon*sum;
    values[agent_2] = (1-epsilon)*sum;
}


void Transaction(double *values, int agent_1, int agent_2, double lambda){

    double epsilon = RandomNumberGenerator(gen);
    double sum = values[agent_1] + values[agent_2];
    values[agent_1] = lambda*values[agent_1]+ epsilon*(1-lambda)*sum;
    values[agent_2] = lambda*values[agent_2]+ (1-epsilon)*(1-lambda)*sum;
}


int P_with_alpha(double M1, double M2, double alpha){

    double p=pow((fabs(M1-M2)),-alpha);
    double r=RandomNumberGenerator(gen);

    if(p>r){
        return 1;
    }
    else{
        return 2;
    }
}


int P_with_alpha_gamma(double M1, double M2, double alpha, double gamma, int *c, int *c2){

    double p=pow((fabs(M1-M2)),-alpha)*pow((1+*c),gamma);
    double r=RandomNumberGenerator(gen);

    if(p>r){
        return 1;
    }
    else{
        return 2;
    }
}


int cmp(const void *x, const void *y) // Sorteringsfunksjon. Hentet fra https://stackoverflow.com/questions/8448790/sort-arrays-of-double-in-c.
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}


void print_histogram(int N, int total_time, int total_runs, int *histogram, int money_steps_per, int money_steps_total,int alpha, int gamma, double lambda){

    string filename = "Histogram N=" + to_string(N) + " time=" + to_string(int(log10(total_time))) + " runs=" + to_string(int(log10(total_runs)))+ " gamma="+ to_string(gamma)+" alpha="+ to_string(alpha)+" lambda="+ to_string(lambda)+ ".txt";
    ofile.open(filename);

    for(int index = 0; index < money_steps_total; index++){

        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << (double)index/money_steps_per;
        ofile << setw(15) << setprecision(8) << histogram[index];
        ofile << setw(15) << setprecision(8) << ((double)histogram[index] / (double)(total_runs));
        ofile << setw(15) << setprecision(8) << log10(((double)histogram[index] / (double)(total_runs)));
        ofile << endl;
    }
    ofile.close();
}


void print_averages(int N, int total_time, int total_runs, double *average_bank, int alpha, int gamma, double lambda){

    string filename = "Gjennomsnittsformuer N=" + to_string(N) + " time=" + to_string(int(log10(total_time))) + " runs=" + to_string(int(log10(total_runs)))+ " gamma="+ to_string(gamma)+" alpha="+ to_string(alpha)+" lambda="+ to_string(lambda*10)+ ".txt";
    ofile2.open(filename);
    cout<<lambda<<endl;
    for(int agent = 0; agent < N; agent++){
        ofile2 << setiosflags(ios::showpoint | ios::uppercase);
        ofile2 << setw(15) << setprecision(8) << agent + 1;
        ofile2 << setw(15) << setprecision(8) << average_bank[agent] << endl;

    }
    ofile2.close();
}
