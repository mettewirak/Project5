#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
using namespace std;

ofstream ofile;

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

void Transaction(double *values, int agent_1, int agent_2, double lambda);

int main()
{
    int N = 500;
    int initial_money = 1; //Expressed in multiples of M_0
    double bank[N];
    double total_bank[N];
    double lambda=0;
    for(int agent = 0; agent < N; agent++){
        bank[agent] = initial_money;
        total_bank[agent] = 0;
    }

    int agent_1 = 0, agent_2 = 1;

    int total_runs = pow(10,0);

    int total_time = pow(10,7);

    int money_steps_per = 100; // Per unit aka d_money = 0.1
    int maximum_money = initial_money*10; // Guessed value, to inlcude most of the measurements.
    int money_steps_total = money_steps_per*maximum_money;

    cout << "Number of possible values: " << money_steps_total << endl;

    int histogram[money_steps_total];
    for(int index = 0; index < money_steps_total; index++){
        histogram[index] = 0;
    }

    ofile.open("Histogram time=7, runs=0, N=500, dm=0.01.txt");


    for(int current_run = 0; current_run < total_runs; current_run++){

        for(int agent = 0; agent < N; agent ++){
            bank[agent] = initial_money;
        }

        for(int current_time = 0; current_time < total_time; current_time++){

                agent_1 = RandomNumberGenerator(gen)*N;
                agent_2 = RandomNumberGenerator(gen)*N;

                while(agent_1 == agent_2){
                    agent_2 = RandomNumberGenerator(gen)*N;
                }

                Transaction(bank, agent_1, agent_2, lambda);
        }

        cout << "\n Run " << (current_run+1) << " finished.\n";

        for(int agent = 0; agent < N; agent++){
            total_bank[agent] += bank[agent];
        }

        for(int agent = 0; agent < N; agent++){
            histogram[(int)(money_steps_per*bank[agent])] ++;
        }
    }

    cout << "\n Mean values over " << total_runs << "runs:\n";
    for(int agent = 0; agent < N; agent++){
        cout << "Agent " << (agent+1) << ": " << total_bank[agent]/total_runs << endl;
    }

    for(int index = 0; index < money_steps_total; index++){
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << (double)index/money_steps_per;
        ofile << setw(15) << setprecision(8) << histogram[index];
        ofile << setw(15) << setprecision(8) << (double)histogram[index] / (double)(total_runs*N);
        ofile << endl;
    }

    ofile.close();
    return 0;
}



void Transaction(double *values, int agent_1, int agent_2, double lambda){

    double epsilon = RandomNumberGenerator(gen);
    double sum = values[agent_1] + values[agent_2];
    values[agent_1] = lambda*values[agent_1]+ epsilon*(1-lambda)*sum;
    values[agent_2] = lambda*values[agent_2]+ (1-epsilon)*(1-lambda)*sum;
}
