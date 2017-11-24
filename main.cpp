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
int P_with_alpha(double M1,double M2,double alpha);
int P_with_alpha_gamma(double M1,double M2,double alpha, double gamma,int *c, int *c2);

int main()
{
    int N = 500;
    int c[N*N]; //number of interactions between two agents.
    int *C1;
    int *C2;

    int initial_money = 1; //Expressed in multiples of M_0
    double bank[N];
    double total_bank[N];
    double lambda=0;
    double alpha=2;
    double gamma=1;
    int v=0;
    for(int agent = 0; agent < N; agent++){
        bank[agent] = initial_money;
        total_bank[agent] = 0;
    }

    int agent_1 = 0, agent_2 = 1;

    int total_runs = pow(10,1);

    int total_time = pow(10,7);

    int money_steps_per = 100; // Per unit aka d_money = 0.1
    int maximum_money = initial_money*100; // Guessed value, to inlcude most of the measurements.
    int money_steps_total = money_steps_per*maximum_money;

    cout << "Number of possible values: " << money_steps_total << endl;

    int histogram[money_steps_total];
    for(int index = 0; index < money_steps_total; index++){
        histogram[index] = 0;
    }

    //ofile.open("Histogram time=7, runs=3, N=500, dm=0.01.txt");
    ofile.open("Test2.txt");

    for(int current_run = 0; current_run < total_runs; current_run++){

        for(int agent = 0; agent < N; agent ++){
            bank[agent] = initial_money;
            for(int agent2=0; agent2<N;agent2++){
            c[agent+agent2*N]=0;}
        }

        for(int current_time = 0; current_time < total_time; current_time++){

            agent_1 = RandomNumberGenerator(gen)*N;
            agent_2 = RandomNumberGenerator(gen)*N;

            while(agent_1 == agent_2){
                agent_2 = RandomNumberGenerator(gen)*N;
            }
            C1=&c[agent_1+agent_2*N];
            C2=&c[agent_2+agent_1*N];
            //if(P_with_alpha(bank[agent_1],bank[agent_2], alpha)==1){Transaction(bank, agent_1, agent_2, lambda);}
            if(P_with_alpha_gamma(bank[agent_1],bank[agent_2], alpha, gamma,C1,C2) ==1){c[agent_1+agent_2*N]=c[agent_1+agent_2*N]+1;c[agent_2+agent_1*N]=c[agent_2+agent_1*N]+1;Transaction(bank, agent_1, agent_2, lambda);}
            //Transaction(bank, agent_1, agent_2, lambda);
        }

        cout << "\n Run " << (current_run+1) << " finished.\n"<<v;

        for(int agent = 0; agent < N; agent++){
            total_bank[agent] += bank[agent];
        }
        double sum=0;
        for(int agent = 0; agent < N; agent++){
            sum+=bank[agent];
            histogram[(int)(money_steps_per*bank[agent])] ++;
        }
    cout<<"penger i systemet= "<<sum<<"i snitt= "<<sum/N<<endl;

    }

    cout << "\n Mean values over " << total_runs << "runs:\n";

    for(int agent = 0; agent < N; agent++){
        cout << "Agent " << (agent+1) << ": " << total_bank[agent]/total_runs << endl;

    }

    for(int index = 0; index < money_steps_total; index++){

        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << (double)index/money_steps_per;
        ofile << setw(15) << setprecision(8) << histogram[index];
        ofile << setw(15) << setprecision(8) << (double)histogram[index] / (double)(total_runs);
        ofile << endl;
    }

    for(int i=0;i<N;i++){
    for(int k=0;k<N;k++){
        cout<<c[i+k*N]<<" ";
    }
    cout<<endl;

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

int P_with_alpha(double M1,double M2,double alpha){

    double p=pow((fabs(M1-M2)),-alpha);
    double r=RandomNumberGenerator(gen);

    if(p>r){return 1;}
    else{return 2;}
}
int P_with_alpha_gamma(double M1,double M2,double alpha, double gamma,int *c, int *c2){

    double p=pow((fabs(M1-M2)),-alpha)*pow((1+*c),gamma);
    double r=RandomNumberGenerator(gen);

    if(p>r){return 1;}
    else{return 2;}


}


