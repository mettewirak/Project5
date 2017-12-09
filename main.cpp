#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <cstdlib>
#include <string>
#include <mpi.h>
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
void print_histogram(int N, int total_time, int total_runs, double lambda, int alpha, int gamma, int *histogram, int money_steps_per, int money_steps_total);
void print_averages(int N, int total_time, int total_runs, double lambda, int alpha, int gamma, double *average_bank);


int main(int argc, char* argv[])
{
    int numprocs, my_rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int N = 1000;

    int initial_money = 1; // Expressed in multiples of M_0
    double bank[N];
    double average_bank[N];
    double average_bank_MPI[N];
    int number_of_runs[numprocs];

    for(int i=0; i<numprocs; i++){
        number_of_runs[i] = 0;
    }

    for(int agent = 0; agent < N; agent++){
        bank[agent] = initial_money;
        average_bank[agent] = 0;
        average_bank_MPI[agent] = 0;
    }
    int agent_1 = 0, agent_2 = 1;

    int total_runs = pow(10,4);
    int runs_per_rank = total_runs/numprocs;
    int total_time = pow(10,7);

    int money_steps_per = 100; // Per unit aka d_money = 0.01
    int maximum_money = initial_money*200; // Guessed value, to inlcude most of the measurements.
    int money_steps_total = money_steps_per*maximum_money;

    int histogram[money_steps_total];
    int histogram_MPI[money_steps_total];
    for(int index = 0; index < money_steps_total; index++){
        histogram[index] = 0;
        histogram_MPI[index] = 0;
    }

    double lambda = 0;
    int alpha = 2;
    int gamma = 1;
    int v = 0;
    int c[N*N]; // Number of interactions between two agents.
    int *C1;
    int *C2;

    for(int current_run = 0; current_run < runs_per_rank; current_run++){

        for(int agent = 0; agent < N; agent ++){
            bank[agent] = initial_money;
            for(int agent2 = 0; agent2 < N; agent2++){
                c[agent+agent2*N] = 0;
            }
        }

        for(int current_time = 0; current_time < total_time; current_time++){

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
            //if(P_with_alpha(bank[agent_1],bank[agent_2], alpha)==1){
            //    Transaction(bank, agent_1, agent_2, lambda);
            //}


            // 5e:
            C1 = &c[agent_1+agent_2*N];
            C2 = &c[agent_2+agent_1*N];
            if(P_with_alpha_gamma(bank[agent_1],bank[agent_2], alpha, gamma,C1,C2) ==1){
                c[agent_1+agent_2*N] = c[agent_1+agent_2*N]+1;
                c[agent_2+agent_1*N] = c[agent_2+agent_1*N]+1;
                Transaction(bank, agent_1, agent_2, lambda);
            }
        }

        number_of_runs[my_rank]++;

        if(my_rank==0){
            cout << (current_run+1) << "/" << runs_per_rank << " finished\n";
        }

        // Sort the array, so that the averages of highest/lowerst(...) moneys can be found
        qsort(bank, sizeof(bank)/sizeof(bank[0]), sizeof(bank[0]), cmp);
        for(int agent = 0; agent < N; agent++){
            average_bank[agent] += (bank[agent]/(double)runs_per_rank);
        }

        for(int agent = 0; agent < N; agent++){
            histogram[(int)(money_steps_per*bank[agent])] ++;
        }

    }

    MPI_Allreduce(histogram, histogram_MPI, money_steps_total, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(average_bank, average_bank_MPI, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int agent = 0; agent<N; agent++){
        average_bank_MPI[agent] /= 8;
    }

    // Checking that the amount of money in the system is conserved, as well as updating histogram.
    double sum=0;
    for(int agent = 0; agent < N; agent++){
        sum += bank[agent];
    }
    //cout << "Gjennomsnittsformue: " << sum/N << endl;

    /*for(int i = 0; i < N; i++){
        for(int k = 0; k < N; k++){
            cout << c[i+k*N] << " ";
        }
        cout << endl;
    }*/

    if(my_rank==0){
        cout << "Inne i printe-if.\n";
        print_histogram(N, total_time, total_runs, lambda, alpha, gamma, histogram_MPI, money_steps_per, money_steps_total);
        print_averages(N, total_time, total_runs, lambda, alpha, gamma, average_bank_MPI);
    }

    cout << "Number of runs in rank " << my_rank << ": " << number_of_runs[my_rank] << endl;

    MPI_Finalize();

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

    double p = pow((fabs(M1-M2)),-alpha);
    double r = RandomNumberGenerator(gen);

    if(p > r){
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


void print_histogram(int N, int total_time, int total_runs, double lambda, int alpha, int gamma, int *histogram, int money_steps_per, int money_steps_total){

    string filename = "Histogram N=" + to_string(N) + " time=" + to_string(int(log10(total_time))) + " runs=" + to_string(int(log10(total_runs))) + " lambda=" + to_string(lambda) +" alpha=" + to_string(alpha) + " gamma= " + to_string(gamma) + ".txt";
    ofile.open(filename);

    for(int index = 0; index < money_steps_total; index++){

        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << (double)index/money_steps_per;
        ofile << setw(15) << setprecision(8) << histogram[index];
        ofile << setw(15) << setprecision(8) << ((double)(histogram[index]*money_steps_per) / (double)(N*total_runs));
        ofile << setw(15) << setprecision(8) << log10(((double)(histogram[index]*money_steps_per) / (double)(N*total_runs)));
        ofile << endl;
    }
    ofile.close();
}


void print_averages(int N, int total_time, int total_runs, double lambda, int alpha, int gamma, double *average_bank){

    string filename = "Gjennomsnittsformuer N=" + to_string(N) + " time=" + to_string(int(log10(total_time))) + " runs=" + to_string(int(log10(total_runs))) +  " lambda=" + to_string(lambda) + " alpha=" + to_string(alpha) + " gamma= " + to_string(gamma) + ".txt";
    ofile2.open(filename);

    for(int agent = (N-1); agent > -1; agent--){
        ofile2 << setiosflags(ios::showpoint | ios::uppercase);
        ofile2 << setw(15) << setprecision(8) << agent + 1;
        ofile2 << setw(15) << setprecision(8) << average_bank[agent];
        ofile2 << endl;
    }
    ofile2.close();
}
