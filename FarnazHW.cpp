#include <iostream>
#include <math.h>
#include <vector>
#include "../eigen-git-mirror/Eigen/Dense"
#include "../eigen-git-mirror/Eigen/LU"
#include <fstream>
#include <stdlib.h>

//#include "matplotlib-cpp/matplotlibcpp.h"
using namespace std;

//this function calculates the current energy of the matrix
double Cal_Canonical_E(Eigen::MatrixXd& mat1, double mu, double V0, double vnn, int size)
{
	double CE=V0;
for (int i=0;i<size;i++)
	for(int j=0; j<size;j++)
		CE=CE+vnn*mat1(i,j)*(mat1(i,(j+1)%size)+mat1((i+1)%size,j));
return CE;
}


//this function calculates the change in E if an elemnt in position (i,j) is flipped
double deltaE(Eigen::MatrixXd& mat1, double Vnn, int i, int j, int size)
{
	double dE;
    int up,down,left,right;
    if (i-1<0)
       up =  mat1(size-1,j);
    else
        up = mat1(i-1,j);

    if (j-1<0)
        left = mat1(i,size-1);
    else
        left = mat1(i,j-1);

    down = mat1((i+1)%size,j);
    right = mat1(i,(j+1)%size);

    dE = -2*Vnn*mat1(i,j)*(up + down + left + right);
return dE;
}

//inja & mizarim chon mikhaym har taghiri k to mat1 dade mishe to tabe to barname main ham ghabele royat bashe
double calculating_N(Eigen::MatrixXd& mat1, int size, double Vnn, double mu,double T, double& Cmu)
{
double Kb = 8.61733e-5;
double dE, dOmega,thermalexcite;
int i,j,i_rand,j_rand;

//ER keeps values of E, OR keeps values of Omega, NN keeps values of N during iterations
vector<double> ER;
vector<double> OR;
vector<double> NN;
//ER is keeping the values of Energy over runs and OR is keeping the values of Omega during runs

double V0=-2*Vnn;
double N = 0;
//at the beginning we calculate the number of ones to get the initial value of N
for(i=0;i<size;i++)
    for(j=0;j<size;j++)
        if(mat1(i,j)==1)
            N=N+1;
//let's first calculate the initial value of energy as well. We need it to calculate Omega.
double E=Cal_Canonical_E(mat1, mu, V0, Vnn, size);
ER.push_back(E);
OR.push_back(E-mu*N);
NN.push_back(N);

double dN=0;
int run;
int num_run=3000*size*size;
//at each run we generate 2 random numbers, calculate DE, from that DOmega and check the condifiton
for (int run=1; run<num_run; run++)
	{
	i_rand = int(rand()% size);
	j_rand = int(rand()% size);
	dE = deltaE(mat1, Vnn, i_rand, j_rand, size);
	dN=-mat1(i_rand,j_rand);
	dOmega = dE - (mu*dN) ;
	thermalexcite= exp (-dOmega / (Kb * T));
	if ((dOmega<0)||(thermalexcite>(((double) rand()) / ((double) RAND_MAX))))
		{
		//if the conditions are met, we update the values of E and Omega and push them back
		ER.push_back(ER[run-1]+dE);
		OR.push_back(OR[run-1]+dOmega);
		mat1(i_rand,j_rand)=dN;
		N=N+dN;
		NN.push_back(N);
		}
	else
		{
    //if condition is not met, we push back the previous values of these parameters
		ER.push_back(ER[run-1]);
		OR.push_back(OR[run-1]);
		NN.push_back(NN[run-1]);
		}
	}
//now we can calculate the averages on N on the last second half of Ns and the variance of Omegas
double OmegaBar=0;
double OmegaSquare=0;
N=0;
for(i=round(OR.size()/2);i<OR.size();i++){
    OmegaBar=OmegaBar+OR[i];
    N=N+NN[i];
	}
N=N/(OR.size()/2);
OmegaBar=OmegaBar/(OR.size()/2);
//here we calculate the variance of Omega on the second half of Omegas
for(i=round(OR.size()/2);i<OR.size();i++)
    OmegaSquare=OmegaSquare+pow(OR[i]-OmegaBar,2);

//since Cmu is defined as like pointer, we can change it here and see it in the main
Cmu=OmegaSquare/((Kb*T*T)*(OR.size()/2));

return N/(size*size);
}

int main()
{
    vector<vector<double>> xn_plot;
    double NC;
    int size=100;
    Eigen::MatrixXd mat1;
    double Vnn=0.030;
    double Cmu;

    for(double mu_given=-0.5; mu_given<=0.5; mu_given=mu_given+0.05){
    	//for each mu, we reset the matrix to start over
    	mat1=-1*Eigen::MatrixXd::Ones(size,size);
        vector<double> xn_temp;
        for (double T=100;T<=1100;T=T+100)
        	{
        	//for each temperature, it uses the previous mat1 sinse it is used as pointer,
        	//it gets Cmu as input but updates it again because it is pointer.
        	NC=calculating_N(mat1,size,Vnn,mu_given,T,Cmu);
        	//first for each mu, we continue pusshing back T, NC, and Cmu per each T in a vector,
        	//the push this vector in another vector one row per each mu.
        	xn_temp.push_back(T);
        	xn_temp.push_back(NC);
        	xn_temp.push_back(Cmu);
        	cout<<NC<<endl;
        	}
        	xn_plot.push_back(xn_temp);
        }

//we write this vector of vector to file to plot with MATLAB
ofstream results;
results.open("C:/Users/Farnood/Desktop/homeworkmatlab/results3.txt");

//writing results to file
for (int i = 0; i < xn_plot.size(); i++)
{
        for (int j = 0; j < xn_plot[i].size(); j++)
                results << xn_plot[i][j]<<'\t';
        results<<endl;
}
results.close();
}


