/*C++ CODE - MANGEAT MATTHIEU - 2026*/
/*HYDRODYNAMIC THEORY OF ACTIVE POTTS MODEL UNDER A BIDIRECTIONAL FIELD*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <omp.h>

using namespace std;

//Personal libraries.
#include "lib/special_functions.cpp"

//////////////////////////////
///// MESHGRID FUNCTIONS /////
//////////////////////////////

void meshInit(vector<vector<double> > &RHO1, vector<vector<double> > &RHO2, vector<vector<double> > &RHO3, vector<vector<double> > &RHO4, const int &init, const double &rhog, const double &rhol, const double &rho0, const double &beta, const double &r, const int &NX, const int &NY)
{
	double magl=0.8, dr=0.25;
	if (init==1)
	{
		magl=0;
	}
	if (init==2)
	{
		dr=0.;
	}
	
	#pragma omp parallel for default(shared)
	for (int x=0; x<NX; x++)
	{
		for (int y=0; y<NY; y++)
		{
			if (x<NX/2)
			{
				double rho=rho0*(1+dr*sin(8*M_PI*y/NY));
				RHO1[x][y]=0.25*rho*(1+3*magl);
				RHO3[x][y]=0.25*rho*(1-magl);
				RHO2[x][y]=0.25*rho*(1-magl);
				RHO4[x][y]=0.25*rho*(1-magl);
			}
			else
			{
				double rho=rho0*(1-dr*sin(8*M_PI*y/NY));
				RHO3[x][y]=0.25*rho*(1+3*magl);
				RHO1[x][y]=0.25*rho*(1-magl);
				RHO2[x][y]=0.25*rho*(1-magl);
				RHO4[x][y]=0.25*rho*(1-magl);
			}
			
		}
	}
}

double Iflip(const double &rj, const double ri, const double &deltaH, const double &rho, const double &beta, const double &r)
{
	double mag=4*beta*(rj-ri)/rho+4*beta*deltaH;
	return exp(r/(2*rho))*((rj+ri-r/(4*beta))*sinh(mag)-(rj-ri)*cosh(mag));
}


void finiteDiff(vector<vector<double> > &RHO1, vector<vector<double> > &RHO2, vector<vector<double> > &RHO3, vector<vector<double> > &RHO4, const double &Dpara, const double &Dperp, const double &v0, const double &gamma0, const double &beta, const vector<vector<double> > &HH, const double &r, const int &NX, const int &NY)
{
	static vector<double> DRHO1(NX*NY,0), DRHO2(NX*NY,0), DRHO3(NX*NY,0), DRHO4(NX*NY,0);
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int x=0; x<NX; x++)
		{
			const int xm=x-1+(x==0)*NX, xp=x+1-(x==NX-1)*NX;
			
			const int reg=(x<NX/2) ? 0 : 1;
			const vector<double> HLOC=HH[reg];
			
			#pragma omp simd
			for (int y=0; y<NY; y++)
			{
				const int ym=y-1+(y==0)*NY, yp=y+1-(y==NY-1)*NY;
				const double rho=RHO1[x][y]+RHO2[x][y]+RHO3[x][y]+RHO4[x][y];
				
				const double I12=Iflip(RHO1[x][y],RHO2[x][y],HLOC[0]-HLOC[1],rho,beta,r);
				const double I13=Iflip(RHO1[x][y],RHO3[x][y],HLOC[0]-HLOC[2],rho,beta,r);
				const double I14=Iflip(RHO1[x][y],RHO4[x][y],HLOC[0]-HLOC[3],rho,beta,r);
				const double I23=Iflip(RHO2[x][y],RHO3[x][y],HLOC[1]-HLOC[2],rho,beta,r);
				const double I24=Iflip(RHO2[x][y],RHO4[x][y],HLOC[1]-HLOC[3],rho,beta,r);
				const double I34=Iflip(RHO3[x][y],RHO4[x][y],HLOC[2]-HLOC[3],rho,beta,r);
				
				DRHO1[x*NY+y]=Dpara*(RHO1[xm][y]+RHO1[xp][y]-2*RHO1[x][y])+Dperp*(RHO1[x][ym]+RHO1[x][yp]-2*RHO1[x][y]) - v0*(RHO1[xp][y]-RHO1[xm][y]) + gamma0*(I12+I13+I14);
				DRHO2[x*NY+y]=Dperp*(RHO2[xm][y]+RHO2[xp][y]-2*RHO2[x][y])+Dpara*(RHO2[x][ym]+RHO2[x][yp]-2*RHO2[x][y]) - v0*(RHO2[x][yp]-RHO2[x][ym]) + gamma0*(-I12+I23+I24);
				DRHO3[x*NY+y]=Dpara*(RHO3[xm][y]+RHO3[xp][y]-2*RHO3[x][y])+Dperp*(RHO3[x][ym]+RHO3[x][yp]-2*RHO3[x][y]) + v0*(RHO3[xp][y]-RHO3[xm][y]) + gamma0*(-I13-I23+I34);
				DRHO4[x*NY+y]=Dperp*(RHO4[xm][y]+RHO4[xp][y]-2*RHO4[x][y])+Dpara*(RHO4[x][ym]+RHO4[x][yp]-2*RHO4[x][y]) + v0*(RHO4[x][yp]-RHO4[x][ym]) + gamma0*(-I14-I24-I34);
			}
		}
		
		#pragma omp for schedule(static)
		for (int x=0; x<NX; x++)
		{
			#pragma omp simd
			for (int y=0; y<NY; y++)
			{
				RHO1[x][y]+=DRHO1[x*NY+y];
				RHO2[x][y]+=DRHO2[x*NY+y];
				RHO3[x][y]+=DRHO3[x*NY+y];
				RHO4[x][y]+=DRHO4[x*NY+y];
			}
		}
	}
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Total average on the all space.
double average(const vector< vector<double> > &RHO, const int &LX, const int &LY)
{
	double rhoAv=0;
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			rhoAv+=RHO[x][y];
		}
	}
	return double(rhoAv)/(LX*LY);
}

//Average along the y-axis.
vector<double> average1d_along_y(const vector< vector<double> > &RHO, const int &LX, const int &LY)
{
	vector<double> rhoAv(LX,0.);
	for (int x0=0; x0<LX; x0++)
	{
		for (int y0=0; y0<LY; y0++)
		{
			rhoAv[x0]+=double(RHO[x0][y0])/LY;
		}
	}
	return rhoAv;
}

//Average along the x-axis.
vector<double> average1d_along_x(const vector< vector<double> > &RHO, const int &LX, const int &LY)
{
	vector<double> rhoAv(LY,0.);
	for (int x0=0; x0<LX; x0++)
	{
		for (int y0=0; y0<LY; y0++)
		{
			rhoAv[y0]+=double(RHO[x0][y0])/LX;
		}
	}
	return rhoAv;
}

//State of the node
int max_state(const double &rho1, const double &rho2, const double &rho3, const double &rho4, double &rho)
{
	static double eps=1e-2;
	rho=rho1+rho2+rho3+rho4;
	
	if (rho1>(0.25+eps)*rho)
	{
		return 1;
	}
	else if (rho2>(0.25+eps)*rho)
	{
		return 2;
	}
	else if (rho3>(0.25+eps)*rho)
	{
		return 3;
	}
	else if (rho4>(0.25+eps)*rho)
	{
		return 4;
	}
	else
	{
		return 0;
	}
}

//Export density and state.
void exportDensityState(const vector< vector<double> > &RHO1, const vector< vector<double> > &RHO2, const vector< vector<double> > &RHO3, const vector< vector<double> > &RHO4, const double &beta, const double &epsilon, const double &rho0, const double &h, const int &LX, const int &LY, const int &NX, const int &NY, const int &init, const double &t)
{
	const int a=NX/LX;
	static vector<float> rho(LX*LY,0.);
	static vector<uint8_t> state(LX*LY,0);
	for (int y=0;y<LY;y++)
	{
		for (int x=0; x<LX; x++)
		{
			double RHO;
			state[y*LX+x]=max_state(RHO1[a*x][a*y],RHO2[a*x][a*y],RHO3[a*x][a*y],RHO4[a*x][a*y],RHO);
			rho[y*LX+x]=RHO;
		}
	}
	
	//Creation of the density file.
	int returnSystem=system("mkdir -p data_RFAPM_dynamics2d/");
	stringstream ssDensity;
	ssDensity << "./data_RFAPM_dynamics2d/RFAPM_density_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_t=" << t << ".bin";
	string nameDensity = ssDensity.str();
	ofstream fileDensity(nameDensity.c_str(),ios::binary);
	fileDensity.write(reinterpret_cast<const char*>(rho.data()), rho.size()*sizeof(float));
	fileDensity.close();
	
	//Creation of the state file.
	stringstream ssState;
	ssState << "./data_RFAPM_dynamics2d/RFAPM_state_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_t=" << t << ".bin";
	string nameState = ssState.str();
	ofstream fileState(nameState.c_str(),ios::binary);
	fileState.write(reinterpret_cast<const char*>(state.data()), state.size()*sizeof(uint8_t));
	fileState.close();
}


//Export averaged densities along y-axis.
void exportProfilesAlongX(const vector< vector<double> > &RHO1, const vector< vector<double> > &RHO2, const vector< vector<double> > &RHO3, const vector< vector<double> > &RHO4, const double &beta, const double &epsilon, const double &rho0, const double &h, const int &LX, const int &LY, const int &NX, const int &NY, const int &init, const double &t)
{
	const vector<double> r0=average1d_along_y(RHO1,NX,NY), r1=average1d_along_y(RHO2,NX,NY), r2=average1d_along_y(RHO3,NX,NY), r3=average1d_along_y(RHO4,NX,NY);
	
	int returnSystem=system("mkdir -p data_RFAPM_dynamics1d/");
	//Création du fichier de données.
	stringstream ss;
	ss << "./data_RFAPM_dynamics1d/RFAPM_profile_along_x_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_t=" << t << ".txt";
	string nameProfile = ss.str();
	
	ofstream fileProfile(nameProfile.c_str(),ios::trunc);
	fileProfile.precision(6);
	fileProfile << "#site\tRHO\tRHO0\tRHO1\tRHO2\tRHO3" << endl;
	
	for (int x0=0; x0<NX; x0++)
	{
		fileProfile << x0 << "\t" << r0[x0]+r1[x0]+r2[x0]+r3[x0] << "\t" << r0[x0] << "\t" << r1[x0] << "\t" << r2[x0] << "\t" << r3[x0] << endl;
	}
	fileProfile.close();
}

//Export averaged densities along x-axis.
void exportProfilesAlongY(const vector< vector<double> > &RHO1, const vector< vector<double> > &RHO2, const vector< vector<double> > &RHO3, const vector< vector<double> > &RHO4, const double &beta, const double &epsilon, const double &rho0, const double &h, const int &LX, const int &LY, const int &NX, const int &NY, const int &init, const double &t)
{
	const vector<double> r0=average1d_along_x(RHO1,NX,NY), r1=average1d_along_x(RHO2,NX,NY), r2=average1d_along_x(RHO3,NX,NY), r3=average1d_along_x(RHO4,NX,NY);
	
	int returnSystem=system("mkdir -p data_RFAPM_dynamics1d/");
	//Création du fichier de données.
	stringstream ss;
	ss << "./data_RFAPM_dynamics1d/RFAPM_profile_along_y_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_t=" << t << ".txt";
	string nameProfile = ss.str();
	
	ofstream fileProfile(nameProfile.c_str(),ios::trunc);
	fileProfile.precision(6);
	fileProfile << "#site\tRHO\tRHO0\tRHO1\tRHO2\tRHO3" << endl;
	
	for (int y0=0; y0<NY; y0++)
	{
		fileProfile << y0 << "\t" << r0[y0]+r1[y0]+r2[y0]+r3[y0] << "\t" << r0[y0] << "\t" << r1[y0] << "\t" << r2[y0] << "\t" << r3[y0] << endl;
	}
	fileProfile.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &epsilon, double &rho0, double &h, int &LX, int &LY, int &init, double &rhog, double &rhol, double &dx, double &dt, double &tmax, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-h=" ))
		{
			h=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-rhog=" ))
		{
			rhog=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-rhol=" ))
		{
			rhol=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-dt=" ))
		{
			dt=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-dx=" ))
		{
			dx=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, epsilon=self-propulsion parameter, rho0=average density, h=external field strenght, LX*LY=size of the box.
	const double D=1., J=1., gamma=1., r=1.;
	double beta=1., epsilon=0.9, rho0=3., h=0.3;
	int LX=512, LY=512;
	
	//Numerical parameters: init=initial condition, tmax=maximal time.
	double dt=0.001, tmax=5000., dx=0.25, rhog=0.9, rhol=3;
	int init=1, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,epsilon,rho0,h,LX,LY,init,rhog,rhol,dx,dt,tmax,THREAD_NUM);
	
	//OpenMP.
	const int OMP_MAX_THREADS=omp_get_max_threads();
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;
	
	const int NX=int(LX/dx), NY=int(LY/dx);
	
	//Global variables.
	const double Dpara=D*(1+epsilon/3.)*dt/square(dx), Dperp=D*(1-epsilon/3.)*dt/square(dx), v0=2*D*epsilon*dt/(3*dx), gamma0=gamma*dt;

	//Density and magnetization on each sites.
	vector<vector<double> > RHO1(NX,vector<double>(NY,0)), RHO2(NX,vector<double>(NY,0)), RHO3(NX,vector<double>(NY,0)), RHO4(NX,vector<double>(NY,0));
	const vector<vector<double> > HH={{h,0,0,0},{0,0,h,0}};
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_RFAPM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_RFAPM_averages/RFAPM_averages_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	cout.precision(8);
	
	const int Nsteps=int(tmax/dt), DeltaT=int(1./dt);
	meshInit(RHO1,RHO2,RHO3,RHO4,init,rhog,rhol,rho0,beta,r,NX,NY);
	
	//Time evolution.
	for(int t=0;t<=Nsteps;t++)
	{
		//Export data.
		if (t%DeltaT==0 or t==Nsteps)
		{
			exportProfilesAlongX(RHO1,RHO2,RHO3,RHO4,beta,epsilon,rho0,h,LX,LY,NX,NY,init,t*dt);
			exportDensityState(RHO1,RHO2,RHO3,RHO4,beta,epsilon,rho0,h,LX,LY,NX,NY,init,t*dt);
			
			const double n1=average(RHO1,NX,NY), n2=average(RHO2,NX,NY), n3=average(RHO3,NX,NY), n4=average(RHO4,NX,NY);
			cout << "time=" << t*dt << " -rho=" << n1+n2+n3+n4 << " -n1=" << n1 << " -n2=" << n2 << " -n3=" << n3 << " -n4=" << n4 << running_time.TimeRun(" ") << endl;
			fileAverages <<  t*dt << " " << n1+n2+n3+n4 << " " << n1 << " " << n2 << " " << n3 << " " << n4 << endl;
		}
		
		//At each time-step update RHO and MAG arrays.
		finiteDiff(RHO1,RHO2,RHO3,RHO4,Dpara,Dperp,v0,gamma0,beta,HH,r,NX,NY);
	}
	return 0;
}
