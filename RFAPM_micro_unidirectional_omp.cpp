/*C++ CODE - MANGEAT MATTHIEU - 2026*/
/*ACTIVE POTTS MODEL UNDER A HOMOGENEOUS UNIDIRECTIONAL FIELD*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <numeric>
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
#include "lib/random_OMP.cpp"
#include "lib/special_functions.cpp"

double modulo(const double &x, const double &LX);
void modulo(int &x, const int &LX);

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

//Position and spin defining the Potts particle
class particle
{
	public:
	
	int x,y;
	int dx,dy;
	int spin;
	
	particle(const int &LX, const int &LY, const int &init);
	int band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX);
	void move(const double &r, const double &epsilon, const int &LX, const int &LY);
	void flip(const int &sigma);
};

//Creation of bands
int particle::band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX)
{
	const double r=ran();
	double XX=0;
	if (r<kappa)
	{
		XX=modulo(x0+alpha*r/kappa,1);
	}
	else
	{
		XX=modulo(x0+(1-alpha)*(r-kappa)/(1-kappa)+alpha,1);
		spin=int(4*ran()); //Gaseous phase -> change the spin to random.
	}
	return int(XX*LX);
}

//Creation of the spin
particle::particle(const int &LX, const int &LY, const int &init)
{
	static const double alpha=0.25, kappa=0.8;
	dx=0; dy=0;
	if (init==0) //Gas phase.
	{
		spin=int(4*ran());
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else if (init==1) //Tranverse band orthogonal to the field.
	{
		spin=1;
		x=int(ran()*LX);
		y=band_formation(0.5,alpha,kappa,LY);
		//y=int(3*LY/8+ran()*LY/4);
	}
	else if (init==2) //Longitudinal band orthogonal to the field.
	{
		spin=1;
		x=band_formation(0.5,alpha,kappa,LX);
		//x=int(3*LX/8+ran()*LX/4);
		y=int(ran()*LY);
	}
	else if (init==3) //Liquid phase orthogonal to the field.
	{
		spin=1;
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else if (init==4) //Tranverse band along the field.
	{
		spin=0;
		x=band_formation(0.5,alpha,kappa,LY);
		//x=int(3*LY/8+ran()*LY/4);
		y=int(ran()*LX);
	}
	else if (init==5) //Longitudinal band along the field.
	{
		spin=0;
		x=int(ran()*LX);
		y=band_formation(0.5,alpha,kappa,LX);
		//y=int(3*LX/8+ran()*LX/4);
	}
	else if (init==6) //Liquid phase along the field.
	{
		spin=0;
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
}

//Update the particle position.
void particle::move(const double &r, const double &epsilon, const int &LX, const int &LY)
{
	int p=spin; //Probability to hop in the favored direction (r<epsilon/3).
	if (r>epsilon/3.) //Probability to hop to a random direction (r>epsilon/3).
	{
		p=int(4*ran()); //Random initial direction (0,1,2,3)
	}
	
	if (p==0) //Move right.
	{
		x++;
		dx++;
	}
	else if (p==1) //Move up.
	{
		y++;
		dy++;
	}
	else if (p==2) //Move left.
	{
		x--;
		dx--;
	}
	else //Move down.
	{
		y--;
		dy--;
	}
	modulo(x,LX);
	modulo(y,LY);
}

//Flip the spin.
void particle::flip(const int &sigma)
{
	spin=sigma;
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
double modulo(const double &x, const double &LX)
{
	if (x<0)
	{
		return x+LX;
	}
	else if (x>=LX)
	{
		return x-LX;
	}
	else
	{
		return x;
	}
}

void modulo(int &x, const int &LX)
{
	if (x<0)
	{
		x+=LX;
	}
	else if (x>=LX)
	{
		x-=LX;
	}
}

//Maximum state
int max_state(const vector<int> &rhos)
{
	int rhomax=rhos[0];
	for (int j=1;j<rhos.size();j++)
	{
		if (rhos[j]>rhomax)
		{
			rhomax=rhos[j];
		}
	}
	
	if (rhomax==0)
	{
		cerr << "THE MAX DENSITY IS ZERO -> ERROR" << endl;
		abort();
	}
	
	vector<int> kmax;
	for (int j=0;j<rhos.size();j++)
	{
		if (rhos[j]==rhomax)
		{
			kmax.push_back(j);
		}
	}
	
	return kmax[int(kmax.size()*ran())];
}

//MSD
double msd(const vector<particle> &POTTS, const int &Npart)
{
	double dx2=0, dy2=0;
	for (int i=0; i<Npart; i++)
	{
		dx2+=square(POTTS[i].dx)/Npart;
		dy2+=square(POTTS[i].dy)/Npart;
	}
	return dx2+dy2;
}

//Total average on all space.
double average(const vector< vector<int> > &RHO, const int &LX, const int &LY)
{
	int rhoAv=0;
	for (int x0=0; x0<LX; x0++)
	{
		for (int y0=0; y0<LY; y0++)
		{
			rhoAv+=RHO[x0][y0];
		}
	}
	return double(rhoAv)/(LX*LY);
}

//Average along the y-axis.
vector<double> average1d_along_y(const vector< vector<int> > &RHO, const int &LX, const int &LY)
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
vector<double> average1d_along_x(const vector< vector<int> > &RHO, const int &LX, const int &LY)
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

//Export averaged densities along y-axis.
void exportProfilesAlongX(const vector<vector<vector<int> > > &RHOS, const double &beta, const double &h, const double &D, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	const vector<double> r0=average1d_along_y(RHOS[0],LX,LY), r1=average1d_along_y(RHOS[1],LX,LY), r2=average1d_along_y(RHOS[2],LX,LY), r3=average1d_along_y(RHOS[3],LX,LY);
	
	int returnSystem=system("mkdir -p data_APM_dynamics1d/");
	stringstream ss;
	ss << "./data_APM_dynamics1d/APM_profile_along_x_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameProfile = ss.str();
	
	ofstream fileProfile(nameProfile.c_str(),ios::trunc);
	fileProfile.precision(6);
	fileProfile << "#site\tRHO\tRHO0\tRHO1\tRHO2\tRHO3" << endl;
	
	for (int x0=0; x0<LX; x0++)
	{
		fileProfile << x0 << "\t" << r0[x0]+r1[x0]+r2[x0]+r3[x0] << "\t" << r0[x0] << "\t" << r1[x0] << "\t" << r2[x0] << "\t" << r3[x0] << endl;
	}
	fileProfile.close();
}


//Export averaged densities along x-axis.
void exportProfilesAlongY(const vector<vector<vector<int> > > &RHOS, const double &beta, const double &h, const double &D, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	const vector<double> r0=average1d_along_x(RHOS[0],LX,LY), r1=average1d_along_x(RHOS[1],LX,LY), r2=average1d_along_x(RHOS[2],LX,LY), r3=average1d_along_x(RHOS[3],LX,LY);
	
	int returnSystem=system("mkdir -p data_APM_dynamics1d/");
	stringstream ss;
	ss << "./data_APM_dynamics1d/APM_profile_along_y_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameProfile = ss.str();
	
	ofstream fileProfile(nameProfile.c_str(),ios::trunc);
	fileProfile.precision(6);
	fileProfile << "#site\tRHO\tRHO0\tRHO1\tRHO2\tRHO3" << endl;
	
	for (int y0=0; y0<LY; y0++)
	{
		fileProfile << y0 << "\t" << r0[y0]+r1[y0]+r2[y0]+r3[y0] << "\t" << r0[y0] << "\t" << r1[y0] << "\t" << r2[y0] << "\t" << r3[y0] << endl;
	}
	fileProfile.close();
}

//Export density.
void exportDensity(const vector< vector<int> > &RHO, const double &beta, const double &h, const double &D, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	vector<uint16_t> rho(LX*LY,0);
	for (int y0=0;y0<LY;y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			rho[y0*LX+x0]=RHO[x0][y0];
		}
	}
	
	//Creation of the file.
	int returnSystem=system("mkdir -p data_APM_dynamics2d/");
	stringstream ssDensity;
	ssDensity << "./data_APM_dynamics2d/APM_density_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameDensity = ssDensity.str();
	ofstream fileDensity(nameDensity.c_str(),ios::binary);
	fileDensity.write(reinterpret_cast<const char*>(rho.data()), rho.size()*sizeof(uint16_t));
	fileDensity.close();
}

//Export state.
void exportState(const vector<vector< vector<int> > > &RHOS, const vector< vector<int> > &RHO, const double &beta, const double &h, const double &D, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	vector<uint8_t> state(LX*LY,0);
	for (int y0=0;y0<LY;y0++)
	{
		for (int x0=0; x0<LX; x0++)
		{
			if (RHO[x0][y0]==0)
			{
				state[y0*LX+x0]=0;
			}
			else
			{
				vector<int> rhos={RHOS[0][x0][y0],RHOS[1][x0][y0],RHOS[2][x0][y0],RHOS[3][x0][y0]};
				state[y0*LX+x0]=max_state(rhos)+1;
			}
		}
	}
	
	//Creation of the file.
	int returnSystem=system("mkdir -p data_APM_dynamics2d/");
	stringstream ssState;
	ssState << "./data_APM_dynamics2d/APM_state_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameState = ssState.str();
	ofstream fileState(nameState.c_str(),ios::binary);
	fileState.write(reinterpret_cast<const char*>(state.data()), state.size()*sizeof(uint8_t));
	fileState.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &h, double &D, double &epsilon, double &rho0, int &LX, int &LY, int &tmax, int &init, int &RAN, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-h=" ))
		{
			h=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-D=" ))
		{
			D=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
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
	cout << "-beta=" << beta << " -h=" << h << " -D=" << D << " -epsilon=" << epsilon << " -rho0=" << rho0 << " -LX=" << LX << " -LY=" << LY << " -tmax=" << tmax << " -init=" << init << " -ran=" << RAN << " -threads=" << THREAD_NUM << endl;
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, h=external field, D=diffusion rate, epsilon=self-propulsion parameter, rho0=average density, LX*LY=size of the box.
	double beta=0.5, h=0.03, D=1., epsilon=0.3, rho0=3.;
	int LX=200, LY=200;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG, THREAD_NUM=number of threads.
	int init=1, tmax=150000, RAN=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,h,D,epsilon,rho0,LX,LY,tmax,init,RAN,THREAD_NUM);

	//OpenMP.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;

	//Start the random number generator.
	init_gsl_ran();
	for (int k=0; k<THREAD_NUM; k++)
	{
		gsl_rng_set(GSL_r[k],THREAD_NUM*RAN+k);
	}
	
	//Number of particles.
	int Npart=int(rho0*LX*LY);
	
	//Density and magnetization on each sites.
	vector<vector<vector<int> > > RHOS(4,vector<vector<int> >(LX,vector<int>(LY,0)));
	vector<vector<int> > RHO(LX,vector<int>(LY,0));
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_APM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_APM_averages/APM_averages_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Creation of particles.
	vector<particle> POTTS;	
	for (int i=0;i<Npart;i++)
	{
		particle Potts(LX,LY,init);
		POTTS.push_back(Potts);
		RHOS[Potts.spin][Potts.x][Potts.y]++;
		RHO[Potts.x][Potts.y]++;
	}
	
	//Time increment.
	const double dt=1./(4*D+exp(4*beta*(1+h)));
	
	//Field component (along alpha=0)
	const vector<double> H={h,0,0,0};
	
	//Probability to hop.
	const double proba_hop=4*D*dt;
	
	//Number of particles per threads
	const int Nth=Npart/THREAD_NUM;
	
	//Lock for the update of magnetization
	vector< vector<omp_lock_t> > lock_site(LX,vector<omp_lock_t>(LY));
	
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			omp_init_lock(&lock_site[x][y]);
		}
	}
	
	vector<double> MAG2(0), MSD(0);
	
	//Time evolution.
	for(int t=0;t<=tmax;t++)
	{
		//Export data.
		const double ntot=average(RHO,LX,LY), p1=average(RHOS[0],LX,LY)/ntot, p2=average(RHOS[1],LX,LY)/ntot, p3=average(RHOS[2],LX,LY)/ntot, p4=average(RHOS[3],LX,LY)/ntot;
		const double m1=(4*p1-1)/3., m2=(4*p2-1)/3., m3=(4*p3-1)/3., m4=(4*p4-1)/3.;
		
		const double mag2=sqrt(m1*m1+m2*m2+m3*m3+m4*m4);
		MAG2.push_back(mag2);
		
		const double DR2=msd(POTTS,Npart);
		MSD.push_back(DR2);
		
		if (t%1000==0 or t==tmax)
		{
			fileAverages <<  t << " " << ntot << " " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << mag2 << " " << DR2 << endl;
			cout << "time=" << t << " -rho=" << ntot << " -p1=" << p1 << " -p2=" << p2 << " -p3=" << p3 << " -p4=" << p4 << " -mag2=" << mag2 << " -msd=" << DR2 << running_time.TimeRun(" ") << endl;
		}
		
		if (t%150==0 or t==tmax)
		{
			exportDensity(RHO,beta,h,D,epsilon,rho0,LX,LY,init,RAN,t);
			exportState(RHOS,RHO,beta,h,D,epsilon,rho0,LX,LY,init,RAN,t);
		}
		
		if (t%1000==0 or t==tmax)
		{
			exportProfilesAlongX(RHOS,beta,h,D,epsilon,rho0,LX,LY,init,RAN,t);
			exportProfilesAlongY(RHOS,beta,h,D,epsilon,rho0,LX,LY,init,RAN,t);
		}
		
		//At each time-step move all particles randomly.
		#pragma omp parallel
		{
			const int actual_thread=omp_get_thread_num();
			for (int i=0; i<Nth; i++)
			{
				//Choose a particle randomly (j), at the site (x0,y0).
				const int j=Nth*(actual_thread+ran());
				const int x0=POTTS[j].x, y0=POTTS[j].y, spin0=POTTS[j].spin;
				
				//New possible spin-state
				int spin1=int(3*ran());
				if (spin1>=spin0)
				{
					spin1++;
				}
				
				//Probability to flip.
				omp_set_lock(&lock_site[x0][y0]);
				const int rhos=RHOS[spin0][x0][y0], rhosp=RHOS[spin1][x0][y0], rho=RHO[x0][y0];
				omp_unset_lock(&lock_site[x0][y0]);
				
				const double DH=H[spin0]-H[spin1];
				const double proba_flip=exp(-4*beta*double(rhos-rhosp-1)/double(rho)-4*beta*DH)*dt;
				
				if (proba_hop+proba_flip>1)
				{
					cerr << "THE PROBABILITY TO WAIT IS NEGATIVE: " << 1-(proba_hop+proba_flip) << endl;
					abort();
				}
				
				const double random_number=ran();
				//The particle hops: perform the hopping on the particle and modify the density/magneization (const spin).
				if (random_number<proba_hop)
				{
					POTTS[j].move(random_number/proba_hop,epsilon,LX,LY);
					
					//Update density and magnetization on old site.
					omp_set_lock(&lock_site[x0][y0]);
					RHOS[spin0][x0][y0]--;
					RHO[x0][y0]--;
					omp_unset_lock(&lock_site[x0][y0]);
					
					//Update density and magnetization on new site.
					const int xN=POTTS[j].x, yN=POTTS[j].y;
					omp_set_lock(&lock_site[xN][yN]);
					RHOS[spin0][xN][yN]++;
					RHO[xN][yN]++;
					omp_unset_lock(&lock_site[xN][yN]);
				}
				//The particle flips: perform the flipping on the particle and modify the magnetization (on-site, ind. of epsilon).
				else if (random_number<proba_hop+proba_flip)
				{
					POTTS[j].flip(spin1);
					
					//Update the magnetization.
					omp_set_lock(&lock_site[x0][y0]);
					RHOS[spin0][x0][y0]--;
					RHOS[spin1][x0][y0]++;
					omp_unset_lock(&lock_site[x0][y0]);
				}
				//Else do nothing (proba_wait)
			}			
		}		
	}
	
	//Creation of the density file.
	int returnSystem=system("mkdir -p data_APM_averages_mag/");
	
	stringstream ssMAG2;
	ssMAG2 << "./data_APM_averages_mag/APM_averages_mag2_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameMAG2 = ssMAG2.str();
	ofstream fileMAG2(nameMAG2.c_str(),ios::binary);
	fileMAG2.write(reinterpret_cast<const char*>(MAG2.data()), MAG2.size()*sizeof(double));
	fileMAG2.close();
	
	returnSystem=system("mkdir -p data_APM_averages_msd/");
	
	stringstream ssMSD;
	ssMSD << "./data_APM_averages_msd/APM_averages_msd_beta=" << beta << "_h=" << h << "_D=" << D << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameMSD = ssMSD.str();
	ofstream fileMSD(nameMSD.c_str(),ios::binary);
	fileMSD.write(reinterpret_cast<const char*>(MSD.data()), MSD.size()*sizeof(double));
	fileMSD.close();
	
	return 0;
}
