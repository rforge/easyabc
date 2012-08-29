/************************************************************************************
    Massive simulations of trait-based filtered communities for ABC estimation
    Franck Jabot, 2nd September 2008.

Last modified on 31th december 2008.

*************************************************************************************/

// Libraries
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <string>
#include <math.h>
#include <malloc.h>
#include <list>
using namespace std;

//Random number generator
double genrand2(void);
void sgenrand2(unsigned long);
unsigned long genrand2i(void);
void sgenrand2i(unsigned long);

double calculfitness(double *tabfitness,int S,double **trait,double *h,double *sig,double SS,double ntrait){
    double min=SS;
    double max=0.0;
    for (int j=0;j<ntrait;j++){
        for (int i=0;i<S;i++){
            tabfitness[i]*=exp(log(SS)/ntrait)*exp(-(trait[i][j]-h[j])*(trait[i][j]-h[j])/(2*sig[j]*sig[j]));
        	//cerr<<trait[i]<<" "<<h<<" "<<sig<<" "<<SS<<" "<<(tabfitness[i])<<endl;
        }
    }
    for (int i=0;i<S;i++){
        if (tabfitness[i]>max){
            max=tabfitness[i];
        }
        if (tabfitness[i]<min){
            min=tabfitness[i];
        }
    }
    for (int i=0;i<S;i++){
        tabfitness[i]=((1+tabfitness[i])/(1+min)-1);
    }
    
return (((1.0+max)/(1.0+min))-1.0);
}

int draw(double *abondreg){
	double u=genrand2();
	int k=0;
	while (abondreg[k]<u){
		u-=abondreg[k];
		k++;
	}
return k;
}

int draw2(double *abondreg,double *sumabondfit,double *tabfitness){
	double u=sumabondfit[1]*genrand2();
	int k=0;
	while ((abondreg[k]*(1+tabfitness[k]))<u){
		u-=(abondreg[k]*(1+tabfitness[k]));
		k++;
	}
return k;
}

int drawint (int *abondloc,double *sumabondfit,double *tabfitness){
	double u=sumabondfit[0]*genrand2();
	int k=0;
	while ((abondloc[k]*(1+tabfitness[k]))<u){
		u-=(abondloc[k]*(1+tabfitness[k]));
		k++;
	}
return k;
}

int drawint2 (int *abondloc,int J){
	int u=genrand2i()%J;
	int k=0;
	while (abondloc[k]<(1+u)){
		u-=abondloc[k];
		k++;
	}
return k;
}


void initialiser(int *abondloc,double *abondreg,int J,double *sumabondfit,double *tabfitness){
	for (int i=0;i<J;i++){
		int k=draw2(abondreg,sumabondfit,tabfitness);
		abondloc[k]++;
		sumabondfit[0]+=(1+tabfitness[k]);
	}
}

void stepdyn(int *abondloc,double *tabfitness,double *abondreg,double maxSS,double I,int J,double *sumabondfit){
	int k=drawint2(abondloc,J);
	abondloc[k]--;
	//cerr<<k<<" "<<abondloc[k]<<" "<<sumabondfit[0]<<" ";
	sumabondfit[0]-=(1+tabfitness[k]);
	//cerr<<sumabondfit[0]<<endl;
	//int toto;
	//cin >> toto;
	double u=genrand2();
	if (u<(I/(I+J-1))){
		//int kk=draw(abondreg);
		int kk=draw2(abondreg,sumabondfit,tabfitness);
		abondloc[kk]++;
		sumabondfit[0]+=(1+tabfitness[kk]);
	}
	else{
		int kk=drawint(abondloc,sumabondfit,tabfitness);
		abondloc[kk]++;
		sumabondfit[0]+=(1+tabfitness[kk]);
	}
}

void forwarddynamics(int *abondloc,double *tabfitness,double *abondreg,int S,int J,double I,double maxSS,int lsimul){
	for (int i=0;i<S;i++){
		abondloc[i]=0;
	}
	double *sumabondfit;
	sumabondfit=new double[2];
	sumabondfit[0]=0.0;
	sumabondfit[1]=0.0;
	for (int i=0;i<S;i++){
		sumabondfit[1]+=(abondreg[i]*(1+tabfitness[i]));
	}
	initialiser(abondloc,abondreg,J,sumabondfit,tabfitness);
	for (int i=0;i<(J*lsimul);i++){
		stepdyn(abondloc,tabfitness,abondreg,maxSS,I,J,sumabondfit);
	}
	//for (int i=0;i<S;i++){
	//	cerr<<abondloc[i]<<" ";
	//}
	//int toto;
	//cin >> toto;
	
}

void calculstat(double *stat,int *abondloc, int l, double **trait,double ntrait){
    //STATISTIQUES CALCULEES:
    //J
    //S
    //Shan
    //Var(Ni)
    //Mean(trait)
    //Var(trait)
    //Skewness(trait)
    //Mean(trait)_espece
    //Var(trait)_espece
    //Skewness(trait)_espece
    //Fourth Moment
    //Fourth Moment espece
    
    for (int i=0;i<(4+ntrait);i++){
        stat[i]=0;
    }
    
    for (int i=0;i<l;i++){
        stat[0]+=abondloc[i];
        if (abondloc[i]>0){
            stat[1]+=1;
            stat[2]-=abondloc[i]*log(abondloc[i]);
        }
    }
    stat[2]+=stat[0]*log(stat[0]);
    stat[2]/=stat[0];
    double meanN=stat[0]/stat[1];
    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            stat[3]+=(abondloc[i]-meanN)*(abondloc[i]-meanN);
        }
    }
    stat[3]/=stat[1];
    if (stat[1]>1){
        stat[3]*=(stat[1]/(stat[1]-1));
    }

 for (int k=0;k<ntrait;k++){

    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            stat[(4+8*k)]+=abondloc[i]*trait[i][k];
            stat[(7+8*k)]+=trait[i][k];
        }
    }
    stat[(4+8*k)]/=stat[0];
    stat[(7+8*k)]/=stat[1];
    
    for (int i=0;i<l;i++){
        if (abondloc[i]>0){
            stat[(5+8*k)]+=abondloc[i]*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)]);
            stat[(8+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
            stat[(6+8*k)]+=abondloc[i]*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)]);
            stat[(9+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
	    stat[(10+8*k)]+=abondloc[i]*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)])*(trait[i][k]-stat[(4+8*k)]);
	    stat[(11+8*k)]+=(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)])*(trait[i][k]-stat[(7+8*k)]);
        }
    }
    stat[(5+8*k)]/=stat[0];
    stat[(6+8*k)]/=stat[0];
    stat[(8+8*k)]/=stat[1];
    stat[(9+8*k)]/=stat[1];
    stat[(10+8*k)]/=stat[0];
    stat[(11+8*k)]/=stat[1];
    
    //rendre les estimateurs non biaisés :
    if (stat[1]>1){
        stat[(8+8*k)]*=(stat[1]/(stat[1]-1));
    }
    stat[(5+8*k)]*=(stat[0]/(stat[0]-1));
  }
}



int main(){


//LECTURE DES FICHIERS D'ENTREE 
    char buffer[256];
    ifstream in("input2");
	int nb_simul;
	in >> nb_simul;
	unsigned long *seed1tab;
	seed1tab= new unsigned long[nb_simul];
	unsigned long seed1;
	int *Jtab;
	Jtab= new int[nb_simul];
	int J; 
	double *Itab;
	Itab= new double[nb_simul];
	double I;
	double *SStab;
	SStab= new double[nb_simul];
	double SS;
	double *ntraittab;
	ntraittab= new double[nb_simul];
	double ntrait;
	double **htab; 
	htab= new double*[nb_simul];
	double **sigtab; 
	sigtab= new double*[nb_simul];
for (int is=0;is<nb_simul;is++){
	//in >> buffer; 
	in >> seed1tab[is];
    //in >> buffer; 
	in >> Jtab[is];
    //in >> buffer; int nbsimul; in >> nbsimul;
    //in >> buffer; 
	in >> Itab[is];
    //in >> buffer; double Imax; in >> Imax;
    //in >> buffer;
	in >> SStab[is];
    //in >> buffer; double SSmax; in >> SSmax;
    //in >> buffer; 
	in >> ntraittab[is];
    //in >> buffer; 
	htab[is]= new double[int(ntraittab[0])];
    for (int i=0;i<ntraittab[0];i++){
        in >> htab[is][i];
    }
    //in >> buffer; double *hmax; hmax= new double[int(ntrait)];
    // for (int i=0;i<ntrait;i++){
    //    in >> hmax[i];
    //}
    //in >> buffer; 
	sigtab[is]= new double[int(ntraittab[0])];
    for (int i=0;i<ntraittab[0];i++){
        in >> sigtab[is][i];
    }
    //in >> buffer; double *sigmax; sigmax= new double[int(ntrait)];
    //for (int i=0;i<ntrait;i++){
    //    in >> sigmax[i];
    //}
}
    in.close();
    int lsimul;
    
    ifstream in2("species_traits");
    in2 >> buffer; int S; in2 >> S;
    ntrait=ntraittab[0];
    double *h; h= new double[int(ntrait)];
	double *sig; sig= new double[int(ntrait)];
	double **trait;
    trait= new double*[S];
    for (int i=0;i<S;i++){
        trait[i]=new double[int(ntrait)];
    }
    double *abondreg;
    abondreg=new double[S];
    double sumabondreg=0.0;
    for (int i=0;i<S;i++){
        in2 >> abondreg[i];
	sumabondreg+=abondreg[i];
        for (int j=0;j<ntrait;j++){
            in2 >> trait[i][j];
        }
    }
    in2.close();
    
	for (int i=0;i<S;i++){
        	abondreg[i]/=sumabondreg;
    	}
	
	//for (int i=0;i<S;i++){
	//	cerr<<abondreg[i]<<" "<<trait[i]<<" "<<abondobs[i]<<endl;
	//}
	//int toto;
	//cin>>toto;


//PREPARATION DES FLUX DE SORTIE
    //sprintf(buffer,"ABCpt_Filteredcom2_J_%d_nbsimul_%d_m_%d_h_%d_sig_%d_ss_%d.txt",J,nbsimul,int(floor(Imax)),int(floor(hmax)),int(floor(sigmax)),int(floor(SSmax)));
    ofstream out("output2");
    //out<<"I\t SS\t maxSS\t h\t sig\t J\t S\t Shan\t VarNi\t"
    //for (int i=0;i<ntrait;i++){
    //    out<<"Meantrait"<<(i+1)<<"\t Vartrait"<<(i+1)<<"\t Moment3trait"<<(i+1)<<"\t Meantraitespece"<<(i+1)<<"\t Vartraitespece"<<(i+1)<<"\t Moment3traitespece"<<(i+1)<<"\t Moment4trait"<<(i+1)<<"\t Moment4traitespece"<<(i+1)<<"\t";
    //}
    //out<<endl;  


 	int *abondloc;
	abondloc= new int[S];
	for (int i=0;i<S;i++){
		abondloc[i]=0;
	}
	double maxSS;
	double *tabfitness;
	tabfitness=new double[S];
	for (int i=0;i<S;i++){
        	tabfitness[i]=1;
    }
	ntrait=ntraittab[0];
	double *stat;
	int nstat=4+ntrait*8;
	stat= new double[nstat];
	lsimul=Jtab[0];

    int i=0;
	while(i<nb_simul){
		seed1=seed1tab[i];
		sgenrand2(seed1*100000);
		sgenrand2i(1027+seed1*10000);
		J=Jtab[i];
		I=Itab[i];
		SS=SStab[i];
		for (int ii=0;ii<ntrait;ii++){
			h[ii]=htab[i][ii];
		}
		for (int ii=0;ii<ntrait;ii++){
			sig[ii]=sigtab[i][ii];
		}
		//DRAWING OF PARAMETERS
		I=exp(I);
		SS=exp(SS);
		for (int j=0;j<ntrait;j++){   
        //    h[j]=hmin[j]+(hmax[j]-hmin[j])*genrand2();
		      sig[j]=exp(sig[j]);
        }
		maxSS=calculfitness(tabfitness,S,trait,h,sig,SS,ntrait);


		//SIMULATION
		forwarddynamics(abondloc,tabfitness,abondreg,S,J,I,maxSS,lsimul);
		
		//COMPUTATION OF STATISTICS
		calculstat(stat,abondloc,S,trait,ntrait);


		//OUTPUT
		//if ((stat[1]>Starget-11)&&(stat[1]<Starget+11)){
        		//out<<I<<" "<<SS<<" "<<maxSS<<" ";
				//for (int k=0;k<ntrait;k++){
				//	out<<h[k]<<" "<<sig[k]<<" ";
				//}
        		for (int j=1;j<nstat;j++){
				    out<<stat[j]<<"\t";
			    }
			    out<<endl;
                i++;
		//}
		//if ((i%50)==49)	{
			//cerr<<(i+1)<<" ";
		//}
	}
	return 1;
}


/**** GENERATEUR DE NOMBRES ALEATOIRES ****/
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

//#include<stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void sgenrand2(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double genrand2()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}


/* initializing the array with a NONZERO seed */
void sgenrand2i(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

unsigned long genrand2i()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2i(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y; 
}
