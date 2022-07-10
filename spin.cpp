#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<cmath>


using namespace std;

#define l 64
#define KB 1
//#define T 4
//#define ImageWidth 1000  //image width
//#define ImageHeight 1000 //image height

int pbc(int vsitio, int vrede){
	if(vsitio < 0){return (vsitio+vrede);}
	if(vsitio >= vrede){return (vsitio-vrede);}
	else{return vsitio;}
}//periodic boundarys conditions

double Energ(int i, int j, int *n){
	 double energia = -n[pbc(i+j*l,l*l)] * ( n[pbc((i+1)+j*l,l*l)] + n[pbc((i-1)+j*l,l*l)] + n[pbc(i+(j+1)*l,l*l)] + n[pbc(i+(j-1)*l,l*l)]);	
	
	return energia;
}//calcula a Energia de um sitio com os seus vizinho mais proximos

double deltaE(int i, int j, int *n){
	
	double dE = 2 * n[pbc(i+j*l,l*l)]*(n[pbc((i+1)+j*l,l*l)]+n[pbc((i-1)+j*l,l*l)]+n[pbc(i+(j+1)*l,l*l)]+n[pbc(i+(j-1)*l,l*l)]);
	return dE;
}//calcula a variacao de Energia


double Magn(int *n){
	
	double mag = 0;
	for(int i = 0; i < l*l; i++){ mag += n[i];}
	
	return mag;
}


int main(){

	double tc = 2.0/(log(sqrt(2)+1));
	double nu = 1.0;
	double beta = 1.0/8.0;
	double gama = 7.0/8.0;
		
	int size = l*l;
	int n[l*l];
	int plai, plaj;
	int nr = int((l*l)/5);//#de steps
	
	double E = 0;
	double E_ = 0;
	double E2 = 0;
	double E2_ = 0;
	double varE = 0;
	double sorte = 0;

	int maxam = 10000;

	srand48(660799);//+28*time(0));

	for(int i = 0; i < (l*l); i++){
		n[i] = 1;		
	}
	

	for(int i = 0; i < l; i++){
		for(int j = 0; j < l; j++){
			E_ = Energ(i,j,n);//pra calcular a Etotal
			E += E_/2;
			//cout << E_ << endl;	
		}		
	}
	//cout << "\t" << E << endl;
	
	
	
	
	
	for(double T=0.5;T<5;T += 0.05){
        	       	
    	double mediaE = 0;
		double mediaM = 0;		
		
		double mediaE2 = 0;
		double mediaM2 = 0;

		double mediaM4 = 0;

		double SpecHeat = 0;
		double MagnSus = 0;
		
		double g_L = 0;
				
		for(int amostras = 0; amostras < maxam; amostras++){
		                
		        
		    //MCsweap
		    for(int MCsweap = 0; MCsweap <= 21; MCsweap++){

		        //MCstep
		        for(int MCstep = 0; MCstep < nr ; MCstep++){
		            plai = drand48()*int(l);
		            plaj = drand48()*int(l);
		            
		            varE  = deltaE(plai, plaj,n);
		            //cout << "\t\t" << varE << endl;
		            
		            if(varE < 0){
		                n[plai+plaj*l] = -n[plai+plaj*l]; 
		                E = E + varE; //actualizar E
		                //cout << "deltaE" << deltaE(plai,plaj,n) << endl;
		            }
		            if(varE > 0){
		                sorte = drand48();
		                if(sorte < exp(-varE/(KB*T))){
		                    n[plai+plaj*l] = -n[plai+plaj*l];
		                    E = E + varE; //actualizar E
		                    //cout << "deltaE " << varE << endl;
		                }
		            }
		        }	

		    
		    //cout << T << "\t" << Magn(n) << endl;

		            
		    }//end MCsweap
				    
		    double magnn = Magn(n);
		    //double magnn_ = abs(magnn/double(size));
		    double magnn_ = magnn/double(size);
			
		    //quantidades por spin		
		    mediaE += E;
		    mediaM += magnn;
		
		    mediaE2 += (E*E);                
		    mediaM2 += (magnn*magnn);         

		    mediaM4 += (magnn*magnn*magnn*magnn);

		    //cout << T << "\t" << E << "\t" << magnn << endl;
		    //cout << (magnn*magnn*magnn*magnn)<< "\t" << magnn_ << endl;				
		
		}//end Amostras

		SpecHeat = ( mediaE2/double(maxam) - (mediaE/double(maxam))*(mediaE/double(maxam)) )/double(T*T);
		SpecHeat = SpecHeat/double(size);
		
		MagnSus = ( mediaM2/double(maxam) - (mediaM/double(maxam))*(mediaM/double(maxam)) )/double(T);
		MagnSus = MagnSus/double(size);
		
		g_L = 0.5*(3 - ((mediaM4/double(maxam))/((mediaM2/double(maxam))*(mediaM2/double(maxam)))));
		
		//cout << T << "\t" << mediaE//double(size*maxam) << "\t" << mediaM//double(size*maxam) << "\t" << SpecHeat << "\t" << MagnSus << endl;
		
		//cout << T << "\t" << double(mediaE/maxam) << "\t" << double(mediaE2/maxam) << "\t" <<  SpecHeat << endl;
		//cout << T << "\t" << double(mediaM/maxam)*double(mediaM/maxam) << "\t" << double(mediaM2/maxam) << "\t" <<  MagnSus << endl;
		
		//cout << T << "\t" << SpecHeat << "\t" << MagnSus << endl;
		
		//cout << pow(l,1.0/nu)*(T - tc) << "\t" << double(mediaM/maxam)*pow(l,beta/nu) << "\t" << MagnSus/double(pow(l,(gama/nu))) << endl;
		
		//cout << pow(l,1.0/nu)*(T - tc) << "\t" << MagnSus << "\t" << MagnSus/double(pow(l,(gama/nu))) << endl;
		//algo de errado, o finite size scaling esta a funcionat pra magnsus*l^gamma/nu 
		//emvez de funcionar pra o magnsus/l^gamma/nu1
		
		//cout << (mediaM4/double(maxam)) << "\t" << (mediaM2/double(maxam)) << endl;
		//cout << T << "\t" << g_L << endl;
		double finite_nu = 1.0;
		cout << pow(l,double(1.0/finite_nu))*(T - 2.275) << "\t" << g_L << endl;

	}//end do for das T
	




}//end main
