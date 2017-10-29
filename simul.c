#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define ABS(a) ((a)<0 ? -(a): (a))
#define PI 3.14159265
#define NMAX 20
#define MAXNODESIZE 400

typedef struct nodetype *NODEptr;
typedef struct nodetype {
  int k; // integer identifier
  int i; // integer identifier x-coodinate
  int j; // integer identifier y-coordinate
  double x; 
  double y;
  int occupied;
  NODEptr link[6]; // starting at 3 o'clock
}NODEtype;


NODEtype *makeLATTICE(int N, double a) {
				int k;
				int i,j;
				double h = sqrt(3)/2.0*a;
				NODEtype *ptr;  
				ptr = (NODEtype *) malloc(N*N*sizeof(NODEtype));
				for (k=0;k<N*N;k++) {
				ptr[k].k = k;
				ptr[k].i = k % N;
				ptr[k].j = k / N;
				i = k % N;
				j = k / N;    
				ptr[k].x = (double) i*a + j*0.5*a;
				ptr[k].y = (double) j*h;
				ptr[k].occupied = -1;

				if (i==0) {
					ptr[k].link[3] = NULL;
				} else {
					ptr[k].link[3] = &ptr[k-1];
				}

				if (i==N-1) {
					ptr[k].link[0] = NULL;
				} else {
					ptr[k].link[0] = &ptr[k+1];
				}

				if (j==0) {
					ptr[k].link[2] = NULL;
				} else {
					ptr[k].link[2] = &ptr[k-N];
				}
    
				if (j==N-1) {
					ptr[k].link[5] = NULL;
				} else {
					ptr[k].link[5] = &ptr[k+N];
				}

				if (i==N-1 || j == 0) {
					ptr[k].link[1] = NULL;
				} else {
					ptr[k].link[1] = &ptr[k-N+1];
				}

				if (i==0 || j == N-1) {
					ptr[k].link[4] = NULL;
				} else {
					ptr[k].link[4] = &ptr[k+N-1];
				}

				} 

				return(ptr);
}

int occupyLattice( NODEtype *ptr,int N){
	int Np;//no of particles to initialise.
	int numV=N*N; //number of nodes.
	int j,k;
	//randomly select no of points to initialise.
	Np=floor(numV*drand48()); 
	//randomly arrange particles. 
	j=0;
	while(j<Np){
			k=floor((numV)*drand48());
			if (ptr[k].occupied==-1){
			   (ptr[k].occupied)=1;			  
			   ++j;				   
			   printf("k is : %2d  occupancy is :%2d \n",k,ptr[k].occupied);	  
			}			
    }
   
    return(Np);
}//end occupyLattice

//Funtion doGillespie performs Gillespie algorithm on Lattice.
//N is the length of lattice.
//Np is the no of particles on lattice.
//NODEtype is a struct of the lattice
//randNp is the random particle to move.
void doGillespie(NODEtype *ptr,int N,int Np,int *Nflux){
		void gillMove(NODEtype*,int,int,int*); 
		double rActual, rTotal=0;
		double kmove,dt,t;
		//int *Nflux;
		int nf=0,randNp;
		Nflux=&nf;
		while(nf<10){//change this to *Nflux for periodic BC
			rTotal=6*kmove*Np;//total rate is based on filled positions.
			rActual = rTotal*drand48(); 
			dt=-1.0/rTotal*log(drand48());
			//particle to move?		
			// in this case since they are all the same kind of reaction 
			//we can compute the position satisfying the Gillespie criterion thus:
			randNp=ceil(rActual/6/kmove);
			//perform monte carlo movements.
			gillMove(ptr,randNp,N,Nflux);	
			t=t+dt;
		}//end while 
}//end doGillespie.

//Function gillMove performs the actual movement of particles on the Lattice sites 
//it also implements the necessary boundary conditions if necessary.
//N is the length of lattice.
//Np is the no of particles on lattice.
//NODEtype is a struct of the lattice
void gillMove(NODEtype *ptr, int randNp,int N,int* Nflux){
		int j;
		int count=0;
		int go;
		int k;
		for (j=0;j<=(N*N)-1;j++){
			if(ptr[j].occupied==1){
			count+=1;	   
			}
			if ( count==randNp-1 || j==(N*N)-1  )
			break; 
			//j now holds position of the occupied site that we need.
		}
		//printf("rand Np is :%d and j is : %d  occup is : %d \n",randNp,j,ptr[j].occupied);
		//printf("occupancy at j is : %d  \n",ptr[j].occupied);
		//attempt to move particle at node j;
		go=floor(6*drand48());
		//attempt a move starting from 3 'o clock.
		if ( (j<N*N)  && (ptr[j].occupied==1)  &&  (ptr[j].link[go]!=NULL)  && (ptr[j].link[go]->occupied==-1) ){// double check it is occupied. 
			//a false is actually a logical error, since j must be occupied.	
			//jump to postion link[go]
			ptr[j].link[go]->occupied=1;
			ptr[j].occupied=-1;		
			//++*Nflux;
	
		}
		//implement periodic B.C if the chosen go position does not exist.
		//it is assumed that the dimension of lattice is N by N else these won't strictly hold
		else if ( (j<N*N)  && (ptr[j].link[go])==NULL  &&   (ptr[j].occupied==1)     ){

				//Left Face Periodic BC		
				//The left hand bottom corner node.//region 1 a
				if ( ptr[j].i ==0  &&   ptr[j].j ==0    &&  ptr[j+N-1].occupied==-1){
					ptr[j].occupied =-1;
					ptr[j+N-1].occupied =1;
					++*Nflux;
				}
		
				//left top corner node.//region 2 b
				else if ( ptr[j].i ==0  &&   ptr[j].j ==N-1   && ptr[j+N-1].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j+N-1].occupied =1;
						++*Nflux;
				}		
				//right face top corner node //region 4  c
				else if (ptr[j].i ==N-1   && ptr[j].j ==N-1    && ptr[j-(N-1)].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j-(N-1)].occupied =1;
						++*Nflux;
				}	

				//Right Bottom corner node  periodic BC  //region 5 d
				else if (ptr[j].j ==0 &&   ptr[j].i ==N-1  &&  ptr[j-(N-1)].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j-( N-1 ) ].occupied=1;
						++*Nflux;
				}		
		
				//Right Face Periodic BC exclude corners //region 6 e
				else if (ptr[j].i ==N-1   && ptr[j].j !=N-1   && ptr[j].j !=0 &&     ptr[j-(N-1)].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j-(N-1)].occupied =1;
						++*Nflux;
				}
		
				//Left Face BC.Excluding the corners //region 3 //Labelled as region for ease of identification f
				else if ( ptr[j].i ==0  &&   ptr[j].j !=0  &&   ptr[j].j !=N-1  &&  ptr[j+N-1].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j+N-1].occupied =1;
						++*Nflux;
				}
				
				//Bottom Face periodic BC Excluding the corners //region 7  g
				else if (ptr[j].j ==0   && ptr[j].i !=0  && ptr[j].i !=N-1   && ptr[j+( N*( N-1 ) )].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j+( N*( N-1 ) )].occupied=1;
						++*Nflux;
				}
				
				//Top Face Periodic BC Excluding the corners //region 8
				else if (ptr[j].j ==N-1 && ptr[j].i !=N-1   && ptr[j].i !=0  && ptr[j-( N*( N-1 ) )].occupied==-1){
						ptr[j].occupied =-1;
						ptr[j-( N*( N-1 ) )].occupied =1;
						++*Nflux;
				}		
		
		}// end of periodic BC branches
}//end function gillMove

int main() { 
			int Rseed;
			Rseed = (unsigned) time(NULL);
			srand48(Rseed); 
			int nflux=0;
			int *Nflux=&nflux;
			int k;
			int N = 7;
			int Np,nc=0;
			int i,m;
			NODEtype *ptr;  
			ptr= makeLATTICE(N,1.0);
			Np=occupyLattice(ptr,N);
			printf("no of particle is %d\n",Np);     
			printf(" results after occupying sites: \n");
			for (k=0;k<N*N;k++) {
				printf("%2d ",(ptr[k].occupied+1)/2);//transform -1 to 0 and 1 to 1
				if ( ptr[k].i ==N-1 ){
					printf("\n");
				}
			}
  
			if (Np>1)
				doGillespie(ptr,N,Np,Nflux);
			else printf("Lattice should have at least 2 particles \n");  
			printf(" results after implementing Gillespie Algorithm: \n");
			for (k=0;k<N*N;k++) {
				printf("%2d ",((ptr[k].occupied)+1)/2);//retransformed values.
				if ( ptr[k].i ==N-1 ){
					printf("\n");
				}
			}  
			printf("\n");
			printf("printing adjacency matrix \n");
			for (k=0;k<N*N;k++) {	 
				for (i=0;i<N*N;i++) {		   
						for (m=0;m<6;m++) {					
								if ( ptr[k].link[m] !=NULL && ptr[k].link[m]->k == ptr[i].k ){
									printf(" %d ",1);
									nc++;
									break;												
								}							
						}
						if (nc==0)
						printf(" %d ",0);
						nc=0;		 
				}
				printf("\n"); 
			} 
		printf("\n ");
		printf("printing x and y co-ordinates \n ");	
		printf("x         y  \n ");
		printf("----------------\n ");	
		for (k=0;k<N*N;k++) {
            printf("|%4.3f |  %4.3f|  \n ",ptr[k].x,ptr[k].y);			
			printf("----------------\n ");		
		} 
		printf("\n");     
		return(0);
		
}

