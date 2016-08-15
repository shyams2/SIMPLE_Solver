#include<iostream>
#include<curses.h>
#include<fstream>
#include<math.h>

using namespace std;

//solver version 1.0
//For the sake of simplicity, the variables will be declared as global so all functions may access it
//In future updates, to improve the efficiency they functions will be passed by reference to subfunctions
//Future updates to improve the solver efficiency and capabilites

double u[90601],v[90601],p[90000];
//velocities and pressure

double psix[90601],psiy[90601],psi[90601];
//for stream function calculation

double vorticity[90000];
//stores vorticity values

double uav[90000],vav[90000];
//stores averages

double dudx[90000],dudy[90000],dvdx[90000],dvdy[90000];
//gradients

double m[90000],pc[90000];
//mass and pressure corrections

double dudt,dvdt,dx,dy,dt;
//parameters for solution and time gradients

double rho,nu,con;
//constants

double resmax,psimax;
//max values

int n,k;
//grid count

int n1,n2,n3,n4,n5;
//used for nodal connectivity

//max capacity of 90000 elements for the grid
//dynamic memory allocation to be added

void intro()
{

	cout<<"================================================="<<endl;
	cout<<"|	    Welcome to SIMPLE_Solver	        |"<<endl;
	cout<<"| This code solver solves for lid driven cavity |"<<endl;
	cout<<"|       This solver is under development        |"<<endl;
	cout<<"|	    				        |"<<endl;
	cout<<"|	    				        |"<<endl;
	cout<<"|       Developer:Shyam Sundar Sankaran	        |"<<endl;
	cout<<"================================================="<<endl;
	cout<<"\n13*****14*****15*****16\n";
	cout<<"*		      *"<<endl;
	cout<<"* VII  * VIII  *  IX  *"<<endl;
	cout<<"* 		      *"<<endl;
	cout<<"9*****10******11*****12"<<endl;
	cout<<"*		      *"<<endl;
	cout<<"*  IV *   V   *   VI  *"<<endl;
	cout<<"*		      *"<<endl;
	cout<<"5******6*******7******8"<<endl;
	cout<<"*		      *"<<endl;
	cout<<"*  I   *   II  *  III *"<<endl;
	cout<<"*		      *"<<endl;
	cout<<"1******2*******3******4"<<endl<<endl<<endl;
	cout<<"The above cells indicate how we refer to velocities"<<endl;
}

void initialize()
{
	int i;
	for(i=0;i<(n*n);i++)
		u[i]=v[i]=0;
		//initializing initially all points with zero velocity
	for(i=(n*k);i<(n*n);i++)
	{
		u[i]=1;
		//top points of the container have a velocity of 1
		psix[i]=psiy[i]=0;
	}
	for(i=0;i<n;i++)
		psix[i]=psiy[i]=0;
	for(i=0;i<(n*n);i=i+n)
		psix[i]=psiy[i]=0;
	for(i=n;i<(n*n);i=i+n)
		psix[i]=psiy[i]=0;
		//the above section of code assigns all points on the boundary with a psi value of zero
	for(i=0;i<(k*k);i++)
		p[i]=0;
		//assigning initial guess pressure as zero to all the cells
}

void nodcon(int i,int n)
{
	n1=i+(int)(i/n);
	n2=n1+1;
	n3=n1+n+2;
	n4=n1+n+1;
	n5=n3+1+(int)(i/n);
	//used to generate nodal connectivity
}

void calculate()
{
	int i;
	resmax=m[0];//initial assigning to the first variable, we are to check for maximum
	for(i=0;i<(k*k);i++)
	{
		nodcon(i,k);
		uav[i]=(u[n1]+u[n2]+u[n3]+u[n4])/4;
		vav[i]=(v[n1]+v[n2]+v[n3]+v[n4])/4;
		dudx[i]=(u[n2]+u[n3]-u[n1]-u[n4])/(2*dx);
		dvdx[i]=(v[n2]+v[n3]-v[n1]-v[n4])/(2*dx);
		dudy[i]=(u[n3]+u[n4]-u[n1]-u[n2])/(2*dy);
		dvdy[i]=(v[n3]+v[n4]-v[n1]-v[n2])/(2*dy);
		m[i]=(dudx[i]+dvdy[i])*dx*dy;
		pc[i]=-(con)*m[i];//the constant may be changed appropriately
		if(m[i]>resmax)
		resmax=m[i];
	}
}

void massbal()//this function performs mass balance for mass balance cells
{
	int i;
	for(i=0;i<(k*k);i++)
	{
		if(i<k&&(i%k)==0)
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i]+pc[i+k]+pc[i]);
		else if(i<k&&(i%k)==(k-1))
			pc[i]=0.8*pc[i]+0.05*(pc[i]+pc[i-1]+pc[i+k]+pc[i]);
		else if(i<k)
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i-1]+pc[i+k]+pc[i]);
                else if(i%k==0)
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i]+pc[i+k]+pc[i-k]);
		else if(i%k==(k-1))
			pc[i]=0.8*pc[i]+0.05*(pc[i]+pc[i-1]+pc[i+k]+pc[i-k]);
		else if(i>((k*(k-1))-1)&&(i%k)==0)
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i]+pc[i]+pc[i-k]);
		else if(i>((k*(k-1))-1)&&(i%k)==(k-1))
			pc[i]=0.8*pc[i]+0.05*(pc[i]+pc[i-1]+pc[i]+pc[i-k]);
		else if(i>((k*(k-1))-1))
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i-1]+pc[i]+pc[i-k]);
		else
			pc[i]=0.8*pc[i]+0.05*(pc[i+1]+pc[i-1]+pc[i+k]+pc[i-k]);
		//the above section of code attempts to stabilize the solution obtained, as we are using an approximate solver
		p[i]=p[i]+pc[i];
	}
}

void correction()//this function performs the velocity correction for each of the nodes
{
	int i;
	for(i=0;i<((k-1)*(k-1));i++)
	{
		nodcon(i,(k-1));
		u[n5]=u[n5]-((dt/rho)*((pc[n2]+pc[n3]-pc[n1]-pc[n4])/(2*dx)));
		v[n5]=v[n5]-((dt/rho)*((pc[n4]+pc[n3]-pc[n1]-pc[n2])/(2*dy)));
	}
}

void mombal()//this function performs momentum balance to the momentum balance cells
{
	//here we apply that (storage of momentum)+(mass efflux)=(forces due to shear and normal pressures)
	int i;
	for(i=0;i<((k-1)*(k-1));i++)
	{
		nodcon(i,(k-1));
		double mn,ms,me,mw,un,us,ue,uw,vn,vs,ve,vw,xnormalf,xshearf,xmflux,ynormalf,yshearf,ymflux;
		mn=0.5*rho*(vav[n3]+vav[n4])*dx;
		ms=0.5*rho*(vav[n1]+vav[n2])*dx;
		me=0.5*rho*(uav[n2]+uav[n3])*dy;
		mw=0.5*rho*(uav[n1]+uav[n4])*dy;
		un=0.5*(uav[n3]+uav[n4]);
		us=0.5*(uav[n1]+uav[n2]);
		ue=0.5*(uav[n2]+uav[n3]);
		uw=0.5*(uav[n1]+uav[n4]);
		vn=0.5*(vav[n3]+vav[n4]);
		vs=0.5*(vav[n1]+vav[n2]);
		ve=0.5*(vav[n2]+vav[n3]);
		vw=0.5*(vav[n1]+vav[n4]);
		//the following section applies upwinding
		if((double)(1/(k*nu))>2)// this condition ensures that upwinding would be implemented only when	the cell reynolds number is 						   greater than two
		{
			if(uw>0)
			{
				uw=0.5*uw+0.5*u[n5-1];
				vw=0.5*vw+0.5*v[n5-1];
			}
			else
			{
				uw=0.5*uw+0.5*u[n5];
				vw=0.5*vw+0.5*v[n5];
			}
			if(ue>0){
				ue=0.5*ue+0.5*u[n5];
				ve=0.5*ve+0.5*v[n5];
			}
			else
			{
				ue=0.5*ue+0.5*u[n5+1];
				ve=0.5*ve+0.5*v[n5+1];
			}
			if(vn>0)
			{
				vn=0.5*vn+0.5*v[n5];
				un=0.5*un+0.5*u[n5];
			}
			else
			{
				vn=0.5*vn+0.5*v[n5+n];
				un=0.5*un+0.5*u[n5+n];
			}
			if(vs>0)
			{
				vs=0.5*vs+0.5*v[n5-n];
				us=0.5*us+0.5*u[n5-n];
			}
			else
			{
				vs=0.5*vs+0.5*v[n5];
				us=0.5*us+0.5*u[n5];
			}
		}
		xmflux=(mn*un-ms*us+me*ue-mw*uw);
		ymflux=(mn*vn-ms*vs+me*ve-mw*vw);
		xnormalf=(0.5*(p[n1]+p[n4]-p[n3]-p[n2])+nu*(dudx[n2]+dudx[n3]-dudx[n1]-dudx[n4]))*dy;
		ynormalf=(0.5*(p[n1]+p[n2]-p[n3]-p[n4])+nu*(dvdx[n3]+dvdx[n4]-dvdx[n1]-dvdx[n2]))*dx;
		xshearf=0.5*nu*dx*(dudy[n3]+dudy[n4]+dvdx[n3]+dvdx[n4]-dudy[n1]-dudy[n2]-dvdx[n1]-dvdx[n2]);
		yshearf=0.5*nu*dy*(dudy[n2]+dudy[n3]+dvdx[n2]+dvdx[n3]-dudy[n1]-dudy[n4]-dvdx[n1]-dvdx[n4]);
		dudt=(xnormalf+xshearf-xmflux)/(rho*dx*dy);
		dvdt=(ynormalf+yshearf-ymflux)/(rho*dx*dy);
		u[n5]=u[n5]+(dudt*dt);
		v[n5]=v[n5]+(dvdt*dt);
	}
}

void psical()//this funtion calculates psi for each of the grid points
{
	int i,j;
	for(j=1;j<k;j++)
	{
		for(i=1;i<k;i++)
		{
			psiy[(n*j)+i]=psiy[(n*(j-1))+i]+0.5*(u[(n*j)+i]+u[(n*(j-1))+i])*dy;
			psix[(n*j)+i]=psix[(n*j)+(i-1)]-0.5*(v[(n*j)+i]+v[(n*j)+(i-1)])*dx;
		}
	}
	psimax=-1;
	for(i=0;i<(n*n);i++)
	{
		psi[i]=0.5*(psix[i]+psiy[i]);
		if(psi[i]>psimax)
			psimax=psi[i];
	}
}

void stability()//this function assigns the value for dt, looking at general stability criterion and correction constant
{
	double trial1,trial2;
	trial1=(double)(dx);
	trial2=(double)(0.5*dx*dx*rho/nu);
	if(trial1<trial2)
		dt=(double)0.01*trial1;
	else
		dt=(double)0.01*trial2;
	con=(double)(0.1*rho/dt);
}

int main()
{
intro();
printf("\nFor a square grid, enter the size to be taken(limit-300)\n");
cin>>k;
n=(k+1);
initialize();
double side;
cout<<"\nEnter the value of size of square cavity concerned\n";
cin>>side;
dx=dy=(double)(side/k);//the size for each element of the grid which will be used for solving
cout<<"\nEnter values of Reynold's number which will be used in the program\n";
rho=1;//since we want to solve using reynolds number, we are taking the value of density as
//unity, so that the value of visocity may be assigned as the reciprocal of entered reynolds number
double rey;
cin>>rey;
nu=(double)(1/rey);
stability();
int i;
long int loopcounter;//as name suggests, this variable will hold the number of outer iterations
printf("\nThe value of dt obtained is %f\n",dt);
getchar();
getchar();
loopcounter=0;
while(1)
{
if(loopcounter%500==0)
{
cout<<"=========================================================================="<<endl;
cout<<"MaxResidue\t\t\tLoopCount"<<endl;
cout<<"=========================================================================="<<endl;
}
double old;//this variable stores the value of psimax a fixed number of iterations before
calculate();
mombal();
while(1)
{
int checker=0;
calculate();
massbal();
correction();
if(loopcounter%500==0&&checker%100==0)
cout<<resmax<<"\t\t\t"<<loopcounter<<endl;
checker++;
if(resmax<=1e-6)
break;
}
psical();
loopcounter++;
if(loopcounter==10000)
old=psimax;
if(loopcounter>10000)
{
if(loopcounter%10000==0)
{
if((int)((10000000)*(old-psimax))==0)//we are checking if the value of psi has stabilized to check for convergence
loopcounter=1000000;
else
old=psimax;
}}
if(loopcounter>=1000000)
break;
}
psical();
ofstream myfile;
myfile.open("psi.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey<<endl;
for(i=0;i<(n*n);i++)
{
myfile<<psi[i]<<" ";
if((i+1)%n==0)
myfile<<endl;
}
myfile.close();
myfile.open("vorticity.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey<<endl;
for(i=0;i<(k*k);i++)
{
calculate();
vorticity[i]=dvdx[i]-dudy[i];
myfile<<vorticity[i]<<" ";
if((i+1)%k==0)
myfile<<endl;
}
myfile.close();
myfile.open("u.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey<<endl;
for(i=0;i<(n*n);i++)
{
myfile<<u[i]<<" ";
if((i+1)%n==0)
myfile<<endl;
}
myfile.close();
myfile.open("v.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey<<endl;
for(i=0;i<(n*n);i++)
{
myfile<<v[i]<<" ";
if((i+1)%n==0)
myfile<<endl;
}
myfile.close();
myfile.open("ucenter.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey;
for(i=(int)(k/2);i<(n*n);i=i+n)
myfile<<u[i]<<" ";
myfile.close();
myfile.open("vcenter.txt");
myfile<<"The value we are printing are for"<<k<<"X"<<k<<"grid, for a reynolds number of"<<rey;
for(i=(n*k*.5);i<(n*k*.5+n);i++)
myfile<<v[i]<<" ";
myfile.close();
return(0);
}
