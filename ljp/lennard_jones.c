#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double lennard_jones_energy(double pos[],int N);

int main(void)
{
	FILE *fp1,*fp2,*fp3;
	fp1=fopen("positions.dat","r");

	//reading positions into a local array

	double pos[20];
	int i,j,k;
	double energy=0;
	double dx,dy,r;

	for(i=0;i<20;i++)
	{
		fscanf(fp1,"%lf",&pos[i]);
	}

	//calculating lennardjones energy
	energy=lennard_jones_energy(pos,10);
	printf("%lf\n",energy);
}

double lennard_jones_energy(double pos[],int N)
{
	int i,j;
	double energy,dx,dy,r;
	
	for(i=0;i<N;i++)
	{
		printf("%lf %lf --\n",pos[2*i],pos[2*i+1]);
		for(j=0;j<N;j++)
		{
			//we want all the pairs of atoms to add to energy
			if(i!=j)
			{
				//double r=(pos[2*i-1]-pos[2*j-1])*(pos[2*i-1]-pos[2*j-1])+	(pos[2*i]-pos[2*j])*(pos[2*i]-pos[2*j]);
				dx=fabs(pos[2*i]-pos[2*j]);
				dy=fabs(pos[2*i+1]-pos[2*j+1]);
				if(dx>0.5)dx=1-dx;
				if(dy>0.5)dy=1-dy;
				r=dx*dx+dy*dy;

				//printf("%lf %lf %lf %lf %d \n",pos[2*j-1],pos[2*j],dx,dy,j);
				double invr6=1/pow(r,3);
				energy=energy+invr6*(invr6-2);
			}
		}
	}
	return energy;

}
