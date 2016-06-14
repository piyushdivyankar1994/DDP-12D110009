#include<stdio.h>
#include<stdlib.h>
#include<libenergyfunctions.h>

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
	energy=lennard_jones_energy(pos,10,1.0);
	printf("%lf\n",energy);
}
