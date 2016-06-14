#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//This piece of code calculates isotropic lennard jones energy
//pos[] takes in an array for positions labelled as (x1,y1,x2,y2,x3,y3,...,xN,yN)
//N is the total number of points (xi,yi) int the array
//box size is needed for the periodic boundary condition implementation

double lennard_jones_energy(double pos[],int N,double box_size)
{
	int i,j;
	double energy,dx,dy,r;
	
	//for total energy all the points are taken pairwise and totol energy is summed 
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
				if(dx>0.5*box_size)dx=1-dx;
				if(dy>0.5*box_size)dy=1-dy;
				r=dx*dx+dy*dy;

				//printf("%lf %lf %lf %lf %d \n",pos[2*j-1],pos[2*j],dx,dy,j);
				double invr6=1/pow(r,3);
				energy=energy+invr6*(invr6-2);
			}
		}
	}
	return energy;

}
