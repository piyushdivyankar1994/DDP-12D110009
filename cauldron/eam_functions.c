#include<stdio.h>
#include<stdlib.h>
#include <gsl/gsl_rng.h>
#include<math.h>
#include "math_functions.h"
#include "eam_functions.h"
#include "utils.h"


double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)
{
    if(field_check_eam(index_field, value_field)==0)return 0;
    int pos=binary_search(ptr,5000*index_field,5000*index_field+4999,index_field_value);    //binary search for position of value in array
    pos=pos-index_field*5000;                                                               //pos value needs to be from 0-4999 to make things easy
    return linear_interpolator(ptr[5000*index_field+pos],ptr[5000*value_field+pos],ptr[5000*index_field+pos+1],ptr[5000*value_field+pos+1],index_field_value);
}
//error check function to see if eam_interpolation_func() has been used correctly
int field_check_eam(int index_field,int value_field)
{
    int ret_val=1;

    if(!(index_field != RADIUS || index_field!=ELECTRON_DENSITY_INDEX)){printf("ERROR:incorrect eam index_feild,returning 0\n");ret_val= 0;}
    if(!((value_field>=PAIRING_POTENTIAL_Al_Al && value_field<=ELECTRON_DENSITY_Ni) || value_field!=EMBEDDING_FUNCTION_Ni || value_field!=EMBEDDING_FUNCTION_Al)){printf("ERROR:incorrect eam value feild returning 0\n");ret_val=0;}
    if(index_field == value_field){printf("ERROR:same eam value and eam index field,returning 0\n");ret_val=0;}
    return ret_val;
}

//function to read a particular type of parameter file
//in this first line contains no of parameters that it will contain
//then next lines shall be filled with these parameters and they can only be real numbers
int read_prameter_file(char *filename,double ret_val[])
{
    FILE *fp;
    fp=fopen(filename,"r");
    int no_of_parameters;

    fscanf(fp,"%d ",&no_of_parameters);
    //ret_val=(double*)malloc(no_of_parameters*sizeof(double));
    for(int i=0;i<no_of_parameters;i++)
    {
        fscanf(fp,"%le ",&ret_val[i]);
    }
    return no_of_parameters;
}

//function that takes reads neighbour atom coordinates
void neighbour_lattice_sites_read(double* sites)
{
    FILE *fp;
    int a[]={12,6,24,12,24,8,48};
    for(int i=0;i<7;i++)
    {
        int count;
        char fn[10];
        sprintf(fn,"S%dn.mat",i+1);
        fp=fopen(fn,"r");
        if(fp==NULL)printf("unable to open file S%dn.mat\n",i+1);
        for(int j=0;j<a[i]*3 && count<TOTAL_NEIGHBOUR_ATOMS_FCC;j++,count++)
        {
            fscanf(fp,"%le",&sites[count]);
        }
    }
}

//function to read eam data
//file needs to be passed an array in which the function will store data and name of the file(char* file) which contains names of all files containing neighbour atom coordinates
void eam_data_read(double *eam_data,char *file)
{
    FILE *fp_init,*fp1;
    fp_init=fopen(file,"r");
    for(int i=0;i<7;i++)
    {
       char str[25];
       fscanf(fp_init,"%s\n",str);
       fp1=fopen(str,"r");
       double val1,val2;
       if(fp1==NULL)printf("flag Unable to open %s\n",str);

       if(i==0){
          for(int j=0;j<5000;j++)
          {
              fscanf(fp1,"%le %le ",&val1,&val2);
              eam_data[j+5000]=val2;
              eam_data[j]=val1;
          }
          continue;
       }
       else if(i<5){
           for(int j=0;j<5000;j++)
           {
               fscanf(fp1,"%le %le ",&val1,&val2);
              eam_data[5000*(i+1)+j]=val2;
           }
           continue;
       }
       else if(i==5){
           for(int j=0;j<5000;j++)
           {
               fscanf(fp1,"%le %le ",&val1,&val2);
              eam_data[30000+j]=val1;
              eam_data[35000+j]=val2;
           }
           continue;
           }
       else if(i==6){
           for(int j=0;j<5000;j++)
           {
               fscanf(fp1,"%le %le ",&val1,&val2);
              eam_data[40000+j]=val2;
           }
       }
    }
}
//Montecarlo simulation function this takes as input eam data table , neighbour atom locations(sites), no of montecarlo steps, parameter file name, crystal data input file name, crystal data output file name
void eam_monte_carlo_simulation(double *eam_data,double *sites,int *neighbour_lattice_sites_number,int upto_n_neighbours,int N_MCS, char *parameter_file_name, char *input_file_name,char *output_file_name)
{
    FILE *f_parameter;                         //parmeter file pointer
    f_parameter=fopen(parameter_file_name,"r");//opening parameter file
    int no_of_parameters;                      //for taking first value from file to know size of array we need to store all the data
    double *parameter;

    fscanf(f_parameter,"%d ",&no_of_parameters);  //reading
    parameter=malloc(no_of_parameters*sizeof(double*));  //allocating memory
    read_prameter_file("parameters.txt",parameter);       //reading function called

    double Nx=parameter[0];//size in x direction
    double Ny=parameter[1];//---- -- y ---------
    double Nz=parameter[2];//---- -- z ---------

    int *atomic_matrix;             //array to store atomic species and location
    int total_no_of_atoms=Nx*Ny*Nz*4;//total no of atoms calculated
    atomic_matrix=malloc(4*Nx*Ny*Nz*sizeof(int*));//allocating memory

    FILE *fp_input;         //file for crystal data input
    //reading input data
    fp_input=fopen(input_file_name,"r");
    for(int i=0;i<total_no_of_atoms;i++)
        fscanf(fp_input,"%d",&atomic_matrix[i]);

    //Following initializes the random number generator
    //We are using gsl random number generator for GNU Public Library
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng * k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    k = gsl_rng_alloc (T);

    double lattice_parameter=parameter[6];      //lattice parameter of the Al-Ni lattice in question
    //begin montecarlo steps
    for(int monte_carlo_steps=0;monte_carlo_steps<N_MCS;monte_carlo_steps++)
    {
        double u=gsl_rng_uniform(r);//random number for site selection gives result in (0,1)
        int init_pos=u*total_no_of_atoms;//scaling the result to number total number of atoms
        //A little about how we store the fcc crystal
        //Every FCC lattice can be thought of a simple cubic lattice with a four atom motiff
        //Motiff atoms have following positions
        //(0.0 0.0 0.0),(0.0 0.5 0.5),(0.5 0.0 0.5),(0.5 0.5 0.0)
        //We store the species information in a linear array such that atoms associated with each lattice points are stored sucessively
        //the motiff atoms are needs to acessed uniquely so a coressponds to which atom we are taking about
        //a=0 , 1 , 2 , 3 represent the 4 types of location that an atom can have

        int motiff_position=init_pos%4;          //type of atom
        int atom_x_index=(init_pos/4)%(int)Nx;                     //x-index
        int atom_y_index=(init_pos/(4*(int)Nx))%(int)Ny;           //y-index
        int atom_z_index=(init_pos/(4*(int)Nx*(int)Ny))%(int)Nz;   //z-index


        double x_atom,y_atom,z_atom;                                //represent coordinates of the site in question in global frame
                                                        //x_atom,yc,zc are the miller indcies coordinate
        if(motiff_position==0)                                        //if atom at first position(0,0,0) indicies are the global frame coordinates
        {
            x_atom=(double)atom_x_index;
            y_atom=(double)atom_y_index;
            z_atom=(double)atom_z_index;
        }
        else                                            //in other cases
        {
            x_atom=(double)atom_x_index+sites[3*motiff_position-3];                  //the first 9 entries in sites file are the motiff translations
            y_atom=(double)atom_y_index+sites[3*motiff_position-2];                  //NOTE: This information about motiff needs to be in seperate file
            z_atom=(double)atom_z_index+sites[3*motiff_position-1];
        }

        //0 if aluminium, 1 if Nickel
        int A=atomic_matrix[init_pos];                          //variable to store information on main site
        int neighbour_atom;                                  //variable to contain information about neighbour site
        //NOTE: A,site_atom are confusing variable names to fixed, reponded 22/06/16, change to neighbour atom

        double xi,yi,zi,dxi,dyi,dzi;         //xi,yi,zi are the coordinates of the i'th lattice site in gloabl frame, and dxi,dyi,dzi are variables used to compare fractional parts of their global coordinates
        double r_new,r_old=0;
        int as,xs,ys,zs;                       //array index coordinates for neighbour sites

        double electron_density_at_atom=0;           //stores electron density without changing atom

        double EAM_energy_initial=0;            //stores Energy calculated using without changing atom
        double EAM_energy_switched=0;           //stores Energy calculated using with atom switched
        double pairwise_energy_d=0;             //stores pairwise elctrostatic term of EAM table without changing atom
        double pairwise_energy_s=0;             //stores pairwise elctrostatic term of EAM table atom switched
        double energy_change=0;                 //stores change in energy for the switch
        double pair_e_Ni_Ni, pair_e_Ni_Al, pair_e_Al_Al, eDen_Ni, eDen_Al;
        int neighbour_counter=0;
        for(int i=0;i<402;i+=3)
        {
            xi=x_atom+sites[i];       //finding nighbour coordinates using sites table
            yi=y_atom+sites[i+1];
            zi=z_atom+sites[i+2];


            //Implementing periodic boundary condition
            //if it were simple cubic lattice we would have gobal coordinates go from 0 to Nx-1 in x-direction 0-Ny-1 in y-direction so on.
            //so SCC with our motiff translation may get out of range of array
            //maximum legal value in any direction is Nx-0.5,Ny-0.5 etc.
            //so xi can only be in [0,Nx-0.5] similarly for y and z
            //following 6 lines make sure these conditions are satisfied
            if(xi<0)xi=xi+Nx;
            if(xi>Nx-0.5)xi=xi-Nx;
            if(yi<0)yi=yi+Ny+0.5;
            if(yi>Ny-0.5)yi=yi-Ny;
            if(zi<0)xi=zi+Nz;
            if(zi>Nz-0.5)zi=zi-Nz;

            dxi=xi-floor(xi);//fractional part xi
            dyi=yi-floor(yi);//---------
            dzi=zi-floor(zi);//-----------
            if(dxi==0 && dzi==0 && dyi==0) as=0;//if all 0 then as=1
            else{
                for(int j=0;j<3;j++)//NOTE: Here also we are taking motiff translations from the site matrix, this needs to be changed as mentioned above
                {
                    if((dxi-sites[3*j]) ==0 &&(dyi-sites[3*j+1]) ==0 &&(dzi-sites[3*j+2]) ==0)as=j+1;
                }
            }
            xs=(int)(xi-dxi);//int type casting to make it feasible to pass in array and retrive site information
            ys=(int)(yi-dyi);
            zs=(int)(zi-dzi);
            int pos=as+4*xs+4*(int)Nx*ys+4*(int)Nx*(int)Ny*zs;
            neighbour_atom=atomic_matrix[pos];

            //NOTE:Wondering if indexing tasks should be written as function,
//double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)


/*
*follwing are the index fields and value feilds
*Values in EAM data field are stored in a linear array, searching is only possible i for index_fields, and interpolation only for value_fields
*entering anything other than written below will cause function to return 0
*index_feild=0---> RADIUS
*value_feild=1---> PAIRING POTENTIAL Al-Al
*value_feild=2---> PAIRING POTENTIAL Ni-Al
*value_feild=3---> PAIRING POTENTIAL Ni-Ni
*value_feild=4---> ELECTRON DENSITY Al
*value_feild=5---> ELCETRON DENSITY Ni
*index_feild=6---> ELECTRON DENSITY INDEX FOR EMBEDDING FUNCTION
*value_feild=7---> EMBEDDING FUNCTION Ni
*value_feild=8---> EMBEDDING FUNCTION Al
*anything other than this is considered illegal and function
 */

            r_new=dist_3d(x_atom,y_atom,z_atom,xi,yi,zi)*lattice_parameter; //distance calculation from neighbour site to atomic site
            if(r_new!=r_old) {
              pair_e_Ni_Ni=eam_data_interpolation_func(eam_data,RADIUS,PAIRING_POTENTIAL_Ni_Ni,r_new);
              pair_e_Ni_Al=eam_data_interpolation_func(eam_data,RADIUS,PAIRING_POTENTIAL_Ni_Al,r_new);
              pair_e_Al_Al=eam_data_interpolation_func(eam_data,RADIUS,PAIRING_POTENTIAL_Al_Al,r_new);
              eDen_Ni=eam_data_interpolation_func(eam_data,RADIUS,ELECTRON_DENSITY_Ni,r_new);
              eDen_Al=eam_data_interpolation_func(eam_data,RADIUS,ELECTRON_DENSITY_Al,r_new);
            }
            if(neighbour_atom*A == 4)//If both neighbour and site are Nickel atoms
            {
                pairwise_energy_d=pairwise_energy_d+pair_e_Ni_Ni;
                pairwise_energy_s=pairwise_energy_s+pair_e_Ni_Al;
                electron_density_at_atom=electron_density_at_atom+eDen_Ni;
            }
            else if(neighbour_atom*A == 1)//If both neighbour and site are Aluminum atoms
            {
                pairwise_energy_d=pairwise_energy_d+pair_e_Al_Al;
                pairwise_energy_s=pairwise_energy_s+pair_e_Ni_Al;
                electron_density_at_atom=electron_density_at_atom+eDen_Al;
            }
            else if(neighbour_atom==Al && A==Ni)//If neighbours nickel and site is Aluminum atoms
            {
                pairwise_energy_d=pairwise_energy_d+pair_e_Ni_Al;
                pairwise_energy_s=pairwise_energy_s+pair_e_Ni_Ni;
                electron_density_at_atom=electron_density_at_atom+eDen_Ni;
            }
             else if(neighbour_atom==Ni && A==Al)//If neighbours Aluminium and site is Nickel atoms
            {
                pairwise_energy_d=pairwise_energy_d+pair_e_Ni_Al;
                pairwise_energy_s=pairwise_energy_s+pair_e_Al_Al;
                electron_density_at_atom=electron_density_at_atom+eDen_Al;
            }
        }
        //Now depending on whether the current atom is Ni or aluminium will change the EAM energy calculation result
        if(neighbour_atom==1){//if Ni atom at site
            EAM_energy_initial=pairwise_energy_d+eam_data_interpolation_func(eam_data,ELECTRON_DENSITY_INDEX ,EMBEDDING_FUNCTION_Ni,electron_density_at_atom);//Energy with atom being Ni
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,ELECTRON_DENSITY_INDEX ,EMBEDDING_FUNCTION_Al,electron_density_at_atom);//Energy with atom being Al at the same location
        }
        else{//Same things but for aluminium
            EAM_energy_initial=pairwise_energy_d+eam_data_interpolation_func(eam_data,ELECTRON_DENSITY_INDEX ,EMBEDDING_FUNCTION_Al,electron_density_at_atom);
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,ELECTRON_DENSITY_INDEX ,EMBEDDING_FUNCTION_Ni,electron_density_at_atom);
        }
        energy_change=EAM_energy_switched-EAM_energy_initial;//energy change required for transition probability

        double chance=gsl_rng_uniform(k);//a random number
        double energy_norm_const=(BOLTZMANN_CONST*NORMAL_TEMPERATURE)/ELECTRONIC_CHARGE;//normalization of energy using erature and electronic charge TEMP is the highest temperature stored as a #define constant
        double probablity=exp(-energy_change/energy_norm_const);//probability calculation using custom exponentiation function
        if(probablity>1){atomic_matrix[init_pos]=!atomic_matrix[init_pos];}             //checking if energy was negetive in that case there is sure probability,switching atoms, '!'(logical NOT) used as the numbers are 1 and 0 so it works like a charm
        else if(chance<probablity){atomic_matrix[init_pos]=!atomic_matrix[init_pos];}   //Metrpolis Implementaion-If probability is a number form 0 to 1 say p, and the probablity that gsl_rng_uniform() is a uniform random genrator so the probality that it generates a number from (0,p) is p. This is our stochastic in this step there might be a situation when energy change was positive but transition still took place.

        r_old=r_new;
    }
    file_write(output_file_name,atomic_matrix,total_no_of_atoms);
}

void random_crystal_generator(char *parameter_file_name, char *output_file_name) {
    const gsl_rng_type * T;
    gsl_rng * r;
    FILE *fp,*fw;
    double *parameter;
    fw=fopen(output_file_name,"w");
    char *ch;
    ch=parameter_file_name;
    fp=fopen(ch,"r");
    int no_of_parameters;

    fscanf(fp,"%d ",&no_of_parameters);
    parameter=malloc(no_of_parameters*sizeof(double*));
    printf("%s\n",output_file_name);
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    read_prameter_file(ch,parameter);
    double u=parameter[0]* parameter[1]* parameter[2];

    double conc_g_p=parameter[5]/100;
    int lim_l=(int)parameter[2]*(100-conc_g_p)/200;
    int lim_u=(int)parameter[2]*(100+conc_g_p)/200;
    for(int i=0;i<u;i++)
    {
        //double Nx=i%(int)parameter[0];                      unused variables may  be used in the future to develop a more general function
        //double Ny=(i/(int)parameter[0])%(int)parameter[1];
        double Nz=(i/((int)parameter[0]*(int)parameter[1]))%(int)parameter[2];
       //NOTE Nickel Atoms reprented by 1
       //NOTE Aluminium atons represented by 0
        if(Nz>=lim_l && Nz<=lim_u){
            fprintf(fw,"0 1 1 1\n");
        }
        else{
            for(int s=0;s<4;s++)
            {
                double k=gsl_rng_uniform(r);
                if(k<=conc_g_p)fprintf(fw,"1 ");
                else fprintf(fw,"0 ");

            }
            fprintf(fw,"\n");
        }

    }

    gsl_rng_free(r);
    fclose(fp);
}
