#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#define size 402

#define TEMP 2000

/*
 *follwing are the index fields and value feilds
*index_feild=0--->RADIUS
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
//function declarations
int field_check_eam(int ,int);
int binary_search(double*, int, int, double);
double linear_interpolator(double, double, double,double,double);
double eam_data_interpolation_func(double* , int , int , double );
double dist_3d(double,double,double,double,double,double);
void eam_monte_carlo_simulation(double *,double *,int, char *, char *,char *);
void neighbour_lattice_sites_read(double* );
int read_parameter(char*,double*);
void file_write(char *,int *,int );
double eam_data_interpolation_func(double *,int, int, double);

//Following function takes in the EAM table array and radius or electron density value and give coressponding potential or e- density or embedding energy values linearly interpolated between the tow values the index_field value is located
double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)
{
    if(field_check_eam(index_field, value_field)==0)return 0;
    int pos=binary_search(ptr,5000*index_field,5000*index_field+4999,index_field_value);    //binary search for position of value in array
    pos=pos-index_field*5000;                                                               //pos value needs to be from 0-4999 to make things easy
    return linear_interpolator(ptr[5000*index_field+pos],ptr[5000*value_field+pos],ptr[5000*index_field+pos+1],ptr[5000*value_field+pos+1],index_field_value);
}

//Interpolates value of a variable given two points of a straight line
double linear_interpolator(double x1,double y1,double x2,double y2,double xc)
{
    double m=(y1-y2)/(x1-x2);   //slope
    return y1+m*(xc-x1);        
}

//Function to calculate euclidean distance between 2 points in 3D
double dist_3d(double x1,double y1,double z1,double x2,double y2,double z2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

//exponentian function so we avoid using the default one ans save time
double exponential(double x)
{
    double sum=1;
    for(int i=20;i>0;--i)
        sum=1+x*sum/i;
}

//binary search function, searches the position of value that is lower limit of the smallest interval
int binary_search(double *ptr,int i,int j,double val)
{
    if((j-i)>1)
    {
        int mid=(j+i)/2;
        if(ptr[mid]==val)return mid;
        else if(ptr[mid]>val)binary_search(ptr,i,mid,val);
        else binary_search(ptr,mid,j,val);
    }
    else
    {
        return i;
    }
}

//error check function to see if eam_interpolation_func() has been used correctly
int field_check_eam(int index_field,int value_field)
{
    int ret_val=1;

    if(!(index_field != 0 || index_field!=6)){printf("ERROR:incorrect eam index_feild,returning 0\n");ret_val= 0;}
    if(!((value_field>=1 && value_field<=5) || value_field!=7 || value_field!=8)){printf("ERROR:incorrect eam value feild returning 0\n");ret_val=0;}
    if(index_field == value_field){printf("ERROR:same eam value and eam index field,returning 0\n");ret_val=0;}
    return ret_val;
}

//function to read a particular type of parameter file
//in this first line contains no of parameters that it will contain 
//then next lines shall be filled with these parameters and they can only be real numbers
int read_parameter(char *filename,double ret_val[])
{
    FILE *fp;
    fp=fopen(filename,"r");
    int no_of_params;

    fscanf(fp,"%d ",&no_of_params);
    //ret_val=(double*)malloc(no_of_params*sizeof(double));
    for(int i=0;i<no_of_params;i++)
    {
        fscanf(fp,"%le ",&ret_val[i]);
    }
    return no_of_params;
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
        printf("Reading from file %s\n",fn);
        fp=fopen(fn,"r");
        if(fp==NULL)printf("unable to open file S%dn.mat\n",i+1);
        for(int j=0;j<a[i]*3 && count<size;j++,count++)
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
       printf("Reading from file %s\n",str);
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

//main function to wrap all the functions 
int main(void)
{

    double eam_data[45000];
   // int *site[4];
    
    double sites[size];
    int count=0;
    neighbour_lattice_sites_read(sites);

    eam_data_read(eam_data,"file_list.txt");

    eam_monte_carlo_simulation(eam_data,sites,1000, "parameters.txt","input_data.txt","Result_lattice.txt");

    return 0;
}

//Montecarlo simulation function this takes as input eam data table , neighbour atom locations(sites), no of montecarlo steps, parameter file name, crystal data input file name, crystal data output file name
void eam_monte_carlo_simulation(double *eam_data,double *sites,int N_MCS, char *parameter_file_name, char *input_file_name,char *output_file_name)
{
    FILE *f_param;                         //parmeter file pointer
    f_param=fopen(parameter_file_name,"r");//opening parameter file
    int no_of_params;                      //for taking first value from file to know size of array we need to store all the data
    double *parameter;

    fscanf(f_param,"%d ",&no_of_params);  //reading 
    parameter=(double*)malloc(no_of_params*sizeof(double));  //allocating memory
    int n=read_parameter("parameters.txt",parameter);       //reading function called
    
    double Nx=parameter[0];//size in x direction
    double Ny=parameter[1];//---- -- y ---------
    double Nz=parameter[2];//---- -- z ---------
    
    int *atoms;             //array to store atomic species and location
    int no_of_atoms=Nx*Ny*Nz*4;//total no of atoms calculated
    atoms=(int *)malloc(4*Nx*Ny*Nz*sizeof(int));//allocating memory
    
    FILE *fp_input;         //file for crystal data input
    //reading input data
    fp_input=fopen(input_file_name,"r");
    printf("Reading file input_data.txt\n");
    for(int i=0;i<no_of_atoms;i++)
        fscanf(fp_input,"%d",&atoms[i]);

    //Following initializes the random number generator 
    //We are using gsl random number generator for GNU Public Library
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng * k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    k = gsl_rng_alloc (T);

    double lattice_parameter=4.00;      //lattice parameter of the Al-Ni lattice in question
    //begin montecarlo steps
    for(int monte_carlo_steps=0;monte_carlo_steps<N_MCS;monte_carlo_steps++)
    {
        double u=gsl_rng_uniform(r);//random number for site selection gives result in (0,1)
        int init_pos=u*no_of_atoms;//scaling the result to number total number of atoms
        //A little about how we store the fcc crystal 
        //Every FCC lattice can be thought of a simple cubic lattice with a four atom motiff
        //Motiff atoms have following positions 
        //(0.0 0.0 0.0),(0.0 0.5 0.5),(0.5 0.0 0.5),(0.5 0.5 0.0)
        //We store the species information in a linear array such that atoms associated with each lattice points are stored sucessively 
        //the motiff atoms are needs to acessed uniquely so a coressponds to which atom we are taking about
        //a=0 , 1 , 2 , 3 represent the 4 types of location that an atom can have
        
        int a=init_pos%4;          //type of atom
        int x=(init_pos/4)%(int)Nx;                     //x-index
        int y=(init_pos/(4*(int)Nx))%(int)Ny;           //y-index
        int z=(init_pos/(4*(int)Nx*(int)Ny))%(int)Nz;   //z-index


        double xc,yc,zc;                                //represent coordinates of the site in question in global frame
                                                        //xc,yc,zc are the miller indcies coordinate
        if(a==0)                                        //if atom at first position(0,0,0) indicies are the global frame coordinates
        {
            xc=(double)x; 
            yc=(double)y;
            zc=(double)z;
        }
        else                                            //in other cases
        {
            xc=(double)x+sites[3*a-3];                  //the first 9 entries in sites file are the motiff translations
            yc=(double)y+sites[3*a-2];                  //NOTE: This information about motiff needs to be in seperate file 
            zc=(double)x+sites[3*a-1];
        }
        
        //0 if aluminium, 1 if Nickel
        int A=atoms[init_pos];                          //variable to store information on main site                      
        int site_atom;                                  //variable to contain information about neighbour site
        //NOTE: A,site_atom are confusing variable names to fixed
        
        double xi,yi,zi,r,dxi,dyi,dzi;         //xi,yi,zi are the coordinates of the i'th lattice site in gloabl frame, and dxi,dyi,dzi are variables used to compare fractional parts of their global coordinates 
        int as,xs,ys,zs;                       //array index coordinates for neighbour sites

        double electron_density_d=0;           //stores electron density without changing atom 
        double electron_density_s=0;           //----------------------- with atom switched(at xc,yc,zc; Al replaced with Ni, or Ni replaced with Ni
        double EAM_energy_default=0;            //stores Energy calculated using without changing atom
        double EAM_energy_switched=0;           //stores Energy calculated using with atom switched
        double pairwise_energy_d=0;             //stores pairwise elctrostatic term of EAM table without changing atom
        double pairwise_energy_s=0;             //stores pairwise elctrostatic term of EAM table atom switched
        double energy_change=0;                 //stores change in energy for the switch
        for(int i=0;i<402;i+=3)
        {
            xi=xc+sites[i];       //finding nighbour coordinates using sites table
            yi=yc+sites[i+1];
            zi=zc+sites[i+2];

            r=dist_3d(xc,yc,zc,xi,yi,zi)*lattice_parameter; //distance calculation from neighbour site to atomic site
           
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
            site_atom=atoms[pos];

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
            if(site_atom==A && A==1)//If both neighbour and site are Nickel atoms
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,3,r);//we access the Ni-Ni pairwise energy in value_field 3 at radius 'r'
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,2,r);//we switch site so we access Al-Ni potential value_field=2
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,5,r);//this the electron density function we acess value_field 4 and key is still radius this will be further used to evaluate Embedding energy
                //All follwing statements contain same order of statements
            }
            else if(site_atom==A && A==0)//If both neighbour and site are Aluminum atoms
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,1,r);//value_field=1 Al-Al
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,2,r);//value_field=2 Al-Ni
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,4,r);//value_field=4 Ni e- Density
            }
            else if(site_atom==0 && A==1)//If neighbours nickel and site is Aluminum atoms
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,2,r);//value_field=2 Al-Ni
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,3,r);//value_field=3 Ni-Ni
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,5,r);//value field=5 Al e- density
            }
             else if(site_atom==1 && A==0)//If neighbours Aluminium and site is Nickel atoms
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,2,r);
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,1,r);
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,4,r);
            }  
        }
        //Now depending on whether the current atom is Ni or aluminium will change the EAM energy calculation result
        if(site_atom==1){//if Ni atom at site
            EAM_energy_default=pairwise_energy_d+eam_data_interpolation_func(eam_data,6,7,electron_density_d);//Energy with atom being Ni
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,6,8,electron_density_d);//Energy with atom being Al at the same location
        }
        else{//Same things but for aluminium 
            EAM_energy_default=pairwise_energy_d+eam_data_interpolation_func(eam_data,6,8,electron_density_d);
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,6,7,electron_density_d);
        }
        energy_change=EAM_energy_switched-EAM_energy_default;//energy change required for transition probability

        double chance=gsl_rng_uniform(k);//a random number
        double energy_norm_const=(1.38e-23*TEMP)/1.602e-19;//normalization of energy using temperature and electronic charge TEMP is the highest temperature stored as a #define constant
        double probablity=exponential(-energy_change/energy_norm_const);//probability calculation using custom exponentiation function 
        if(probablity>1){atoms[init_pos]=!atoms[init_pos];}             //checking if energy was negetive in that case there is sure probability,switching atoms, '!'(logical NOT) used as the numbers are 1 and 0 so it works like a charm 
        else if(chance<probablity){atoms[init_pos]=!atoms[init_pos];}   //Metrpolis Implementaion-If probability is a number form 0 to 1 say p, and the probablity that gsl_rng_uniform() is a uniform random genrator so the probality that it generates a number from (0,p) is p. This is our stochastic in this step there might be a situation when energy change was positive but transition still took place.
        

    }
    file_write(output_file_name,atoms,no_of_atoms);
}

void file_write(char *filename,int *data,int data_points)
{
    FILE *fp;
    fp=fopen(filename,"w");
    for(int i=0;i<data_points;i+=4)
    {
        fprintf(fp,"%d %d %d %d \n",data[i],data[i+1],data[i+2],data[i+3]);
    }
    printf("%s is the resulting lattice arrangement\n",filename);
}
