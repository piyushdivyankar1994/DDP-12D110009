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
int field_check_eam(int ,int);
int binary_search(double*, int, int, double);
double linear_interpolator(double, double, double,double,double);
double eam_data_interpolation_func(double* , int , int , double );
double dist_3d(double,double,double,double,double,double);

int read_parameter(char*,double*);
void file_write(char *,int *,int );
double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)
{
    if(field_check_eam(index_field, value_field)==0)return 0;
    int pos=binary_search(ptr,5000*index_field,5000*index_field+4999,index_field_value);
    pos=pos-index_field*5000;
    return linear_interpolator(ptr[5000*index_field+pos],ptr[5000*value_field+pos],ptr[5000*index_field+pos+1],ptr[5000*value_field+pos+1],index_field_value);
}

double linear_interpolator(double x1,double y1,double x2,double y2,double xc)
{
    double m=(y1-y2)/(x1-x2);
    return y1+m*(xc-x1);
}


double dist_3d(double x1,double y1,double z1,double x2,double y2,double z2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}


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
int field_check_eam(int index_field,int value_field)
{
    int ret_val=1;

    if(!(index_field != 0 || index_field!=6)){printf("ERROR:incorrect eam index_feild,returning 0\n");ret_val= 0;}
    if(!((value_field>=1 && value_field<=5) || value_field!=7 || value_field!=8)){printf("ERROR:incorrect eam value feild returning 0\n");ret_val=0;}
    if(index_field == value_field){printf("ERROR:same eam value and eam index field,returning 0\n");ret_val=0;}
    return ret_val;
}

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
                        

int main(void)
{
    FILE *fp_init,*fp1,*fp;
    FILE *f_param;

    double eam_data[45000];
   // int *site[4];
    
    float sites[size];
    int count=0;
    int a[]={12,6,24,12,24,8,48};
    for(int i=0;i<7;i++)
    {
        char fn[10];
        sprintf(fn,"S%dn.mat",i+1);
        printf("Reading from file %s\n",fn);
        fp=fopen(fn,"r");
        if(fp==NULL)printf("unable to open file S%dn.mat\n",i+1);
        for(int j=0;j<a[i]*3 && count<size;j++,count++)
        {
            fscanf(fp,"%f",&sites[count]);
        } 
    }

    fp_init=fopen("file_list.txt","r");

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
    //example usage=..[ double a=eam_data_interpolation_func(eam_data,0,2,1.6); ]
    

    f_param=fopen("parameters.txt","r");
    int no_of_params;
    double *parameter;

    fscanf(f_param,"%d ",&no_of_params);
    parameter=(double*)malloc(no_of_params*sizeof(double));
    int n=read_parameter("parameters.txt",parameter);
    
    double conc_g_p=parameter[5]/100;
    double Nx=parameter[0];
    double Ny=parameter[1];
    double Nz=parameter[2];
    
    int *atoms;
    int no_of_atoms=Nx*Ny*Nz*4;
    atoms=(int *)malloc(4*Nx*Ny*Nz*sizeof(int));
    
    FILE *fp_input;
    fp_input=fopen("input_data.txt","r");
    printf("Reading file input_data.txt\n");
    for(int i=0;i<no_of_atoms;i++)
        fscanf(fp_input,"%d",&atoms[i]);

    //At this point all the needed data has been extracted from the various data files 
    //testing 21:09 06/13/2016 
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng * k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    k = gsl_rng_alloc (T);
    int N_MCS=1000;
    double lattice_parameter=4.00;

    for(int monte_carlo_steps=0;monte_carlo_steps<N_MCS;monte_carlo_steps++)
    {
        double u=gsl_rng_uniform(r);
        int init_pos=u*no_of_atoms;
        int a=init_pos%4;
        int x=(init_pos/4)%(int)Nx;
        int y=(init_pos/(4*(int)Nx))%(int)Ny;
        int z=(init_pos/(4*(int)Nx*(int)Ny))%(int)Nz;


        double xc,yc,zc;
        //xc,yc,zc are the miller indcies coordinate
        if(a==0)
        {
            xc=(double)x; 
            yc=(double)y;
            zc=(double)z;
        }
        else
        {
            xc=(double)x+sites[3*a-3];
            yc=(double)y+sites[3*a-2];
            zc=(double)x+sites[3*a-1];
        }

        int A=atoms[init_pos];          //0 if aluminium, 1 if Nickel
        int site_atom;
        //printf("%d %d %d %d %d %f %f %f\n",init_pos,a,x,y,z,xc,yc,zc);
        double dx,dy,dz,xi,yi,zi,r,dxi,dyi,dzi;
        int as,xs,ys,zs;

        double electron_density_d=0;
        double electron_density_s=0;
        double EAM_energy_default=0;
        double EAM_energy_switched=0;
        double pairwise_energy_d=0;
        double pairwise_energy_s=0;
        double energy_change=0;
        for(int i=0;i<402;i+=3)
        {
            dx=sites[i];
            dy=sites[i+1];
            dz=sites[i+2];
            xi=xc+dx;
            yi=yc+dy;
            zi=zc+dz;

            r=dist_3d(xc,yc,zc,xi,yi,zi)*lattice_parameter;
           
            if(xi<0)xi=xi+Nx;
            if(xi>Nx-0.5)xi=xi-Nx;
            if(yi<0)yi=yi+Ny+0.5;
            if(yi>Ny-0.5)yi=yi-Ny;
            if(zi<0)xi=zi+Nz;
            if(zi>Nz-0.5)zi=zi-Nz;

            dxi=xi-floor(xi);
            dyi=yi-floor(yi);
            dzi=zi-floor(zi);
            if(dxi==0 && dzi==0 && dyi==0) as=0;
            else{
                for(int j=0;j<3;j++)
                {
                    if((dxi-sites[3*j]) ==0 &&(dyi-sites[3*j+1]) ==0 &&(dzi-sites[3*j+2]) ==0)as=j+1;
                }
            }
            xs=(int)(xi-dxi);
            ys=(int)(yi-dyi); 
            zs=(int)(zi-dzi);
            int pos=as+4*xs+4*(int)Nx*ys+4*(int)Nx*(int)Ny*zs;
            site_atom=atoms[pos];
//double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)
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
            if(site_atom==A && A==1)
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,3,r);
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,2,r);
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,5,r);
            }
            else if(site_atom==A && A==0)
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,1,r);
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,2,r);
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,4,r);
            }
            else if(site_atom==0 && A==1)
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,2,r);
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,3,r);
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,5,r);
            }
             else if(site_atom==1 && A==0)
            {
                pairwise_energy_d=pairwise_energy_d+eam_data_interpolation_func(eam_data,0,2,r);
                pairwise_energy_s=pairwise_energy_s+eam_data_interpolation_func(eam_data,0,1,r);
                electron_density_d=electron_density_d+eam_data_interpolation_func(eam_data,0,4,r);
            }  
        }
        //segementation fault occured 
        //possible error in eam_interpolation function 
        //relating to inablity to handle
        //---resolved--fixed the error in eam_interpolation_func function 
        if(site_atom==1){
            EAM_energy_default=pairwise_energy_d+eam_data_interpolation_func(eam_data,6,7,electron_density_d);
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,6,8,electron_density_d);
        }
        else{
            EAM_energy_default=pairwise_energy_d+eam_data_interpolation_func(eam_data,6,8,electron_density_d);
            EAM_energy_switched=pairwise_energy_s+eam_data_interpolation_func(eam_data,6,7,electron_density_d);
        }
        energy_change=EAM_energy_switched-EAM_energy_default;

        double chance=gsl_rng_uniform(k);
        double energy_norm_const=(1.38e-23*TEMP)/1.602e-19;
        double probablity=exp(-energy_change/energy_norm_const);
        if(probablity>1){atoms[init_pos]=!atoms[init_pos];}
        else if(chance<probablity){atoms[init_pos]=!atoms[init_pos];}
        

    }
    //printf("%d\n",!0);
    file_write("Result_lattice.txt",atoms,no_of_atoms);
    return 0;
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
