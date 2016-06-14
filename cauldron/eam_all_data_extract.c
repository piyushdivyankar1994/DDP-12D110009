#include<stdio.h>
#include<stdlib.h>

#define size 402
#define N 10 

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


int read_parameter(char*,double*);

double eam_data_interpolation_func(double *ptr, int index_field, int value_field, double index_field_value)
{
    if(field_check_eam(index_field, value_field)==0)return 0;
    int pos=binary_search(ptr,5000*index_field,5000*index_field+4999,index_field_value);
    printf("\n%le\t%le\t%le\t%le\t%le\n",ptr[5000*index_field+pos],ptr[5000*value_field+pos],ptr[5000*index_field+pos+1],ptr[5000*value_field+pos+1],index_field_value);
    return linear_interpolator(ptr[5000*index_field+pos],ptr[5000*value_field+pos],ptr[5000*index_field+pos+1],ptr[5000*value_field+pos+1],index_field_value);
}

double linear_interpolator(double x1,double y1,double x2,double y2,double xc)
{
    double m=(y1-y2)/(x1-x2);
    return y1+m*(xc-x1);
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
 
    return 0;
}
