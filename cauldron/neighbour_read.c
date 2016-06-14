#include<stdio.h>
#include<stdlib.h>
#include<string.h>


int main(void)
{
    FILE *fp;
    int size=(12+6+24+12+24+8+48)*3;
    int a[]={12,6,24,12,24,8,48};
    double *neighbours;
    neighbours=(double*)malloc(size*sizeof(double));
    int count_g=0;
    for(int i=1;i<=7;i++)
    {
        char *fname;
        int count_l=0;
        sprintf(fname,"S%dn.mat",i);
        fp=fopen(fname,"r");
        if(fp==NULL)printf("Unable to open file S%dn.mat\n",i);
        printf("1 file taken in\n");
        while(count_l<3*a[i-1])
        {
            fscanf(fp,"%le",&neighbours[count_g]);
            count_g++;
            count_l++;
        }
        fclose(fp);
        printf("1 file taken in\n");
    }
    count_g=0;
    while(count_g<size)
    {
        printf("%le\t%le%le\n",neighbours[count_g],neighbours[count_g+1],neighbours[count_g+2]);
        count_g+=3;
    }
}
