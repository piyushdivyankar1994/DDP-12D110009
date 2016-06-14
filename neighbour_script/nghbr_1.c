#include<stdio.h>
#include<stdlib.h>

#define  size 402
int main(void)
{
    FILE *fp;
    fp=fopen("S1n.mat","r");
    float sites[size];
    int count=0;
    int a[]={12,6,24,12,24,8,48};
    for(int i=0;i<7;i++)
    {
        char fn[10];
        sprintf(fn,"S%dn.mat",i+1);
        printf("%s\n",fn);
        fp=fopen(fn,"r");
        if(fp==NULL)printf("unable to open file S%dn.mat\n",i+1);
        for(int j=0;j<a[i]*3 && count<size;j++,count++)
        {
            fscanf(fp,"%f",&sites[count]);
            printf("%f\n",sites[count]);
        }
    }
    return 0;
}
