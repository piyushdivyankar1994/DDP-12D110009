#include <stdio.h>

void file_write(char *filename,int *data,int data_points)
{
    FILE *fp;
    fp=fopen(filename,"w");
    for(int i=0;i<data_points;i+=4)
    {
        fprintf(fp,"%d %d %d %d \n",data[i],data[i+1],data[i+2],data[i+3]);
    }
    printf("%s",filename);
}
