#include<stdio.h>
#include<stdlib.h>
#include "utils.h"

/// scratched for now -- 19 July ---
void parameterFileRead(double** data, char* fileName, int *size) {
  FILE *fp;
  fp = fopen(fileName, "r");
  int arraySize;
  fscanf(fp,"%d",&arraySize);
  *data = (double*)malloc(sizeof(double) * (arraySize+1));
  for(int i = 0; i < arraySize; i++) {
    fscanf(fp, "%le", data[i]);
  }
  return;
}

int main() {
  double *data;
  int size;
  parameterFileRead(data, "parameters.txt");
  printArrayDouble(data);
}
