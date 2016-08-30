#include "datatypes.h"

#ifndef _PARAMETERS_H
#define _PARAMETERS_H 1
  void print_parameters(parameter* a);
  parameter* new_parameters(char* filename);
  void createParameterFileFromInput();
  parameter* _defaultFCCparameter();
#endif
