/*
   This file is part of DDP-12D110009.

    DDP-12D110009 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DDP-12D110009 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
   \file parameters.h
   \brief Contains functions that operate on ::parameter .
   \author Piyush Divyankar
   \date 31/08/2016
 */

#include "datatypes.h"

#ifndef _PARAMETERS_H
#define _PARAMETERS_H 1
void print_parameters(parameter * a);
parameter * new_parameters(char * fileName);
parameter * createParameterFileFromInput();
void parameterWriteToFile(parameter *);
parameter * _defaultFCCparameter();
parameter * parameterReadFromFile(char * fileName);
void parameterDefaultFile();
#endif
