//============================================================================
// Name        : LSGFD.cpp
// Author      : Tongfei Guo
// Version     :
// Copyright   : Your copyright notice
//============================================================================

#include <iostream>
#include "Controller.h"
using namespace std;
int main(int argc, char *argv[]) {
	PetscInitialize(&argc,&argv,NULL,NULL);
	Controller control;
	control.Start();
	return 0;
}
