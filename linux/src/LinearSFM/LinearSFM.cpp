// LinearSLAM.cpp : Defines the entry point for the console application.
//

#include <iostream>
using namespace std;
#include "LinearSFMImp.h"


int main(int argc, char* argv[])
{
	 CLinearSFMImp ptr;
	 
	 //ptr.runStereo(argc, argv);

	 ptr.run(argc, argv);

	return 0;
}
