// LinearSLAM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "LinearSFMImp.h"


int _tmain(int argc, _TCHAR* argv[])
{
	 CLinearSFMImp ptr;
	 int nMapCount;

	 char* szSt = NULL, *szPose = NULL, *szFeature = NULL, *szData = NULL;

    //========================================================================================
	//Select Dataset

	//Stereo Dataset
	//New College Dataset
	//szData = "../../../DataForC/NC3500_C";			nMapCount = 3499;
	//szPose = "../../../Pose_NC3500.txt";
	//szFeature = "../../../Feature_NC3500.txt";
	
	//Monocular Dataset
	//Aerial Photogrammetric Village
	szData = "../../../DataForC/RS90_C";			nMapCount = 88;
	szPose = "../../../Pose_RS90.txt";
	szFeature = "../../../Feature_RS90.txt";

	//Aerial Photogrammetric College
	//szData = "../../../DataForC/RS468_C";			nMapCount = 466;
	//szPose = "../../../Pose_RS468.txt";
	//szFeature = "../../../Feature_RS468.txt";




	//Run Linear SFM

	//Stereo
	//ptr.runStereo( szSt, szPose, szFeature, szData, nMapCount );

	//Monocular
	ptr.runMono( szSt, szPose, szFeature, szData, nMapCount );

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	system("pause");
	return 0;
}

