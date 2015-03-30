//==============================================================================================================
// Imperial College London, United Kingdom
// University of Technology, Sydney, Australia
// 
// Authors:  Liang Zhao         -- liang.zhao@imperial.ac.uk 
// 		     Shoudong Huang     -- Shoudong.Huang@uts.edu.au
// 		     Gamini Dissanayake -- Gamini.Dissanayake@uts.edu.au
// 
//			 Hamlyn Centre for Robotic Surgery
//           Department of Computing
//           Faculty of Engineering
//           Imperial College London, United Kingdom
//
// 		     Centre for Autonomous Systems
// 		     Faculty of Engineering and Information Technology 
// 		     University of Technology, Sydney
// 		     NSW 2007, Australia
// 
// 		  License
// 
// 		  Linear SFM by Liang Zhao, Shoudong Huang, Gamini Dissanayake is licensed under a 
// 		  Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// 
// 		  Please contact Liang Zhao {liang.zhao@imperial.ac.uk} if you have any questions/comments about the code.
//==============================================================================================================

#pragma once
#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>
using namespace Eigen;

#include "suitesparse/cholmod.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>  
using namespace std;

#define PI 3.1415926

typedef enum LinearSFM_DataType
{
	LinearSFM_Mono = 0,
	LinearSFM_Stereo = 1,
}LinearSFM_DataType;

//use sba_crsm structure from SBA (http://www.ics.forth.gr/~lourakis/sba/) to store sparse matrix
struct sba_crsm	
{
	int nr, nc;   
	int nnz;      
	int *val;     
	int *colidx;  
	int *rowptr;  
};

typedef struct LocalMapInfoStereo
{
public:
	int r, Ref;
	int	*stno;

	LocalMapInfoStereo(){};
	LocalMapInfoStereo( int row ){  r = row;	}

	void setDimension( int row ) { r = row; }

	LocalMapInfoStereo& operator = (const LocalMapInfoStereo& lm) 	
	{
		stno = lm.stno;
		r = lm.r;	Ref = lm.Ref;	

		stVal = lm.stVal;
		U = lm.U;
		V = lm.V;
		W = lm.W;
		nW = lm.nW;
		nU = lm.nU;
		Ui = lm.Ui;
		Uj = lm.Uj;
		photo = lm.photo;
		feature = lm.feature;
		FBlock = lm.FBlock;
		m = lm.m;
		n = lm.n;

		FRef = lm.FRef;

		return *this;
	}

	double* stVal;
	double* U, *W, *V;
	int nU, nW;
	int* Ui, *Uj;
	int* photo, *feature;
	int* FBlock;
		
	int m, n;

	int FRef;

} LocalMapInfoStereo;


typedef struct LocalMapInfo
{
public:
	int r, Ref;
	int	*stno;

	LocalMapInfo(){};
	LocalMapInfo( int row ){  r = row;	}

	void setDimension( int row ) { r = row; }

	LocalMapInfo& operator = (const LocalMapInfo& lm) 	
	{
		stno = lm.stno;
		r = lm.r;	Ref = lm.Ref;	

		stVal = lm.stVal;
		U = lm.U;
		V = lm.V;
		W = lm.W;
		nW = lm.nW;
		nU = lm.nU;
		Ui = lm.Ui;
		Uj = lm.Uj;
		photo = lm.photo;
		feature = lm.feature;
		FBlock = lm.FBlock;
		m = lm.m;
		n = lm.n;

		FRef = lm.FRef;
		FScaP = lm.FScaP;
		FFix = lm.FFix;

		ScaP = lm.ScaP;
		Sign = lm.Sign;
		Fix = lm.Fix;

		return *this;
	}

	double* stVal;
	double* U, *W, *V;
	int nU, nW;
	int* Ui, *Uj;
	int* photo, *feature;
	int* FBlock;
		
	int m, n;

	int FRef, FScaP, FFix;

	int ScaP, Sign, Fix;

} LocalMapInfo;


class CLinearSFMImp
{
public:
	CLinearSFMImp(void);
	~CLinearSFMImp(void);

	//========================================================================================
	void	runStereo( char* szState, char* szPose, char* szFeature, char* szData, int nMapCount );
	void	runMono( char* szState, char* szPose, char* szFeature, char* szData, int nMapCount );

	bool	lmj_parseArgs( int argc, char* argv[] );
	bool    run( int argc, char** argv );
	void	lmj_printHelp();
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

private:
		
	//========================================================================================	
	void	lmj_SaveStateVector( char* szSt, double* st, int* stno, int n );
	void	lmj_SavePoses_3DPF( char* szPose, char* szFea, int* stno, double* st, int n );

public:
	//Stereo
	void	lmj_loadLocalMapsStereo( int nMapCount );
	void	lmj_PF3D_Divide_ConquerStereo( int nLocalMapCount );
	void	lmj_Transform_PF3DStereo( LocalMapInfoStereo& GMap_End, int Ref );
	void	lmj_LinearLS_PF3DStereo( LocalMapInfoStereo& GMap_End, LocalMapInfoStereo& GMap_Cur );
	void	lmj_readInformationStereo( LocalMapInfoStereo& GMap_End, char* szPath );
	void    lmj_solveLinearSFMStereo( double* stVal, double* eb, double* ea,  double* U, double*W, double* V, int* Ui, int* Uj, int* photo, int* feature, int m, int n, int nU, int nW );
	void	pba_constructAuxCSSLM( int *Ap, int *Aii, int m, char* smask );
	void	pba_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, sba_crsm& Sidxij, bool init, int m, char* smask);
	bool	pba_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering, cholmod_common m_cS, int m_nS, int m );
	void	pba_inverseV( double* V, int m, int n );
	void	pba_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb, int m, int n, int* mapCor, int* photo );

	
	//Monocular
	void	lmj_Transform_PF3DMono( LocalMapInfo& GMap_End, int Ref, int ScaP, int Fix );
	void	lmj_loadLocalMapsMono( int nMapCount );
	void	lmj_PF3D_Divide_ConquerMono( int nLocalMapCount );
	void	lmj_LinearLS_PF3DMono( LocalMapInfo& GMap_End, LocalMapInfo& GMap_Cur );
	void	lmj_readInformationMono( LocalMapInfo& GMap_End, char* szPath );
	void    lmj_solveLinearSFMMono( double* stVal, double* eb, double* ea,  double* U, double*W, double* V, int* Ui, int* Uj, int* photo, int* feature, int m, int n, int nU, int nW, int Ref, int ScaP, int Fix, int Sign, int FixBlk );
	void	pba_constructAuxCSSGN( int *Ap, int *Aii, int m, char* smask, int Ref );
	void	pba_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, sba_crsm& Sidxij, bool init, int m, char* smask, int Ref, int ScaP, int Fix);
	bool	pba_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering, cholmod_common m_cS, int m_nS, int m, int FixBlk );
	

private:
	char*			m_szPath;
	char*			m_szInfo;
	char*			m_szSt;
	char*			m_szPose;
	char*			m_szFeature;
	LocalMapInfoStereo*	m_LMsetS;
	LocalMapInfoStereo	m_GMapS;
	LocalMapInfo*	m_LMset;
	LocalMapInfo	m_GMap;
	int				m_nMapCount;
	LinearSFM_DataType m_data;

private:
	//use cholmod mathmatic software package to solve linear equation
	cholmod_sparse *m_sparseS; 
	cholmod_factor *m_factorS; 
	cholmod_common m_cS; 
	cholmod_dense  *m_sparseR, *m_sparseE;


private:
	char *m_umask, *m_wmask;

};

