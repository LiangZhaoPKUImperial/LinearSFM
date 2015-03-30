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

#include "stdafx.h"
#include "LinearSFMImp.h"

//cite from SBA code in order to search sub matrix of S according to (i,j)
static void sba_crsm_alloc(struct sba_crsm *sm, int nr, int nc, int nnz)
{
	int msz;
	sm->nr=nr;
	sm->nc=nc;
	sm->nnz=nnz;
	msz=2*nnz+nr+1;
	sm->val=(int *)malloc(msz*sizeof(int));  /* required memory is allocated in a single step */
	if(!sm->val){
		fprintf(stderr, "memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
		exit(1);
	}
	sm->colidx=sm->val+nnz;
	sm->rowptr=sm->colidx+nnz;
}

static void sba_crsm_free(struct sba_crsm *sm)
{
	sm->nr=sm->nc=sm->nnz=-1;
	free(sm->val);
	sm->val=sm->colidx=sm->rowptr=NULL;
}


/* returns the index of the (i, j) element. No bounds checking! */
static int sba_crsm_elmidx(struct sba_crsm *sm, int i, int j)
{
	register int low, high, mid, diff;

	low=sm->rowptr[i];
	high=sm->rowptr[i+1]-1;

	/* binary search for finding the element at column j */
	while(low<=high)
	{
		mid=(low+high)>>1; //(low+high)/2;
		diff=j-sm->colidx[mid];
		if(diff<0)
			high=mid-1;
		else if(diff>0)
			low=mid+1;
		else
			return mid;
	}

	return -1; /* not found */
}



CLinearSFMImp::CLinearSFMImp(void)
{
	cholmod_start (&m_cS);			
}


CLinearSFMImp::~CLinearSFMImp(void)
{
	cholmod_finish (&m_cS);	
}

bool IsOdd (int i) {  return (i<0);}

void CLinearSFMImp::runStereo( char* szState, char* szPose, char* szFeature, char* szData, int nMapCount )
{
	m_szPath = szData;
	m_szSt   = szState; 
	m_szPose = szPose;
	m_szFeature = szFeature;

	m_LMsetS = new LocalMapInfoStereo[nMapCount];
	lmj_loadLocalMapsStereo( nMapCount );
	
	lmj_PF3D_Divide_ConquerStereo( nMapCount );
	
	free( m_LMsetS );

	system("pause");
}

void CLinearSFMImp::lmj_loadLocalMapsStereo( int nMapCount )
{
	FILE* fpSt;
	int i, rows, tmp, n;

	for ( i = 0; i < nMapCount; i++)
	{
		char szSt[200];
		strcpy( szSt, m_szPath );	
		char szPathExt[20];

		sprintf( szPathExt, "/localmap_%d.txt", i+1 );
		strcat( szSt, szPathExt );
		
		lmj_readInformationStereo( m_LMsetS[i], szSt );			
	}
}

static void lmj_RMatrixYPR22( double* R, double Alpha, double Beta, double Gamma )
{
	R[0] = cos(Beta)*cos(Alpha);
	R[1] = cos(Beta)*sin(Alpha);
	R[2] = -sin(Beta);
	R[3] = sin(Gamma)*sin(Beta)*cos(Alpha)-cos(Gamma)*sin(Alpha);
	R[4] = sin(Gamma)*sin(Beta)*sin(Alpha)+cos(Gamma)*cos(Alpha);
	R[5] = sin(Gamma)*cos(Beta);
	R[6] = cos(Gamma)*sin(Beta)*cos(Alpha)+sin(Gamma)*sin(Alpha);
	R[7] = cos(Gamma)*sin(Beta)*sin(Alpha)-sin(Gamma)*cos(Alpha);
	R[8] = cos(Gamma)*cos(Beta);
}

static void lmj_InvRotMatrixYPR22T(double*R, double &alpha, double &beta, double &gamma )
{
	beta = atan2( -R[6], sqrt(R[0]*R[0]+R[3]*R[3]) );

	if ( cos(beta)==0 )
	{
		alpha = 0; 
		beta = PI/2;
		gamma = atan2( R[3], R[4]);
	}
	else
	{
		alpha = atan2( R[3]/cos(beta), R[0]/cos(beta) );
		gamma = atan2( R[7]/cos(beta), R[8]/cos(beta) );
	}
}

static void lmj_InvRotMatrixYPR22(double*R, double &alpha, double &beta, double &gamma )
{
	beta = atan2( -R[2], sqrt(R[0]*R[0]+R[1]*R[1]) );

	if ( cos(beta)==0 )
	{
		alpha = 0; 
		beta = PI/2;
		gamma = atan2( R[1], R[4]);
	}
	else
	{
		alpha = atan2( R[1]/cos(beta), R[0]/cos(beta) );
		gamma = atan2( R[5]/cos(beta), R[8]/cos(beta) );
	}
}

static void lmj_Rderivation( double Alpha, double Beta, double Gamma, double* matR, double* dRA, double* dRB, double*dRG )
{
	matR[0] = cos(Beta) * cos(Alpha);
	matR[1] = cos(Beta) * sin(Alpha);
	matR[2] = -sin(Beta);
	matR[3] = sin(Gamma)*sin(Beta)*cos(Alpha)-cos(Gamma)*sin(Alpha);
	matR[4] = sin(Gamma)*sin(Beta)*sin(Alpha)+cos(Gamma)*cos(Alpha);
	matR[5] = sin(Gamma)*cos(Beta);
	matR[6] = cos(Gamma)*sin(Beta)*cos(Alpha) + sin(Gamma)*sin(Alpha);
	matR[7] = cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha);
	matR[8] = cos(Gamma)*cos(Beta);

	double matRG[9], matRB[9], matRA[9];
	matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
	matRG[3] = 0;		matRG[4] = cos(Gamma);		matRG[5] = sin(Gamma);
	matRG[6] = 0;		matRG[7] = -sin(Gamma);		matRG[8] = cos(Gamma);

	matRB[0] = cos(Beta);		matRB[1] = 0;		matRB[2] = -sin(Beta);
	matRB[3] = 0;				matRB[4] = 1;		matRB[5] = 0;
	matRB[6] = sin(Beta);		matRB[7] = 0;		matRB[8] = cos(Beta);

	matRA[0] = cos(Alpha);		matRA[1] = sin(Alpha);			matRA[2] = 0;
	matRA[3] = -sin(Alpha);		matRA[4] = cos(Alpha);			matRA[5] = 0;
	matRA[6] = 0;				matRA[7] = 0;					matRA[8] = 1;

	double matDRG[9], matDRB[9], matDRA[9];
	matDRG[0] = 0;		matDRG[1] = 0;				matDRG[2] = 0;
	matDRG[3] = 0;		matDRG[4] = -sin(Gamma);	matDRG[5] = cos(Gamma);
	matDRG[6] = 0;		matDRG[7] = -cos(Gamma);	matDRG[8] = -sin(Gamma);

	matDRB[0] = -sin(Beta);		matDRB[1] = 0;		matDRB[2] = -cos(Beta);
	matDRB[3] = 0;				matDRB[4] = 0;		matDRB[5] = 0;
	matDRB[6] = cos(Beta);		matDRB[7] = 0;		matDRB[8] = -sin(Beta);

	matDRA[0] = -sin(Alpha);		matDRA[1] = cos(Alpha);			matDRA[2] = 0;
	matDRA[3] = -cos(Alpha);		matDRA[4] = -sin(Alpha);		matDRA[5] = 0;
	matDRA[6] = 0;					matDRA[7] = 0;					matDRA[8] = 0;

	//dRG
	double tmp1[9];
	tmp1[0] = matDRG[0]*matRB[0]+matDRG[1]*matRB[3]+matDRG[2]*matRB[6];
	tmp1[1] = matDRG[0]*matRB[1]+matDRG[1]*matRB[4]+matDRG[2]*matRB[7];
	tmp1[2] = matDRG[0]*matRB[2]+matDRG[1]*matRB[5]+matDRG[2]*matRB[8];
	tmp1[3] = matDRG[3]*matRB[0]+matDRG[4]*matRB[3]+matDRG[5]*matRB[6];
	tmp1[4] = matDRG[3]*matRB[1]+matDRG[4]*matRB[4]+matDRG[5]*matRB[7];
	tmp1[5] = matDRG[3]*matRB[2]+matDRG[4]*matRB[5]+matDRG[5]*matRB[8];
	tmp1[6] = matDRG[6]*matRB[0]+matDRG[7]*matRB[3]+matDRG[8]*matRB[6];
	tmp1[7] = matDRG[6]*matRB[1]+matDRG[7]*matRB[4]+matDRG[8]*matRB[7];
	tmp1[8] = matDRG[6]*matRB[2]+matDRG[7]*matRB[5]+matDRG[8]*matRB[8];

	dRG[0] = tmp1[0]*matRA[0]+tmp1[1]*matRA[3]+tmp1[2]*matRA[6];
	dRG[1] = tmp1[0]*matRA[1]+tmp1[1]*matRA[4]+tmp1[2]*matRA[7];
	dRG[2] = tmp1[0]*matRA[2]+tmp1[1]*matRA[5]+tmp1[2]*matRA[8];
	dRG[3] = tmp1[3]*matRA[0]+tmp1[4]*matRA[3]+tmp1[5]*matRA[6];
	dRG[4] = tmp1[3]*matRA[1]+tmp1[4]*matRA[4]+tmp1[5]*matRA[7];
	dRG[5] = tmp1[3]*matRA[2]+tmp1[4]*matRA[5]+tmp1[5]*matRA[8];
	dRG[6] = tmp1[6]*matRA[0]+tmp1[7]*matRA[3]+tmp1[8]*matRA[6];
	dRG[7] = tmp1[6]*matRA[1]+tmp1[7]*matRA[4]+tmp1[8]*matRA[7];
	dRG[8] = tmp1[6]*matRA[2]+tmp1[7]*matRA[5]+tmp1[8]*matRA[8];

	//dRB
	tmp1[0] = matRG[0]*matDRB[0]+matRG[1]*matDRB[3]+matRG[2]*matDRB[6];
	tmp1[1] = matRG[0]*matDRB[1]+matRG[1]*matDRB[4]+matRG[2]*matDRB[7];
	tmp1[2] = matRG[0]*matDRB[2]+matRG[1]*matDRB[5]+matRG[2]*matDRB[8];
	tmp1[3] = matRG[3]*matDRB[0]+matRG[4]*matDRB[3]+matRG[5]*matDRB[6];
	tmp1[4] = matRG[3]*matDRB[1]+matRG[4]*matDRB[4]+matRG[5]*matDRB[7];
	tmp1[5] = matRG[3]*matDRB[2]+matRG[4]*matDRB[5]+matRG[5]*matDRB[8];
	tmp1[6] = matRG[6]*matDRB[0]+matRG[7]*matDRB[3]+matRG[8]*matDRB[6];
	tmp1[7] = matRG[6]*matDRB[1]+matRG[7]*matDRB[4]+matRG[8]*matDRB[7];
	tmp1[8] = matRG[6]*matDRB[2]+matRG[7]*matDRB[5]+matRG[8]*matDRB[8];

	dRB[0] = tmp1[0]*matRA[0]+tmp1[1]*matRA[3]+tmp1[2]*matRA[6];
	dRB[1] = tmp1[0]*matRA[1]+tmp1[1]*matRA[4]+tmp1[2]*matRA[7];
	dRB[2] = tmp1[0]*matRA[2]+tmp1[1]*matRA[5]+tmp1[2]*matRA[8];
	dRB[3] = tmp1[3]*matRA[0]+tmp1[4]*matRA[3]+tmp1[5]*matRA[6];
	dRB[4] = tmp1[3]*matRA[1]+tmp1[4]*matRA[4]+tmp1[5]*matRA[7];
	dRB[5] = tmp1[3]*matRA[2]+tmp1[4]*matRA[5]+tmp1[5]*matRA[8];
	dRB[6] = tmp1[6]*matRA[0]+tmp1[7]*matRA[3]+tmp1[8]*matRA[6];
	dRB[7] = tmp1[6]*matRA[1]+tmp1[7]*matRA[4]+tmp1[8]*matRA[7];
	dRB[8] = tmp1[6]*matRA[2]+tmp1[7]*matRA[5]+tmp1[8]*matRA[8];

	//dRA
	tmp1[0] = matRG[0]*matRB[0]+matRG[1]*matRB[3]+matRG[2]*matRB[6];
	tmp1[1] = matRG[0]*matRB[1]+matRG[1]*matRB[4]+matRG[2]*matRB[7];
	tmp1[2] = matRG[0]*matRB[2]+matRG[1]*matRB[5]+matRG[2]*matRB[8];
	tmp1[3] = matRG[3]*matRB[0]+matRG[4]*matRB[3]+matRG[5]*matRB[6];
	tmp1[4] = matRG[3]*matRB[1]+matRG[4]*matRB[4]+matRG[5]*matRB[7];
	tmp1[5] = matRG[3]*matRB[2]+matRG[4]*matRB[5]+matRG[5]*matRB[8];
	tmp1[6] = matRG[6]*matRB[0]+matRG[7]*matRB[3]+matRG[8]*matRB[6];
	tmp1[7] = matRG[6]*matRB[1]+matRG[7]*matRB[4]+matRG[8]*matRB[7];
	tmp1[8] = matRG[6]*matRB[2]+matRG[7]*matRB[5]+matRG[8]*matRB[8];

	dRA[0] = tmp1[0]*matDRA[0]+tmp1[1]*matDRA[3]+tmp1[2]*matDRA[6];
	dRA[1] = tmp1[0]*matDRA[1]+tmp1[1]*matDRA[4]+tmp1[2]*matDRA[7];
	dRA[2] = tmp1[0]*matDRA[2]+tmp1[1]*matDRA[5]+tmp1[2]*matDRA[8];
	dRA[3] = tmp1[3]*matDRA[0]+tmp1[4]*matDRA[3]+tmp1[5]*matDRA[6];
	dRA[4] = tmp1[3]*matDRA[1]+tmp1[4]*matDRA[4]+tmp1[5]*matDRA[7];
	dRA[5] = tmp1[3]*matDRA[2]+tmp1[4]*matDRA[5]+tmp1[5]*matDRA[8];
	dRA[6] = tmp1[6]*matDRA[0]+tmp1[7]*matDRA[3]+tmp1[8]*matDRA[6];
	dRA[7] = tmp1[6]*matDRA[1]+tmp1[7]*matDRA[4]+tmp1[8]*matDRA[7];
	dRA[8] = tmp1[6]*matDRA[2]+tmp1[7]*matDRA[5]+tmp1[8]*matDRA[8];
}

static void lmj_dRi( double* dRid, double* dRi, double* Ri )
{
	double F1, F2, F3, F4, F5, dAdF1, dBdF2, dGdF3, dF1d, dF2d, dF3d, dF4d, dF5d, dF4dF5;

	F1 = Ri[1]/Ri[0];
	F3 = Ri[5]/Ri[8];
	F5 = Ri[0]*Ri[0] + Ri[1]*Ri[1];
	F4 = sqrt(F5);
	F2 = -Ri[2]/F4;

	dAdF1 = 1.0/(1+F1*F1);
	dBdF2 = 1.0/(1+F2*F2);
	dGdF3 = 1.0/(1+F3*F3);

	dF1d = (dRi[1]*Ri[0]-Ri[1]*dRi[0])/(Ri[0]*Ri[0]);
	dF3d = (dRi[5]*Ri[8]-Ri[5]*dRi[8])/(Ri[8]*Ri[8]);

	dF4dF5 = 1.0/(2*sqrt(F5));
	dF5d = 2*Ri[0]*dRi[0] + 2*Ri[1]*dRi[1];
	dF4d = dF4dF5 *dF5d;

	dF2d = (-dRi[2]*F4+Ri[2]*dF4d)/F5;
	dRid[0] = dAdF1*dF1d;
	dRid[1] = dBdF2*dF2d;
	dRid[2] = dGdF3*dF3d;
}

static void lmj_dRiTT( double* dRid, double* dRi, double* Ri )
{
	double F1, F2, F3, F4, F5, dAdF1, dBdF2, dGdF3, dF1d, dF2d, dF3d, dF4d, dF5d, dF4dF5;

	F1 = Ri[3]/Ri[0];
	F3 = Ri[7]/Ri[8];
	F5 = Ri[0]*Ri[0] + Ri[3]*Ri[3];
	F4 = sqrt(F5);
	F2 = -Ri[6]/F4;

	dAdF1 = 1.0/(1+F1*F1);
	dBdF2 = 1.0/(1+F2*F2);
	dGdF3 = 1.0/(1+F3*F3);

	dF1d = (dRi[3]*Ri[0]-Ri[3]*dRi[0])/(Ri[0]*Ri[0]);
	dF3d = (dRi[7]*Ri[8]-Ri[7]*dRi[8])/(Ri[8]*Ri[8]);

	dF4dF5 = 1.0/(2*sqrt(F5));
	dF5d = 2*Ri[0]*dRi[0] + 2*Ri[3]*dRi[3];
	dF4d = dF4dF5 *dF5d;

	dF2d = (-dRi[6]*F4+Ri[6]*dF4d)/F5;
	dRid[0] = dAdF1*dF1d;
	dRid[1] = dBdF2*dF2d;
	dRid[2] = dGdF3*dF3d;
}

static void lmj_TimesRRT( double*R3, double* R1, double* R2)
{
	R3[0] = R1[0]*R2[0] + R1[1]*R2[1] + R1[2]*R2[2];
	R3[1] = R1[0]*R2[3] + R1[1]*R2[4] + R1[2]*R2[5];
	R3[2] = R1[0]*R2[6] + R1[1]*R2[7] + R1[2]*R2[8];
	R3[3] = R1[3]*R2[0] + R1[4]*R2[1] + R1[5]*R2[2];
	R3[4] = R1[3]*R2[3] + R1[4]*R2[4] + R1[5]*R2[5];
	R3[5] = R1[3]*R2[6] + R1[4]*R2[7] + R1[5]*R2[8];
	R3[6] = R1[6]*R2[0] + R1[7]*R2[1] + R1[8]*R2[2];
	R3[7] = R1[6]*R2[3] + R1[7]*R2[4] + R1[8]*R2[5];
	R3[8] = R1[6]*R2[6] + R1[7]*R2[7] + R1[8]*R2[8];
}

void CLinearSFMImp::lmj_Transform_PF3DStereo( LocalMapInfoStereo& GMap_End, int Ref )
{

	if (m_GMapS.Ref == Ref )
    {
         GMap_End = m_GMapS;
     }
    else
    {	
		int pos, i, N; 
		double t[3], Alpha, Beta, Gamma, R[9], R2[9], R3[9], dRA[9], dRB[9], dRG[9];
		double dA[3], dB[3], dG[3], dRA2[9], dRB2[9], dRG2[9];
		double t2[3],  Alpha2, Beta2, Gamma2, Ri[9]; //, tmp1[3], tmp2[3], tmp3[3];
		double tmp1[6], tmp2[6], tmp3[6], tmp4[6], tmp5[6], tmp6[6], ttmp1[6], ttmp2[6], ttmp3[6], ttmp4[6], ttmp5[6], ttmp6[6];
		double ddA2[3], ddB2[3], ddG2[3], ddA[3], ddB[3], ddG[3];
		double* ptr1, *ptr2;
		int* iter;
		int* stno;
		N = 6*m_GMapS.m + 3*m_GMapS.n;
		GMap_End.r = m_GMapS.r;
        
		int m = m_GMapS.m;
		int n  = m_GMapS.n;
		double* m_U, *m_W, *m_V;
		int* m_Ui, *m_Uj, *m_photo, *m_feature;
		int m_nU, m_nW, m_nFea;	

		m_U = m_GMapS.U;
		m_V = m_GMapS.V;
		m_W = m_GMapS.W;
		m_Ui= m_GMapS.Ui;
		m_Uj= m_GMapS.Uj;
		m_photo = m_GMapS.photo;
		m_feature = m_GMapS.feature;
		m_nU = m_GMapS.nU;
		m_nW = m_GMapS.nW;
		m_nFea = m_GMapS.n;
        ptr1 = m_GMapS.stVal;
        stno = m_GMapS.stno;

        iter = find( stno, stno+N, -Ref );
        pos = iter - stno;

        t[0] = ptr1[pos];
        t[1] = ptr1[pos+1];
        t[2] = ptr1[pos+2];

        Alpha  = ptr1[pos+3];
        Beta   = ptr1[pos+4];
        Gamma  = ptr1[pos+5];

        lmj_RMatrixYPR22( R, Alpha, Beta, Gamma );

		GMap_End.m = m_GMapS.m;
		GMap_End.n = m_GMapS.n;
		GMap_End.stno = (int*)malloc( (6*m+3*n)*sizeof(int) );
		GMap_End.stVal = (double*)malloc( (6*m+3*n)*sizeof(double) );

		GMap_End.FBlock = (int*)malloc(n*sizeof(int));
		memset( GMap_End.FBlock, -1, n * sizeof(int) ); 
		GMap_End.Ref = Ref;

		GMap_End.FRef = m_GMapS.FRef;

		ptr2 = GMap_End.stVal;
        memcpy( GMap_End.stno, stno, (6*m+3*n)*sizeof(int) );

        GMap_End.stno[pos] = GMap_End.stno[pos+1] = GMap_End.stno[pos+2] =
        GMap_End.stno[pos+3] = GMap_End.stno[pos+4] = GMap_End.stno[pos+5] = -m_GMapS.Ref;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================
        for ( i=0; i<N; i++ )
        {
			if ( stno[i] <=0 )
            {
				if ( i == pos )
                {
					ptr2[i]   = -( R[0]*t[0]+R[1]*t[1]+R[2]*t[2] );
                    ptr2[i+1] = -( R[3]*t[0]+R[4]*t[1]+R[5]*t[2] );
                    ptr2[i+2] = -( R[6]*t[0]+R[7]*t[1]+R[8]*t[2] );

                    lmj_InvRotMatrixYPR22T( R, ptr2[i+3], ptr2[i+4], ptr2[i+5] );                                   
				}
                else
				{
					ptr2[i]   = ( R[0]*(ptr1[i]-t[0])+R[1]*(ptr1[i+1]-t[1])+R[2]*(ptr1[i+2]-t[2]) );
					ptr2[i+1] = ( R[3]*(ptr1[i]-t[0])+R[4]*(ptr1[i+1]-t[1])+R[5]*(ptr1[i+2]-t[2]) );
					ptr2[i+2] = ( R[6]*(ptr1[i]-t[0])+R[7]*(ptr1[i+1]-t[1])+R[8]*(ptr1[i+2]-t[2]) );

					lmj_RMatrixYPR22( R2, ptr1[i+3], ptr1[i+4], ptr1[i+5] );

					lmj_TimesRRT( R3, R2, R );
					lmj_InvRotMatrixYPR22( R3, ptr2[i+3], ptr2[i+4], ptr2[i+5] );
				}

					i+=5;
			}
            else
			{
				ptr2[i]   = ( R[0]*(ptr1[i]-t[0])+R[1]*(ptr1[i+1]-t[1])+R[2]*(ptr1[i+2]-t[2]) );
				ptr2[i+1] = ( R[3]*(ptr1[i]-t[0])+R[4]*(ptr1[i+1]-t[1])+R[5]*(ptr1[i+2]-t[2]) );
				ptr2[i+2] = ( R[6]*(ptr1[i]-t[0])+R[7]*(ptr1[i+1]-t[1])+R[8]*(ptr1[i+2]-t[2]) );

				i = i+2;
			}
		}
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================
        t[0] =  ptr2[pos];
        t[1] =  ptr2[pos+1];
        t[2]  = ptr2[pos+2];

        Alpha = ptr2[pos+3];
        Beta  = ptr2[pos+4];
        Gamma = ptr2[pos+5];

        lmj_Rderivation( Alpha, Beta, Gamma, R, dRA, dRB, dRG );

        lmj_dRiTT( dA, dRA, R );
        lmj_dRiTT( dB, dRB, R );
        lmj_dRiTT( dG, dRG, R );
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================
        int size1 = m*6*6 + n*3*3;
        double *J1 = (double*)malloc( size1*sizeof(double) );
        memset( J1, 0, size1*sizeof(double) );
        int size2 = m*6*6 + n*3*6;
        double *J2 = (double*)malloc( size2*sizeof(double) );
        memset( J2, 0, size2*sizeof(double) );

        double* pJ1, *pJ2;
        int Id;

        for ( i = 0; i < N; i++ )
		{
			if (GMap_End.stno[i]<=0)
			{
				if ( i == pos )
				{
					tmp1[0] = (dRA[0]*t[0]+dRA[1]*t[1]+dRA[2]*t[2]);
					tmp1[1] = (dRA[3]*t[0]+dRA[4]*t[1]+dRA[5]*t[2]);
					tmp1[2] = (dRA[6]*t[0]+dRA[7]*t[1]+dRA[8]*t[2]);

					tmp2[0] = (dRB[0]*t[0]+dRB[1]*t[1]+dRB[2]*t[2]);
					tmp2[1] = (dRB[3]*t[0]+dRB[4]*t[1]+dRB[5]*t[2]);
					tmp2[2] = (dRB[6]*t[0]+dRB[7]*t[1]+dRB[8]*t[2]);

					tmp3[0] = (dRG[0]*t[0]+dRG[1]*t[1]+dRG[2]*t[2]);
					tmp3[1] = (dRG[3]*t[0]+dRG[4]*t[1]+dRG[5]*t[2]);
					tmp3[2] = (dRG[6]*t[0]+dRG[7]*t[1]+dRG[8]*t[2]);

					Id = i/6;

					pJ1 = J1 + Id * 6*6; 

					pJ1[0]  += -R[0];                               //triI[sum] = i;                        triJ[sum] = pos;        triVal[sum] = -R[0];    
					pJ1[6]  += -R[3];                               //triI[sum] = i+1;              triJ[sum] = pos;        triVal[sum] = -R[3];    
					pJ1[12] += -R[6];                               //triI[sum] = i+2;              triJ[sum] = pos;        triVal[sum] = -R[6];    
					pJ1[1]  += -R[1];                               //triI[sum] = i;                        triJ[sum] = pos+1;      triVal[sum] = -R[1];    
                    pJ1[7]  += -R[4];                               //triI[sum] = i+1;              triJ[sum] = pos+1;      triVal[sum] = -R[4];    
                    pJ1[13] += -R[7];                               //triI[sum] = i+2;              triJ[sum] = pos+1;      triVal[sum] = -R[7];
                    pJ1[2]  += -R[2];                               //triI[sum] = i;                        triJ[sum] = pos+2;      triVal[sum] = -R[2];            
                    pJ1[8]  += -R[5];                               //triI[sum] = i+1;              triJ[sum] = pos+2;      triVal[sum] = -R[5];
                    pJ1[14] += -R[8];                               //triI[sum] = i+2;              triJ[sum] = pos+2;      triVal[sum] = -R[8];    
                    pJ1[3]  += -tmp1[0];                    //triI[sum] = i;                        triJ[sum] = pos+3;      triVal[sum] = -tmp1[0]; 
                    pJ1[9]  += -tmp1[1];                    //triI[sum] = i+1;              triJ[sum] = pos+3;      triVal[sum] = -tmp1[1]; 
                    pJ1[15] += -tmp1[2];                    //triI[sum] = i+2;              triJ[sum] = pos+3;      triVal[sum] = -tmp1[2]; 
                    pJ1[4]  += -tmp2[0];                    //triI[sum] = i;                        triJ[sum] = pos+4;      triVal[sum] = -tmp2[0]; 
                    pJ1[10] += -tmp2[1];                    //triI[sum] = i+1;              triJ[sum] = pos+4;      triVal[sum] = -tmp2[1]; 
                    pJ1[16] += -tmp2[2];                    //triI[sum] = i+2;              triJ[sum] = pos+4;      triVal[sum] = -tmp2[2];
                    pJ1[5]  += -tmp3[0];                    //triI[sum] = i;                        triJ[sum] = pos+5;      triVal[sum] = -tmp3[0];         
                    pJ1[11] += -tmp3[1];                    //triI[sum] = i+1;              triJ[sum] = pos+5;      triVal[sum] = -tmp3[1];
                    pJ1[17] += -tmp3[2];                    //triI[sum] = i+2;              triJ[sum] = pos+5;      triVal[sum] = -tmp3[2]; 
                    pJ1[21] += dA[0];                               //triI[sum] = i+3;              triJ[sum] = pos+3;      triVal[sum] = dA[0];    
                    pJ1[27] += dA[1];                               //triI[sum] = i+4;              triJ[sum] = pos+3;      triVal[sum] = dA[1];    
                    pJ1[33] += dA[2];                               //triI[sum] = i+5;              triJ[sum] = pos+3;      triVal[sum] = dA[2];    
                    pJ1[22] += dB[0];                               //triI[sum] = i+3;              triJ[sum] = pos+4;      triVal[sum] = dB[0];    
                    pJ1[28] += dB[1];                               //triI[sum] = i+4;              triJ[sum] = pos+4;      triVal[sum] = dB[1];    
                    pJ1[34] += dB[2];                               //triI[sum] = i+5;              triJ[sum] = pos+4;      triVal[sum] = dB[2];
                    pJ1[23] += dG[0];                               //triI[sum] = i+3;              triJ[sum] = pos+5;      triVal[sum] = dG[0];            
                    pJ1[29] += dG[1];                               //triI[sum] = i+4;              triJ[sum] = pos+5;      triVal[sum] = dG[1];
                    pJ1[35] += dG[2];                               //triI[sum] = i+5;              triJ[sum] = pos+5;      triVal[sum] = dG[2];    

				}
                else
				{
					t2[0] = ptr2[i];
					t2[1] = ptr2[i+1];
					t2[2] = ptr2[i+2];

					Alpha2 = ptr2[i+3];
					Beta2  = ptr2[i+4];
					Gamma2 = ptr2[i+5];

					lmj_Rderivation( Alpha2, Beta2, Gamma2, R2, dRA2, dRB2, dRG2 );

					lmj_TimesRRT( Ri, R2, R );

					double dRidA2[9], dRidB2[9], dRidG2[9];

					lmj_TimesRRT( dRidA2, dRA2, R );
					lmj_TimesRRT( dRidB2, dRB2, R );
					lmj_TimesRRT( dRidG2, dRG2, R );


					double dRidA[9], dRidB[9], dRidG[9];
					lmj_TimesRRT( dRidA, R2, dRA );
					lmj_TimesRRT( dRidB, R2, dRB );
					lmj_TimesRRT( dRidG, R2, dRG );


					lmj_dRi( ddA2, dRidA2, Ri );
					lmj_dRi( ddB2, dRidB2, Ri );
					lmj_dRi( ddG2, dRidG2, Ri );

					lmj_dRi( ddA, dRidA, Ri );
					lmj_dRi( ddB, dRidB, Ri );
					lmj_dRi( ddG, dRidG, Ri );

					tmp1[0] = dRA[0]*(t2[0]-t[0]) + dRA[1]*(t2[1]-t[1]) + dRA[2]*(t2[2]-t[2]);
					tmp1[1] = dRA[3]*(t2[0]-t[0]) + dRA[4]*(t2[1]-t[1]) + dRA[5]*(t2[2]-t[2]);
					tmp1[2] = dRA[6]*(t2[0]-t[0]) + dRA[7]*(t2[1]-t[1]) + dRA[8]*(t2[2]-t[2]);

					tmp2[0] = dRB[0]*(t2[0]-t[0]) + dRB[1]*(t2[1]-t[1]) + dRB[2]*(t2[2]-t[2]);
					tmp2[1] = dRB[3]*(t2[0]-t[0]) + dRB[4]*(t2[1]-t[1]) + dRB[5]*(t2[2]-t[2]);
					tmp2[2] = dRB[6]*(t2[0]-t[0]) + dRB[7]*(t2[1]-t[1]) + dRB[8]*(t2[2]-t[2]);

					tmp3[0] = dRG[0]*(t2[0]-t[0]) + dRG[1]*(t2[1]-t[1]) + dRG[2]*(t2[2]-t[2]);
					tmp3[1] = dRG[3]*(t2[0]-t[0]) + dRG[4]*(t2[1]-t[1]) + dRG[5]*(t2[2]-t[2]);
					tmp3[2] = dRG[6]*(t2[0]-t[0]) + dRG[7]*(t2[1]-t[1]) + dRG[8]*(t2[2]-t[2]);

					Id = i/6;
					pJ2 = J2 + Id*6*6;
                    pJ2[0]  += -R[0];                                       //triI[sum] = i;                        triJ[sum] = pos;        triVal[sum] = -R[0];    
                    pJ2[6]  += -R[3];                                       //triI[sum] = i+1;              triJ[sum] = pos;        triVal[sum] = -R[3];    
                    pJ2[12]  += -R[6];                                      //triI[sum] = i+2;              triJ[sum] = pos;        triVal[sum] = -R[6];    
                    pJ2[1]  += -R[1];                                       //triI[sum] = i;                        triJ[sum] = pos+1;      triVal[sum] = -R[1];    
                    pJ2[7]  += -R[4];                                       //triI[sum] = i+1;              triJ[sum] = pos+1;      triVal[sum] = -R[4];    
                    pJ2[13]  += -R[7];                                      //triI[sum] = i+2;              triJ[sum] = pos+1;      triVal[sum] = -R[7];
                    pJ2[2]  += -R[2];                                       //triI[sum] = i;                        triJ[sum] = pos+2;      triVal[sum] = -R[2];            
                    pJ2[8]  += -R[5];                                       //triI[sum] = i+1;              triJ[sum] = pos+2;      triVal[sum] = -R[5];
                    pJ2[14]  += -R[8];                                      //triI[sum] = i+2;              triJ[sum] = pos+2;      triVal[sum] = -R[8];    
                    pJ2[3]  += tmp1[0];                                     //triI[sum] = i;                        triJ[sum] = pos+3;      triVal[sum] = tmp1[0];  
                    pJ2[9]  += tmp1[1];                                     //triI[sum] = i+1;              triJ[sum] = pos+3;      triVal[sum] = tmp1[1];  
                    pJ2[15]  += tmp1[2];                            //triI[sum] = i+2;              triJ[sum] = pos+3;      triVal[sum] = tmp1[2];  
                    pJ2[4]  += tmp2[0];                                     //triI[sum] = i;                        triJ[sum] = pos+4;      triVal[sum] = tmp2[0];  
                    pJ2[10]  += tmp2[1];                            //triI[sum] = i+1;              triJ[sum] = pos+4;      triVal[sum] = tmp2[1];  
                    pJ2[16]  += tmp2[2];                            //triI[sum] = i+2;              triJ[sum] = pos+4;      triVal[sum] = tmp2[2];
                    pJ2[5]  += tmp3[0];                                     //triI[sum] = i;                        triJ[sum] = pos+5;      triVal[sum] = tmp3[0];          
                    pJ2[11]  += tmp3[1];                            //triI[sum] = i+1;              triJ[sum] = pos+5;      triVal[sum] = tmp3[1];
                    pJ2[17]  += tmp3[2];                            //triI[sum] = i+2;              triJ[sum] = pos+5;      triVal[sum] = tmp3[2];  
                    pJ2[21]  += ddA[0];                                     //triI[sum] = i+3;                      triJ[sum] = pos+3;      triVal[sum] = ddA[0];   
                    pJ2[27]  += ddA[1];                                     //triI[sum] = i+4;                      triJ[sum] = pos+3;      triVal[sum] = ddA[1];   
                    pJ2[33]  += ddA[2];                                     //triI[sum] = i+5;                      triJ[sum] = pos+3;      triVal[sum] = ddA[2];   
                    pJ2[22]  += ddB[0];                                     //triI[sum] = i+3;                      triJ[sum] = pos+4;      triVal[sum] = ddB[0];   
                    pJ2[28]  += ddB[1];                                     //triI[sum] = i+4;                      triJ[sum] = pos+4;      triVal[sum] = ddB[1];   
                    pJ2[34]  += ddB[2];                                     //triI[sum] = i+5;                      triJ[sum] = pos+4;      triVal[sum] = ddB[2];
                    pJ2[23]  += ddG[0];                                     //triI[sum] = i+3;                      triJ[sum] = pos+5;      triVal[sum] = ddG[0];           
                    pJ2[29]  += ddG[1];                                     //triI[sum] = i+4;                      triJ[sum] = pos+5;      triVal[sum] = ddG[1];
                    pJ2[35]  += ddG[2];                                     //triI[sum] = i+5;                      triJ[sum] = pos+5;      triVal[sum] = ddG[2];   

                    Id = i/6;
                    pJ1 = J1 + Id*6*6;
                    pJ1[0] += R[0];                                         //triI[sum] = i;                        triJ[sum] = i;          triVal[sum] = R[0];     
                    pJ1[6] += R[3];                                         //triI[sum] = i+1;              triJ[sum] = i;          triVal[sum] = R[3];     
                    pJ1[12] += R[6];                                        //triI[sum] = i+2;              triJ[sum] = i;          triVal[sum] = R[6];     
                    pJ1[1] += R[1];                                         //triI[sum] = i;                        triJ[sum] = i+1;        triVal[sum] = R[1];     
                    pJ1[7] += R[4];                                         //triI[sum] = i+1;              triJ[sum] = i+1;        triVal[sum] = R[4];     
                    pJ1[13] += R[7];                                        //triI[sum] = i+2;              triJ[sum] = i+1;        triVal[sum] = R[7];
                    pJ1[2] += R[2];                                         //triI[sum] = i;                        triJ[sum] = i+2;        triVal[sum] = R[2];             
                    pJ1[8] += R[5];                                         //triI[sum] = i+1;              triJ[sum] = i+2;        triVal[sum] = R[5];
                    pJ1[14] += R[8];                                        //triI[sum] = i+2;              triJ[sum] = i+2;        triVal[sum] = R[8];     
                    pJ1[21] += ddA2[0];                                     //triI[sum] = i+3;                      triJ[sum] = i+3;        triVal[sum] = ddA2[0];  
                    pJ1[27] += ddA2[1];                                     //triI[sum] = i+4;                      triJ[sum] = i+3;        triVal[sum] = ddA2[1];  
                    pJ1[33] += ddA2[2];                                     //triI[sum] = i+5;                      triJ[sum] = i+3;        triVal[sum] = ddA2[2];  
                    pJ1[22] += ddB2[0];                                     //triI[sum] = i+3;                      triJ[sum] = i+4;        triVal[sum] = ddB2[0];  
                    pJ1[28] += ddB2[1];                                     //triI[sum] = i+4;                      triJ[sum] = i+4;        triVal[sum] = ddB2[1];  
                    pJ1[34] += ddB2[2];                                     //triI[sum] = i+5;                      triJ[sum] = i+4;        triVal[sum] = ddB2[2];
                    pJ1[23] += ddG2[0];                                     //triI[sum] = i+3;                      triJ[sum] = i+5;        triVal[sum] = ddG2[0];          
                    pJ1[29] += ddG2[1];                                     //triI[sum] = i+4;                      triJ[sum] = i+5;        triVal[sum] = ddG2[1];
                    pJ1[35] += ddG2[2];                                     //triI[sum] = i+5;                      triJ[sum] = i+5;        triVal[sum] = ddG2[2];  
				}
                           i += 5;
			}
            else
			{
				tmp1[0] = dRA[0]*(ptr2[i]-t[0])+dRA[1]*(ptr2[i+1]-t[1])+dRA[2]*(ptr2[i+2]-t[2]);
				tmp1[1] = dRA[3]*(ptr2[i]-t[0])+dRA[4]*(ptr2[i+1]-t[1])+dRA[5]*(ptr2[i+2]-t[2]);
				tmp1[2] = dRA[6]*(ptr2[i]-t[0])+dRA[7]*(ptr2[i+1]-t[1])+dRA[8]*(ptr2[i+2]-t[2]);

				tmp2[0] = dRB[0]*(ptr2[i]-t[0])+dRB[1]*(ptr2[i+1]-t[1])+dRB[2]*(ptr2[i+2]-t[2]);
				tmp2[1] = dRB[3]*(ptr2[i]-t[0])+dRB[4]*(ptr2[i+1]-t[1])+dRB[5]*(ptr2[i+2]-t[2]);
				tmp2[2] = dRB[6]*(ptr2[i]-t[0])+dRB[7]*(ptr2[i+1]-t[1])+dRB[8]*(ptr2[i+2]-t[2]);

				tmp3[0] = dRG[0]*(ptr2[i]-t[0])+dRG[1]*(ptr2[i+1]-t[1])+dRG[2]*(ptr2[i+2]-t[2]);
				tmp3[1] = dRG[3]*(ptr2[i]-t[0])+dRG[4]*(ptr2[i+1]-t[1])+dRG[5]*(ptr2[i+2]-t[2]);
				tmp3[2] = dRG[6]*(ptr2[i]-t[0])+dRG[7]*(ptr2[i+1]-t[1])+dRG[8]*(ptr2[i+2]-t[2]);

				Id = (i - m*6)/3;
				pJ2 = J2 + m*6*6 + Id*3*6;
				pJ2[0]  += -R[0];                               //triI[sum] = i;                        triJ[sum] = pos;        triVal[sum] = -R[0];    
				pJ2[6]  += -R[3];                               //triI[sum] = i+1;              triJ[sum] = pos;        triVal[sum] = -R[3];    
				pJ2[12]  += -R[6];                              //triI[sum] = i+2;              triJ[sum] = pos;        triVal[sum] = -R[6];    
				pJ2[1]  += -R[1];                               //triI[sum] = i;                        triJ[sum] = pos+1;      triVal[sum] = -R[1];    
				pJ2[7]  += -R[4];                               //triI[sum] = i+1;              triJ[sum] = pos+1;      triVal[sum] = -R[4];    
				pJ2[13]  += -R[7];                              //triI[sum] = i+2;              triJ[sum] = pos+1;      triVal[sum] = -R[7];
				pJ2[2]  += -R[2];                               //triI[sum] = i;                        triJ[sum] = pos+2;      triVal[sum] = -R[2];            
				pJ2[8]  += -R[5];                               //triI[sum] = i+1;              triJ[sum] = pos+2;      triVal[sum] = -R[5];
				pJ2[14]  += -R[8];                              //triI[sum] = i+2;              triJ[sum] = pos+2;      triVal[sum] = -R[8];    
				pJ2[3]  += tmp1[0];                             //triI[sum] = i;                        triJ[sum] = pos+3;      triVal[sum] = tmp1[0];  
				pJ2[9]  += tmp1[1];                             //triI[sum] = i+1;              triJ[sum] = pos+3;      triVal[sum] = tmp1[1];  
				pJ2[15]  += tmp1[2];                    //triI[sum] = i+2;              triJ[sum] = pos+3;      triVal[sum] = tmp1[2];  
				pJ2[4]  += tmp2[0];                             //triI[sum] = i;                        triJ[sum] = pos+4;      triVal[sum] = tmp2[0];  
				pJ2[10]  += tmp2[1];                    //triI[sum] = i+1;              triJ[sum] = pos+4;      triVal[sum] = tmp2[1];  
				pJ2[16]  += tmp2[2];                    //triI[sum] = i+2;              triJ[sum] = pos+4;      triVal[sum] = tmp2[2];
				pJ2[5]  += tmp3[0];                             //triI[sum] = i;                        triJ[sum] = pos+5;      triVal[sum] = tmp3[0];          
				pJ2[11]  += tmp3[1];                    //triI[sum] = i+1;              triJ[sum] = pos+5;      triVal[sum] = tmp3[1];
				pJ2[17]  += tmp3[2];                    //triI[sum] = i+2;              triJ[sum] = pos+5;      triVal[sum] = tmp3[2];  

				pJ1 = J1 + m*6*6 + Id*3*3;
                pJ1[0]  += R[0];                                //triI[sum] = i;                        triJ[sum] = i;          triVal[sum] = R[0];     
                pJ1[3]  += R[3];                                //triI[sum] = i+1;              triJ[sum] = i;          triVal[sum] = R[3];     
                pJ1[6]  += R[6];                                //triI[sum] = i+2;              triJ[sum] = i;          triVal[sum] = R[6];     
                pJ1[1]  += R[1];                                //triI[sum] = i;                        triJ[sum] = i+1;        triVal[sum] = R[1];     
                pJ1[4]  += R[4];                                //triI[sum] = i+1;              triJ[sum] = i+1;        triVal[sum] = R[4];     
                pJ1[7]  += R[7];                                //triI[sum] = i+2;              triJ[sum] = i+1;        triVal[sum] = R[7];
                pJ1[2]  += R[2];                                //triI[sum] = i;                        triJ[sum] = i+2;        triVal[sum] = R[2];             
                pJ1[5]  += R[5];                                //triI[sum] = i+1;              triJ[sum] = i+2;        triVal[sum] = R[5];
                pJ1[8]  += R[8];                                //triI[sum] = i+2;              triJ[sum] = i+2;        triVal[sum] = R[8];     
                i += 2;
			}
		}
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                //================================================================
                // Compute newU, the number of elements in newU, as well as the ID of newU
                double* ptrU, *ptrnU; // U 和 newU 指针

                double* pJ1i, *pJ1j, *pJ2i, *pJ2j, *ptmp;

                

                int sizenU, n_newU;// n_newU is the real number of elements in newU, it should be smaller than sizenU 
                sizenU = (m_nU+m);
                
				GMap_End.U = (double*)malloc( sizenU*6*6 * sizeof(double) ); 
				double*  newU = GMap_End.U;
                memset( newU, 0, sizenU*6*6 * sizeof(double) );
				GMap_End.Ui = (int*)malloc( sizenU * sizeof(int) );
                int*  m_nUi = GMap_End.Ui;
				GMap_End.Uj = (int*)malloc( sizenU * sizeof(int) );
                int*  m_nUj =GMap_End.Uj;

                n_newU = m;
                int posID = pos/6;

                ptrU = m_U;
                ptrnU = newU+m*(6*6);

                for ( int i = 0; i < m; i++ )
                {
                        if ( i <= posID)
                        {
                                m_nUi[i] = i;
                                m_nUj[i] = posID;
                        }
                        else
                        {
                                m_nUi[i] = posID;
                                m_nUj[i] = i;
                        }
				}
                
				for ( int i = 0; i < m_nU; i++ )
                {
                        ptrU = m_U + i*6*6;

                        pJ1i = J1 + m_Ui[i]*6*6;
                        pJ1j = J1 + m_Uj[i]*6*6;
                        pJ2i = J2 + m_Ui[i]*6*6;
                        pJ2j = J2 + m_Uj[i]*6*6;

                        //Algorithm Line 4
                        tmp1[0] = pJ2i[0]*ptrU[0]+pJ2i[6]*ptrU[6]+pJ2i[12]*ptrU[12]+pJ2i[18]*ptrU[18]+pJ2i[24]*ptrU[24]+pJ2i[30]*ptrU[30];
                        tmp1[1] = pJ2i[0]*ptrU[1]+pJ2i[6]*ptrU[7]+pJ2i[12]*ptrU[13]+pJ2i[18]*ptrU[19]+pJ2i[24]*ptrU[25]+pJ2i[30]*ptrU[31];
                        tmp1[2] = pJ2i[0]*ptrU[2]+pJ2i[6]*ptrU[8]+pJ2i[12]*ptrU[14]+pJ2i[18]*ptrU[20]+pJ2i[24]*ptrU[26]+pJ2i[30]*ptrU[32];
                        tmp1[3] = pJ2i[0]*ptrU[3]+pJ2i[6]*ptrU[9]+pJ2i[12]*ptrU[15]+pJ2i[18]*ptrU[21]+pJ2i[24]*ptrU[27]+pJ2i[30]*ptrU[33];
                        tmp1[4] = pJ2i[0]*ptrU[4]+pJ2i[6]*ptrU[10]+pJ2i[12]*ptrU[16]+pJ2i[18]*ptrU[22]+pJ2i[24]*ptrU[28]+pJ2i[30]*ptrU[34];
                        tmp1[5] = pJ2i[0]*ptrU[5]+pJ2i[6]*ptrU[11]+pJ2i[12]*ptrU[17]+pJ2i[18]*ptrU[23]+pJ2i[24]*ptrU[29]+pJ2i[30]*ptrU[35];

                        tmp2[0] = pJ2i[1]*ptrU[0]+pJ2i[7]*ptrU[6]+pJ2i[13]*ptrU[12]+pJ2i[19]*ptrU[18]+pJ2i[25]*ptrU[24]+pJ2i[31]*ptrU[30];
                        tmp2[1] = pJ2i[1]*ptrU[1]+pJ2i[7]*ptrU[7]+pJ2i[13]*ptrU[13]+pJ2i[19]*ptrU[19]+pJ2i[25]*ptrU[25]+pJ2i[31]*ptrU[31];
                        tmp2[2] = pJ2i[1]*ptrU[2]+pJ2i[7]*ptrU[8]+pJ2i[13]*ptrU[14]+pJ2i[19]*ptrU[20]+pJ2i[25]*ptrU[26]+pJ2i[31]*ptrU[32];
                        tmp2[3] = pJ2i[1]*ptrU[3]+pJ2i[7]*ptrU[9]+pJ2i[13]*ptrU[15]+pJ2i[19]*ptrU[21]+pJ2i[25]*ptrU[27]+pJ2i[31]*ptrU[33];
                        tmp2[4] = pJ2i[1]*ptrU[4]+pJ2i[7]*ptrU[10]+pJ2i[13]*ptrU[16]+pJ2i[19]*ptrU[22]+pJ2i[25]*ptrU[28]+pJ2i[31]*ptrU[34];
                        tmp2[5] = pJ2i[1]*ptrU[5]+pJ2i[7]*ptrU[11]+pJ2i[13]*ptrU[17]+pJ2i[19]*ptrU[23]+pJ2i[25]*ptrU[29]+pJ2i[31]*ptrU[35];

                        tmp3[0] = pJ2i[2]*ptrU[0]+pJ2i[8]*ptrU[6]+pJ2i[14]*ptrU[12]+pJ2i[20]*ptrU[18]+pJ2i[26]*ptrU[24]+pJ2i[32]*ptrU[30];
                        tmp3[1] = pJ2i[2]*ptrU[1]+pJ2i[8]*ptrU[7]+pJ2i[14]*ptrU[13]+pJ2i[20]*ptrU[19]+pJ2i[26]*ptrU[25]+pJ2i[32]*ptrU[31];
                        tmp3[2] = pJ2i[2]*ptrU[2]+pJ2i[8]*ptrU[8]+pJ2i[14]*ptrU[14]+pJ2i[20]*ptrU[20]+pJ2i[26]*ptrU[26]+pJ2i[32]*ptrU[32];
                        tmp3[3] = pJ2i[2]*ptrU[3]+pJ2i[8]*ptrU[9]+pJ2i[14]*ptrU[15]+pJ2i[20]*ptrU[21]+pJ2i[26]*ptrU[27]+pJ2i[32]*ptrU[33];
                        tmp3[4] = pJ2i[2]*ptrU[4]+pJ2i[8]*ptrU[10]+pJ2i[14]*ptrU[16]+pJ2i[20]*ptrU[22]+pJ2i[26]*ptrU[28]+pJ2i[32]*ptrU[34];
                        tmp3[5] = pJ2i[2]*ptrU[5]+pJ2i[8]*ptrU[11]+pJ2i[14]*ptrU[17]+pJ2i[20]*ptrU[23]+pJ2i[26]*ptrU[29]+pJ2i[32]*ptrU[35];

                        tmp4[0] = pJ2i[3]*ptrU[0]+pJ2i[9]*ptrU[6]+pJ2i[15]*ptrU[12]+pJ2i[21]*ptrU[18]+pJ2i[27]*ptrU[24]+pJ2i[33]*ptrU[30];
                        tmp4[1] = pJ2i[3]*ptrU[1]+pJ2i[9]*ptrU[7]+pJ2i[15]*ptrU[13]+pJ2i[21]*ptrU[19]+pJ2i[27]*ptrU[25]+pJ2i[33]*ptrU[31];
                        tmp4[2] = pJ2i[3]*ptrU[2]+pJ2i[9]*ptrU[8]+pJ2i[15]*ptrU[14]+pJ2i[21]*ptrU[20]+pJ2i[27]*ptrU[26]+pJ2i[33]*ptrU[32];
                        tmp4[3] = pJ2i[3]*ptrU[3]+pJ2i[9]*ptrU[9]+pJ2i[15]*ptrU[15]+pJ2i[21]*ptrU[21]+pJ2i[27]*ptrU[27]+pJ2i[33]*ptrU[33];
                        tmp4[4] = pJ2i[3]*ptrU[4]+pJ2i[9]*ptrU[10]+pJ2i[15]*ptrU[16]+pJ2i[21]*ptrU[22]+pJ2i[27]*ptrU[28]+pJ2i[33]*ptrU[34];
                        tmp4[5] = pJ2i[3]*ptrU[5]+pJ2i[9]*ptrU[11]+pJ2i[15]*ptrU[17]+pJ2i[21]*ptrU[23]+pJ2i[27]*ptrU[29]+pJ2i[33]*ptrU[35];

                        tmp5[0] = pJ2i[4]*ptrU[0]+pJ2i[10]*ptrU[6]+pJ2i[16]*ptrU[12]+pJ2i[22]*ptrU[18]+pJ2i[28]*ptrU[24]+pJ2i[34]*ptrU[30];
                        tmp5[1] = pJ2i[4]*ptrU[1]+pJ2i[10]*ptrU[7]+pJ2i[16]*ptrU[13]+pJ2i[22]*ptrU[19]+pJ2i[28]*ptrU[25]+pJ2i[34]*ptrU[31];
                        tmp5[2] = pJ2i[4]*ptrU[2]+pJ2i[10]*ptrU[8]+pJ2i[16]*ptrU[14]+pJ2i[22]*ptrU[20]+pJ2i[28]*ptrU[26]+pJ2i[34]*ptrU[32];
                        tmp5[3] = pJ2i[4]*ptrU[3]+pJ2i[10]*ptrU[9]+pJ2i[16]*ptrU[15]+pJ2i[22]*ptrU[21]+pJ2i[28]*ptrU[27]+pJ2i[34]*ptrU[33];
                        tmp5[4] = pJ2i[4]*ptrU[4]+pJ2i[10]*ptrU[10]+pJ2i[16]*ptrU[16]+pJ2i[22]*ptrU[22]+pJ2i[28]*ptrU[28]+pJ2i[34]*ptrU[34];
                        tmp5[5] = pJ2i[4]*ptrU[5]+pJ2i[10]*ptrU[11]+pJ2i[16]*ptrU[17]+pJ2i[22]*ptrU[23]+pJ2i[28]*ptrU[29]+pJ2i[34]*ptrU[35];

                        tmp6[0] = pJ2i[5]*ptrU[0]+pJ2i[11]*ptrU[6]+pJ2i[17]*ptrU[12]+pJ2i[23]*ptrU[18]+pJ2i[29]*ptrU[24]+pJ2i[35]*ptrU[30];
                        tmp6[1] = pJ2i[5]*ptrU[1]+pJ2i[11]*ptrU[7]+pJ2i[17]*ptrU[13]+pJ2i[23]*ptrU[19]+pJ2i[29]*ptrU[25]+pJ2i[35]*ptrU[31];
                        tmp6[2] = pJ2i[5]*ptrU[2]+pJ2i[11]*ptrU[8]+pJ2i[17]*ptrU[14]+pJ2i[23]*ptrU[20]+pJ2i[29]*ptrU[26]+pJ2i[35]*ptrU[32];
                        tmp6[3] = pJ2i[5]*ptrU[3]+pJ2i[11]*ptrU[9]+pJ2i[17]*ptrU[15]+pJ2i[23]*ptrU[21]+pJ2i[29]*ptrU[27]+pJ2i[35]*ptrU[33];
                        tmp6[4] = pJ2i[5]*ptrU[4]+pJ2i[11]*ptrU[10]+pJ2i[17]*ptrU[16]+pJ2i[23]*ptrU[22]+pJ2i[29]*ptrU[28]+pJ2i[35]*ptrU[34];
                        tmp6[5] = pJ2i[5]*ptrU[5]+pJ2i[11]*ptrU[11]+pJ2i[17]*ptrU[17]+pJ2i[23]*ptrU[23]+pJ2i[29]*ptrU[29]+pJ2i[35]*ptrU[35];

                        ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12]+tmp1[3]*pJ2j[18]+tmp1[4]*pJ2j[24]+tmp1[5]*pJ2j[30];
                        ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13]+tmp1[3]*pJ2j[19]+tmp1[4]*pJ2j[25]+tmp1[5]*pJ2j[31];
                        ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14]+tmp1[3]*pJ2j[20]+tmp1[4]*pJ2j[26]+tmp1[5]*pJ2j[32];
                        ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15]+tmp1[3]*pJ2j[21]+tmp1[4]*pJ2j[27]+tmp1[5]*pJ2j[33];
                        ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16]+tmp1[3]*pJ2j[22]+tmp1[4]*pJ2j[28]+tmp1[5]*pJ2j[34];
                        ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17]+tmp1[3]*pJ2j[23]+tmp1[4]*pJ2j[29]+tmp1[5]*pJ2j[35];

                        ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12]+tmp2[3]*pJ2j[18]+tmp2[4]*pJ2j[24]+tmp2[5]*pJ2j[30];
                        ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13]+tmp2[3]*pJ2j[19]+tmp2[4]*pJ2j[25]+tmp2[5]*pJ2j[31];
                        ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14]+tmp2[3]*pJ2j[20]+tmp2[4]*pJ2j[26]+tmp2[5]*pJ2j[32];
                        ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15]+tmp2[3]*pJ2j[21]+tmp2[4]*pJ2j[27]+tmp2[5]*pJ2j[33];
                        ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16]+tmp2[3]*pJ2j[22]+tmp2[4]*pJ2j[28]+tmp2[5]*pJ2j[34];
                        ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17]+tmp2[3]*pJ2j[23]+tmp2[4]*pJ2j[29]+tmp2[5]*pJ2j[35];

                        ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12]+tmp3[3]*pJ2j[18]+tmp3[4]*pJ2j[24]+tmp3[5]*pJ2j[30];
                        ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13]+tmp3[3]*pJ2j[19]+tmp3[4]*pJ2j[25]+tmp3[5]*pJ2j[31];
                        ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14]+tmp3[3]*pJ2j[20]+tmp3[4]*pJ2j[26]+tmp3[5]*pJ2j[32];
                        ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15]+tmp3[3]*pJ2j[21]+tmp3[4]*pJ2j[27]+tmp3[5]*pJ2j[33];
                        ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16]+tmp3[3]*pJ2j[22]+tmp3[4]*pJ2j[28]+tmp3[5]*pJ2j[34];
                        ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17]+tmp3[3]*pJ2j[23]+tmp3[4]*pJ2j[29]+tmp3[5]*pJ2j[35];

                        ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12]+tmp4[3]*pJ2j[18]+tmp4[4]*pJ2j[24]+tmp4[5]*pJ2j[30];
                        ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13]+tmp4[3]*pJ2j[19]+tmp4[4]*pJ2j[25]+tmp4[5]*pJ2j[31];
                        ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14]+tmp4[3]*pJ2j[20]+tmp4[4]*pJ2j[26]+tmp4[5]*pJ2j[32];
                        ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15]+tmp4[3]*pJ2j[21]+tmp4[4]*pJ2j[27]+tmp4[5]*pJ2j[33];
                        ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16]+tmp4[3]*pJ2j[22]+tmp4[4]*pJ2j[28]+tmp4[5]*pJ2j[34];
                        ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17]+tmp4[3]*pJ2j[23]+tmp4[4]*pJ2j[29]+tmp4[5]*pJ2j[35];

                        ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12]+tmp5[3]*pJ2j[18]+tmp5[4]*pJ2j[24]+tmp5[5]*pJ2j[30];
                        ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13]+tmp5[3]*pJ2j[19]+tmp5[4]*pJ2j[25]+tmp5[5]*pJ2j[31];
                        ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14]+tmp5[3]*pJ2j[20]+tmp5[4]*pJ2j[26]+tmp5[5]*pJ2j[32];
                        ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15]+tmp5[3]*pJ2j[21]+tmp5[4]*pJ2j[27]+tmp5[5]*pJ2j[33];
                        ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16]+tmp5[3]*pJ2j[22]+tmp5[4]*pJ2j[28]+tmp5[5]*pJ2j[34];
                        ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17]+tmp5[3]*pJ2j[23]+tmp5[4]*pJ2j[29]+tmp5[5]*pJ2j[35];

                        ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12]+tmp6[3]*pJ2j[18]+tmp6[4]*pJ2j[24]+tmp6[5]*pJ2j[30];
                        ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13]+tmp6[3]*pJ2j[19]+tmp6[4]*pJ2j[25]+tmp6[5]*pJ2j[31];
                        ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14]+tmp6[3]*pJ2j[20]+tmp6[4]*pJ2j[26]+tmp6[5]*pJ2j[32];
                        ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15]+tmp6[3]*pJ2j[21]+tmp6[4]*pJ2j[27]+tmp6[5]*pJ2j[33];
                        ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16]+tmp6[3]*pJ2j[22]+tmp6[4]*pJ2j[28]+tmp6[5]*pJ2j[34];
                        ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17]+tmp6[3]*pJ2j[23]+tmp6[4]*pJ2j[29]+tmp6[5]*pJ2j[35];

                        ptmp = newU+posID*(6*6);

                        // ttmp = J2(i)^T*I(j,k)*J2(j)
                        ptmp[0] += ttmp1[0];
                        ptmp[1] += ttmp1[1];
                        ptmp[2] += ttmp1[2];
                        ptmp[3] += ttmp1[3];
                        ptmp[4] += ttmp1[4];
                        ptmp[5] += ttmp1[5];
                        ptmp[6] += ttmp2[0];
                        ptmp[7] += ttmp2[1];
                        ptmp[8] += ttmp2[2];
                        ptmp[9] += ttmp2[3];
                        ptmp[10] += ttmp2[4];
                        ptmp[11] += ttmp2[5];
                        ptmp[12] += ttmp3[0];
                        ptmp[13] += ttmp3[1];
                        ptmp[14] += ttmp3[2];
                        ptmp[15] += ttmp3[3];
                        ptmp[16] += ttmp3[4];
                        ptmp[17] += ttmp3[5];
                        ptmp[18] += ttmp4[0];
                        ptmp[19] += ttmp4[1];
                        ptmp[20] += ttmp4[2];
                        ptmp[21] += ttmp4[3];
                        ptmp[22] += ttmp4[4];
                        ptmp[23] += ttmp4[5];
                        ptmp[24] += ttmp5[0];
                        ptmp[25] += ttmp5[1];
                        ptmp[26] += ttmp5[2];
                        ptmp[27] += ttmp5[3];
                        ptmp[28] += ttmp5[4];
                        ptmp[29] += ttmp5[5];
                        ptmp[30] += ttmp6[0];
                        ptmp[31] += ttmp6[1];
                        ptmp[32] += ttmp6[2];
                        ptmp[33] += ttmp6[3];
                        ptmp[34] += ttmp6[4];
                        ptmp[35] += ttmp6[5];                                   


                        if ( m_Ui[i] != m_Uj[i] )
                        {
                                //ttmp = ttmp^T
                                
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];           

                        }


                        //Algorithm Line 5

                        ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[6]+tmp1[2]*pJ1j[12]+tmp1[3]*pJ1j[18]+tmp1[4]*pJ1j[24]+tmp1[5]*pJ1j[30];
                        ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[7]+tmp1[2]*pJ1j[13]+tmp1[3]*pJ1j[19]+tmp1[4]*pJ1j[25]+tmp1[5]*pJ1j[31];
                        ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[8]+tmp1[2]*pJ1j[14]+tmp1[3]*pJ1j[20]+tmp1[4]*pJ1j[26]+tmp1[5]*pJ1j[32];
                        ttmp1[3] = tmp1[0]*pJ1j[3]+tmp1[1]*pJ1j[9]+tmp1[2]*pJ1j[15]+tmp1[3]*pJ1j[21]+tmp1[4]*pJ1j[27]+tmp1[5]*pJ1j[33];
                        ttmp1[4] = tmp1[0]*pJ1j[4]+tmp1[1]*pJ1j[10]+tmp1[2]*pJ1j[16]+tmp1[3]*pJ1j[22]+tmp1[4]*pJ1j[28]+tmp1[5]*pJ1j[34];
                        ttmp1[5] = tmp1[0]*pJ1j[5]+tmp1[1]*pJ1j[11]+tmp1[2]*pJ1j[17]+tmp1[3]*pJ1j[23]+tmp1[4]*pJ1j[29]+tmp1[5]*pJ1j[35];

                        ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[6]+tmp2[2]*pJ1j[12]+tmp2[3]*pJ1j[18]+tmp2[4]*pJ1j[24]+tmp2[5]*pJ1j[30];
                        ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[7]+tmp2[2]*pJ1j[13]+tmp2[3]*pJ1j[19]+tmp2[4]*pJ1j[25]+tmp2[5]*pJ1j[31];
                        ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[8]+tmp2[2]*pJ1j[14]+tmp2[3]*pJ1j[20]+tmp2[4]*pJ1j[26]+tmp2[5]*pJ1j[32];
                        ttmp2[3] = tmp2[0]*pJ1j[3]+tmp2[1]*pJ1j[9]+tmp2[2]*pJ1j[15]+tmp2[3]*pJ1j[21]+tmp2[4]*pJ1j[27]+tmp2[5]*pJ1j[33];
                        ttmp2[4] = tmp2[0]*pJ1j[4]+tmp2[1]*pJ1j[10]+tmp2[2]*pJ1j[16]+tmp2[3]*pJ1j[22]+tmp2[4]*pJ1j[28]+tmp2[5]*pJ1j[34];
                        ttmp2[5] = tmp2[0]*pJ1j[5]+tmp2[1]*pJ1j[11]+tmp2[2]*pJ1j[17]+tmp2[3]*pJ1j[23]+tmp2[4]*pJ1j[29]+tmp2[5]*pJ1j[35];

                        ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[6]+tmp3[2]*pJ1j[12]+tmp3[3]*pJ1j[18]+tmp3[4]*pJ1j[24]+tmp3[5]*pJ1j[30];
                        ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[7]+tmp3[2]*pJ1j[13]+tmp3[3]*pJ1j[19]+tmp3[4]*pJ1j[25]+tmp3[5]*pJ1j[31];
                        ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[8]+tmp3[2]*pJ1j[14]+tmp3[3]*pJ1j[20]+tmp3[4]*pJ1j[26]+tmp3[5]*pJ1j[32];
                        ttmp3[3] = tmp3[0]*pJ1j[3]+tmp3[1]*pJ1j[9]+tmp3[2]*pJ1j[15]+tmp3[3]*pJ1j[21]+tmp3[4]*pJ1j[27]+tmp3[5]*pJ1j[33];
                        ttmp3[4] = tmp3[0]*pJ1j[4]+tmp3[1]*pJ1j[10]+tmp3[2]*pJ1j[16]+tmp3[3]*pJ1j[22]+tmp3[4]*pJ1j[28]+tmp3[5]*pJ1j[34];
                        ttmp3[5] = tmp3[0]*pJ1j[5]+tmp3[1]*pJ1j[11]+tmp3[2]*pJ1j[17]+tmp3[3]*pJ1j[23]+tmp3[4]*pJ1j[29]+tmp3[5]*pJ1j[35];

                        ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[6]+tmp4[2]*pJ1j[12]+tmp4[3]*pJ1j[18]+tmp4[4]*pJ1j[24]+tmp4[5]*pJ1j[30];
                        ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[7]+tmp4[2]*pJ1j[13]+tmp4[3]*pJ1j[19]+tmp4[4]*pJ1j[25]+tmp4[5]*pJ1j[31];
                        ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[8]+tmp4[2]*pJ1j[14]+tmp4[3]*pJ1j[20]+tmp4[4]*pJ1j[26]+tmp4[5]*pJ1j[32];
                        ttmp4[3] = tmp4[0]*pJ1j[3]+tmp4[1]*pJ1j[9]+tmp4[2]*pJ1j[15]+tmp4[3]*pJ1j[21]+tmp4[4]*pJ1j[27]+tmp4[5]*pJ1j[33];
                        ttmp4[4] = tmp4[0]*pJ1j[4]+tmp4[1]*pJ1j[10]+tmp4[2]*pJ1j[16]+tmp4[3]*pJ1j[22]+tmp4[4]*pJ1j[28]+tmp4[5]*pJ1j[34];
                        ttmp4[5] = tmp4[0]*pJ1j[5]+tmp4[1]*pJ1j[11]+tmp4[2]*pJ1j[17]+tmp4[3]*pJ1j[23]+tmp4[4]*pJ1j[29]+tmp4[5]*pJ1j[35];

                        ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[6]+tmp5[2]*pJ1j[12]+tmp5[3]*pJ1j[18]+tmp5[4]*pJ1j[24]+tmp5[5]*pJ1j[30];
                        ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[7]+tmp5[2]*pJ1j[13]+tmp5[3]*pJ1j[19]+tmp5[4]*pJ1j[25]+tmp5[5]*pJ1j[31];
                        ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[8]+tmp5[2]*pJ1j[14]+tmp5[3]*pJ1j[20]+tmp5[4]*pJ1j[26]+tmp5[5]*pJ1j[32];
                        ttmp5[3] = tmp5[0]*pJ1j[3]+tmp5[1]*pJ1j[9]+tmp5[2]*pJ1j[15]+tmp5[3]*pJ1j[21]+tmp5[4]*pJ1j[27]+tmp5[5]*pJ1j[33];
                        ttmp5[4] = tmp5[0]*pJ1j[4]+tmp5[1]*pJ1j[10]+tmp5[2]*pJ1j[16]+tmp5[3]*pJ1j[22]+tmp5[4]*pJ1j[28]+tmp5[5]*pJ1j[34];
                        ttmp5[5] = tmp5[0]*pJ1j[5]+tmp5[1]*pJ1j[11]+tmp5[2]*pJ1j[17]+tmp5[3]*pJ1j[23]+tmp5[4]*pJ1j[29]+tmp5[5]*pJ1j[35];

                        ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[6]+tmp6[2]*pJ1j[12]+tmp6[3]*pJ1j[18]+tmp6[4]*pJ1j[24]+tmp6[5]*pJ1j[30];
                        ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[7]+tmp6[2]*pJ1j[13]+tmp6[3]*pJ1j[19]+tmp6[4]*pJ1j[25]+tmp6[5]*pJ1j[31];
                        ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[8]+tmp6[2]*pJ1j[14]+tmp6[3]*pJ1j[20]+tmp6[4]*pJ1j[26]+tmp6[5]*pJ1j[32];
                        ttmp6[3] = tmp6[0]*pJ1j[3]+tmp6[1]*pJ1j[9]+tmp6[2]*pJ1j[15]+tmp6[3]*pJ1j[21]+tmp6[4]*pJ1j[27]+tmp6[5]*pJ1j[33];
                        ttmp6[4] = tmp6[0]*pJ1j[4]+tmp6[1]*pJ1j[10]+tmp6[2]*pJ1j[16]+tmp6[3]*pJ1j[22]+tmp6[4]*pJ1j[28]+tmp6[5]*pJ1j[34];
                        ttmp6[5] = tmp6[0]*pJ1j[5]+tmp6[1]*pJ1j[11]+tmp6[2]*pJ1j[17]+tmp6[3]*pJ1j[23]+tmp6[4]*pJ1j[29]+tmp6[5]*pJ1j[35];

                        ptmp = newU+m_Uj[i]*(6*6);
                        if ( m_Uj[i] >= posID )
                        {
                                // ttmp = J2(i)^T*I(j,k)*J1(j)
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];
                        }
                        if ( m_Uj[i] <= posID && m_Ui[i] != m_Uj[i] )
                        {
                                //if ( m_Ui[i] != m_Uj[i] )
                                //{
                                // ttmp = ttmp^T
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                                //}


                        }

                        //Algorithm Line 3

                        tmp1[0] = pJ1i[0]*ptrU[0]+pJ1i[6]*ptrU[6]+pJ1i[12]*ptrU[12]+pJ1i[18]*ptrU[18]+pJ1i[24]*ptrU[24]+pJ1i[30]*ptrU[30];
                        tmp1[1] = pJ1i[0]*ptrU[1]+pJ1i[6]*ptrU[7]+pJ1i[12]*ptrU[13]+pJ1i[18]*ptrU[19]+pJ1i[24]*ptrU[25]+pJ1i[30]*ptrU[31];
                        tmp1[2] = pJ1i[0]*ptrU[2]+pJ1i[6]*ptrU[8]+pJ1i[12]*ptrU[14]+pJ1i[18]*ptrU[20]+pJ1i[24]*ptrU[26]+pJ1i[30]*ptrU[32];
                        tmp1[3] = pJ1i[0]*ptrU[3]+pJ1i[6]*ptrU[9]+pJ1i[12]*ptrU[15]+pJ1i[18]*ptrU[21]+pJ1i[24]*ptrU[27]+pJ1i[30]*ptrU[33];
                        tmp1[4] = pJ1i[0]*ptrU[4]+pJ1i[6]*ptrU[10]+pJ1i[12]*ptrU[16]+pJ1i[18]*ptrU[22]+pJ1i[24]*ptrU[28]+pJ1i[30]*ptrU[34];
                        tmp1[5] = pJ1i[0]*ptrU[5]+pJ1i[6]*ptrU[11]+pJ1i[12]*ptrU[17]+pJ1i[18]*ptrU[23]+pJ1i[24]*ptrU[29]+pJ1i[30]*ptrU[35];

                        tmp2[0] = pJ1i[1]*ptrU[0]+pJ1i[7]*ptrU[6]+pJ1i[13]*ptrU[12]+pJ1i[19]*ptrU[18]+pJ1i[25]*ptrU[24]+pJ1i[31]*ptrU[30];
                        tmp2[1] = pJ1i[1]*ptrU[1]+pJ1i[7]*ptrU[7]+pJ1i[13]*ptrU[13]+pJ1i[19]*ptrU[19]+pJ1i[25]*ptrU[25]+pJ1i[31]*ptrU[31];
                        tmp2[2] = pJ1i[1]*ptrU[2]+pJ1i[7]*ptrU[8]+pJ1i[13]*ptrU[14]+pJ1i[19]*ptrU[20]+pJ1i[25]*ptrU[26]+pJ1i[31]*ptrU[32];
                        tmp2[3] = pJ1i[1]*ptrU[3]+pJ1i[7]*ptrU[9]+pJ1i[13]*ptrU[15]+pJ1i[19]*ptrU[21]+pJ1i[25]*ptrU[27]+pJ1i[31]*ptrU[33];
                        tmp2[4] = pJ1i[1]*ptrU[4]+pJ1i[7]*ptrU[10]+pJ1i[13]*ptrU[16]+pJ1i[19]*ptrU[22]+pJ1i[25]*ptrU[28]+pJ1i[31]*ptrU[34];
                        tmp2[5] = pJ1i[1]*ptrU[5]+pJ1i[7]*ptrU[11]+pJ1i[13]*ptrU[17]+pJ1i[19]*ptrU[23]+pJ1i[25]*ptrU[29]+pJ1i[31]*ptrU[35];

                        tmp3[0] = pJ1i[2]*ptrU[0]+pJ1i[8]*ptrU[6]+pJ1i[14]*ptrU[12]+pJ1i[20]*ptrU[18]+pJ1i[26]*ptrU[24]+pJ1i[32]*ptrU[30];
                        tmp3[1] = pJ1i[2]*ptrU[1]+pJ1i[8]*ptrU[7]+pJ1i[14]*ptrU[13]+pJ1i[20]*ptrU[19]+pJ1i[26]*ptrU[25]+pJ1i[32]*ptrU[31];
                        tmp3[2] = pJ1i[2]*ptrU[2]+pJ1i[8]*ptrU[8]+pJ1i[14]*ptrU[14]+pJ1i[20]*ptrU[20]+pJ1i[26]*ptrU[26]+pJ1i[32]*ptrU[32];
                        tmp3[3] = pJ1i[2]*ptrU[3]+pJ1i[8]*ptrU[9]+pJ1i[14]*ptrU[15]+pJ1i[20]*ptrU[21]+pJ1i[26]*ptrU[27]+pJ1i[32]*ptrU[33];
                        tmp3[4] = pJ1i[2]*ptrU[4]+pJ1i[8]*ptrU[10]+pJ1i[14]*ptrU[16]+pJ1i[20]*ptrU[22]+pJ1i[26]*ptrU[28]+pJ1i[32]*ptrU[34];
                        tmp3[5] = pJ1i[2]*ptrU[5]+pJ1i[8]*ptrU[11]+pJ1i[14]*ptrU[17]+pJ1i[20]*ptrU[23]+pJ1i[26]*ptrU[29]+pJ1i[32]*ptrU[35];

                        tmp4[0] = pJ1i[3]*ptrU[0]+pJ1i[9]*ptrU[6]+pJ1i[15]*ptrU[12]+pJ1i[21]*ptrU[18]+pJ1i[27]*ptrU[24]+pJ1i[33]*ptrU[30];
                        tmp4[1] = pJ1i[3]*ptrU[1]+pJ1i[9]*ptrU[7]+pJ1i[15]*ptrU[13]+pJ1i[21]*ptrU[19]+pJ1i[27]*ptrU[25]+pJ1i[33]*ptrU[31];
                        tmp4[2] = pJ1i[3]*ptrU[2]+pJ1i[9]*ptrU[8]+pJ1i[15]*ptrU[14]+pJ1i[21]*ptrU[20]+pJ1i[27]*ptrU[26]+pJ1i[33]*ptrU[32];
                        tmp4[3] = pJ1i[3]*ptrU[3]+pJ1i[9]*ptrU[9]+pJ1i[15]*ptrU[15]+pJ1i[21]*ptrU[21]+pJ1i[27]*ptrU[27]+pJ1i[33]*ptrU[33];
                        tmp4[4] = pJ1i[3]*ptrU[4]+pJ1i[9]*ptrU[10]+pJ1i[15]*ptrU[16]+pJ1i[21]*ptrU[22]+pJ1i[27]*ptrU[28]+pJ1i[33]*ptrU[34];
                        tmp4[5] = pJ1i[3]*ptrU[5]+pJ1i[9]*ptrU[11]+pJ1i[15]*ptrU[17]+pJ1i[21]*ptrU[23]+pJ1i[27]*ptrU[29]+pJ1i[33]*ptrU[35];

                        tmp5[0] = pJ1i[4]*ptrU[0]+pJ1i[10]*ptrU[6]+pJ1i[16]*ptrU[12]+pJ1i[22]*ptrU[18]+pJ1i[28]*ptrU[24]+pJ1i[34]*ptrU[30];
                        tmp5[1] = pJ1i[4]*ptrU[1]+pJ1i[10]*ptrU[7]+pJ1i[16]*ptrU[13]+pJ1i[22]*ptrU[19]+pJ1i[28]*ptrU[25]+pJ1i[34]*ptrU[31];
                        tmp5[2] = pJ1i[4]*ptrU[2]+pJ1i[10]*ptrU[8]+pJ1i[16]*ptrU[14]+pJ1i[22]*ptrU[20]+pJ1i[28]*ptrU[26]+pJ1i[34]*ptrU[32];
                        tmp5[3] = pJ1i[4]*ptrU[3]+pJ1i[10]*ptrU[9]+pJ1i[16]*ptrU[15]+pJ1i[22]*ptrU[21]+pJ1i[28]*ptrU[27]+pJ1i[34]*ptrU[33];
                        tmp5[4] = pJ1i[4]*ptrU[4]+pJ1i[10]*ptrU[10]+pJ1i[16]*ptrU[16]+pJ1i[22]*ptrU[22]+pJ1i[28]*ptrU[28]+pJ1i[34]*ptrU[34];
                        tmp5[5] = pJ1i[4]*ptrU[5]+pJ1i[10]*ptrU[11]+pJ1i[16]*ptrU[17]+pJ1i[22]*ptrU[23]+pJ1i[28]*ptrU[29]+pJ1i[34]*ptrU[35];

                        tmp6[0] = pJ1i[5]*ptrU[0]+pJ1i[11]*ptrU[6]+pJ1i[17]*ptrU[12]+pJ1i[23]*ptrU[18]+pJ1i[29]*ptrU[24]+pJ1i[35]*ptrU[30];
                        tmp6[1] = pJ1i[5]*ptrU[1]+pJ1i[11]*ptrU[7]+pJ1i[17]*ptrU[13]+pJ1i[23]*ptrU[19]+pJ1i[29]*ptrU[25]+pJ1i[35]*ptrU[31];
                        tmp6[2] = pJ1i[5]*ptrU[2]+pJ1i[11]*ptrU[8]+pJ1i[17]*ptrU[14]+pJ1i[23]*ptrU[20]+pJ1i[29]*ptrU[26]+pJ1i[35]*ptrU[32];
                        tmp6[3] = pJ1i[5]*ptrU[3]+pJ1i[11]*ptrU[9]+pJ1i[17]*ptrU[15]+pJ1i[23]*ptrU[21]+pJ1i[29]*ptrU[27]+pJ1i[35]*ptrU[33];
                        tmp6[4] = pJ1i[5]*ptrU[4]+pJ1i[11]*ptrU[10]+pJ1i[17]*ptrU[16]+pJ1i[23]*ptrU[22]+pJ1i[29]*ptrU[28]+pJ1i[35]*ptrU[34];
                        tmp6[5] = pJ1i[5]*ptrU[5]+pJ1i[11]*ptrU[11]+pJ1i[17]*ptrU[17]+pJ1i[23]*ptrU[23]+pJ1i[29]*ptrU[29]+pJ1i[35]*ptrU[35];

                        if ( m_Ui[i] == posID )
                        {
                                ptmp = newU+m_Uj[i]*(6*6);
                        }
                        else if ( m_Uj[i] == posID )
                        {
                                ptmp = newU+m_Ui[i]*(6*6);
                        }
                        else
                        {
                                ptmp = newU + n_newU*6*6;
                                //ptrnU += 6*6;
                                m_nUi[n_newU] = m_Ui[i];
                                m_nUj[n_newU] = m_Uj[i];
                                n_newU += 1;
                        }

                        // ttmp = J1(i)^T*I(j,k)*J1(j)
                        ptmp[0] += tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[6]+tmp1[2]*pJ1j[12]+tmp1[3]*pJ1j[18]+tmp1[4]*pJ1j[24]+tmp1[5]*pJ1j[30];
                        ptmp[1] += tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[7]+tmp1[2]*pJ1j[13]+tmp1[3]*pJ1j[19]+tmp1[4]*pJ1j[25]+tmp1[5]*pJ1j[31];
                        ptmp[2] += tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[8]+tmp1[2]*pJ1j[14]+tmp1[3]*pJ1j[20]+tmp1[4]*pJ1j[26]+tmp1[5]*pJ1j[32];
                        ptmp[3] += tmp1[0]*pJ1j[3]+tmp1[1]*pJ1j[9]+tmp1[2]*pJ1j[15]+tmp1[3]*pJ1j[21]+tmp1[4]*pJ1j[27]+tmp1[5]*pJ1j[33];
                        ptmp[4] += tmp1[0]*pJ1j[4]+tmp1[1]*pJ1j[10]+tmp1[2]*pJ1j[16]+tmp1[3]*pJ1j[22]+tmp1[4]*pJ1j[28]+tmp1[5]*pJ1j[34];
                        ptmp[5] += tmp1[0]*pJ1j[5]+tmp1[1]*pJ1j[11]+tmp1[2]*pJ1j[17]+tmp1[3]*pJ1j[23]+tmp1[4]*pJ1j[29]+tmp1[5]*pJ1j[35];

                        ptmp[6] += tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[6]+tmp2[2]*pJ1j[12]+tmp2[3]*pJ1j[18]+tmp2[4]*pJ1j[24]+tmp2[5]*pJ1j[30];
                        ptmp[7] += tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[7]+tmp2[2]*pJ1j[13]+tmp2[3]*pJ1j[19]+tmp2[4]*pJ1j[25]+tmp2[5]*pJ1j[31];
                        ptmp[8] += tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[8]+tmp2[2]*pJ1j[14]+tmp2[3]*pJ1j[20]+tmp2[4]*pJ1j[26]+tmp2[5]*pJ1j[32];
                        ptmp[9] += tmp2[0]*pJ1j[3]+tmp2[1]*pJ1j[9]+tmp2[2]*pJ1j[15]+tmp2[3]*pJ1j[21]+tmp2[4]*pJ1j[27]+tmp2[5]*pJ1j[33];
                        ptmp[10] += tmp2[0]*pJ1j[4]+tmp2[1]*pJ1j[10]+tmp2[2]*pJ1j[16]+tmp2[3]*pJ1j[22]+tmp2[4]*pJ1j[28]+tmp2[5]*pJ1j[34];
                        ptmp[11] += tmp2[0]*pJ1j[5]+tmp2[1]*pJ1j[11]+tmp2[2]*pJ1j[17]+tmp2[3]*pJ1j[23]+tmp2[4]*pJ1j[29]+tmp2[5]*pJ1j[35];

                        ptmp[12] += tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[6]+tmp3[2]*pJ1j[12]+tmp3[3]*pJ1j[18]+tmp3[4]*pJ1j[24]+tmp3[5]*pJ1j[30];
                        ptmp[13] += tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[7]+tmp3[2]*pJ1j[13]+tmp3[3]*pJ1j[19]+tmp3[4]*pJ1j[25]+tmp3[5]*pJ1j[31];
                        ptmp[14] += tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[8]+tmp3[2]*pJ1j[14]+tmp3[3]*pJ1j[20]+tmp3[4]*pJ1j[26]+tmp3[5]*pJ1j[32];
                        ptmp[15] += tmp3[0]*pJ1j[3]+tmp3[1]*pJ1j[9]+tmp3[2]*pJ1j[15]+tmp3[3]*pJ1j[21]+tmp3[4]*pJ1j[27]+tmp3[5]*pJ1j[33];
                        ptmp[16] += tmp3[0]*pJ1j[4]+tmp3[1]*pJ1j[10]+tmp3[2]*pJ1j[16]+tmp3[3]*pJ1j[22]+tmp3[4]*pJ1j[28]+tmp3[5]*pJ1j[34];
                        ptmp[17] += tmp3[0]*pJ1j[5]+tmp3[1]*pJ1j[11]+tmp3[2]*pJ1j[17]+tmp3[3]*pJ1j[23]+tmp3[4]*pJ1j[29]+tmp3[5]*pJ1j[35];

                        ptmp[18] += tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[6]+tmp4[2]*pJ1j[12]+tmp4[3]*pJ1j[18]+tmp4[4]*pJ1j[24]+tmp4[5]*pJ1j[30];
                        ptmp[19] += tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[7]+tmp4[2]*pJ1j[13]+tmp4[3]*pJ1j[19]+tmp4[4]*pJ1j[25]+tmp4[5]*pJ1j[31];
                        ptmp[20] += tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[8]+tmp4[2]*pJ1j[14]+tmp4[3]*pJ1j[20]+tmp4[4]*pJ1j[26]+tmp4[5]*pJ1j[32];
                        ptmp[21] += tmp4[0]*pJ1j[3]+tmp4[1]*pJ1j[9]+tmp4[2]*pJ1j[15]+tmp4[3]*pJ1j[21]+tmp4[4]*pJ1j[27]+tmp4[5]*pJ1j[33];
                        ptmp[22] += tmp4[0]*pJ1j[4]+tmp4[1]*pJ1j[10]+tmp4[2]*pJ1j[16]+tmp4[3]*pJ1j[22]+tmp4[4]*pJ1j[28]+tmp4[5]*pJ1j[34];
                        ptmp[23] += tmp4[0]*pJ1j[5]+tmp4[1]*pJ1j[11]+tmp4[2]*pJ1j[17]+tmp4[3]*pJ1j[23]+tmp4[4]*pJ1j[29]+tmp4[5]*pJ1j[35];

                        ptmp[24] += tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[6]+tmp5[2]*pJ1j[12]+tmp5[3]*pJ1j[18]+tmp5[4]*pJ1j[24]+tmp5[5]*pJ1j[30];
                        ptmp[25] += tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[7]+tmp5[2]*pJ1j[13]+tmp5[3]*pJ1j[19]+tmp5[4]*pJ1j[25]+tmp5[5]*pJ1j[31];
                        ptmp[26] += tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[8]+tmp5[2]*pJ1j[14]+tmp5[3]*pJ1j[20]+tmp5[4]*pJ1j[26]+tmp5[5]*pJ1j[32];
                        ptmp[27] += tmp5[0]*pJ1j[3]+tmp5[1]*pJ1j[9]+tmp5[2]*pJ1j[15]+tmp5[3]*pJ1j[21]+tmp5[4]*pJ1j[27]+tmp5[5]*pJ1j[33];
                        ptmp[28] += tmp5[0]*pJ1j[4]+tmp5[1]*pJ1j[10]+tmp5[2]*pJ1j[16]+tmp5[3]*pJ1j[22]+tmp5[4]*pJ1j[28]+tmp5[5]*pJ1j[34];
                        ptmp[29] += tmp5[0]*pJ1j[5]+tmp5[1]*pJ1j[11]+tmp5[2]*pJ1j[17]+tmp5[3]*pJ1j[23]+tmp5[4]*pJ1j[29]+tmp5[5]*pJ1j[35];

                        ptmp[30] += tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[6]+tmp6[2]*pJ1j[12]+tmp6[3]*pJ1j[18]+tmp6[4]*pJ1j[24]+tmp6[5]*pJ1j[30];
                        ptmp[31] += tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[7]+tmp6[2]*pJ1j[13]+tmp6[3]*pJ1j[19]+tmp6[4]*pJ1j[25]+tmp6[5]*pJ1j[31];
                        ptmp[32] += tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[8]+tmp6[2]*pJ1j[14]+tmp6[3]*pJ1j[20]+tmp6[4]*pJ1j[26]+tmp6[5]*pJ1j[32];
                        ptmp[33] += tmp6[0]*pJ1j[3]+tmp6[1]*pJ1j[9]+tmp6[2]*pJ1j[15]+tmp6[3]*pJ1j[21]+tmp6[4]*pJ1j[27]+tmp6[5]*pJ1j[33];
                        ptmp[34] += tmp6[0]*pJ1j[4]+tmp6[1]*pJ1j[10]+tmp6[2]*pJ1j[16]+tmp6[3]*pJ1j[22]+tmp6[4]*pJ1j[28]+tmp6[5]*pJ1j[34];
                        ptmp[35] += tmp6[0]*pJ1j[5]+tmp6[1]*pJ1j[11]+tmp6[2]*pJ1j[17]+tmp6[3]*pJ1j[23]+tmp6[4]*pJ1j[29]+tmp6[5]*pJ1j[35];

                        //Algorithm Line 6

                        ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12]+tmp1[3]*pJ2j[18]+tmp1[4]*pJ2j[24]+tmp1[5]*pJ2j[30];
                        ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13]+tmp1[3]*pJ2j[19]+tmp1[4]*pJ2j[25]+tmp1[5]*pJ2j[31];
                        ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14]+tmp1[3]*pJ2j[20]+tmp1[4]*pJ2j[26]+tmp1[5]*pJ2j[32];
                        ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15]+tmp1[3]*pJ2j[21]+tmp1[4]*pJ2j[27]+tmp1[5]*pJ2j[33];
                        ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16]+tmp1[3]*pJ2j[22]+tmp1[4]*pJ2j[28]+tmp1[5]*pJ2j[34];
                        ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17]+tmp1[3]*pJ2j[23]+tmp1[4]*pJ2j[29]+tmp1[5]*pJ2j[35];

                        ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12]+tmp2[3]*pJ2j[18]+tmp2[4]*pJ2j[24]+tmp2[5]*pJ2j[30];
                        ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13]+tmp2[3]*pJ2j[19]+tmp2[4]*pJ2j[25]+tmp2[5]*pJ2j[31];
                        ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14]+tmp2[3]*pJ2j[20]+tmp2[4]*pJ2j[26]+tmp2[5]*pJ2j[32];
                        ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15]+tmp2[3]*pJ2j[21]+tmp2[4]*pJ2j[27]+tmp2[5]*pJ2j[33];
                        ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16]+tmp2[3]*pJ2j[22]+tmp2[4]*pJ2j[28]+tmp2[5]*pJ2j[34];
                        ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17]+tmp2[3]*pJ2j[23]+tmp2[4]*pJ2j[29]+tmp2[5]*pJ2j[35];

                        ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12]+tmp3[3]*pJ2j[18]+tmp3[4]*pJ2j[24]+tmp3[5]*pJ2j[30];
                        ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13]+tmp3[3]*pJ2j[19]+tmp3[4]*pJ2j[25]+tmp3[5]*pJ2j[31];
                        ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14]+tmp3[3]*pJ2j[20]+tmp3[4]*pJ2j[26]+tmp3[5]*pJ2j[32];
                        ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15]+tmp3[3]*pJ2j[21]+tmp3[4]*pJ2j[27]+tmp3[5]*pJ2j[33];
                        ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16]+tmp3[3]*pJ2j[22]+tmp3[4]*pJ2j[28]+tmp3[5]*pJ2j[34];
                        ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17]+tmp3[3]*pJ2j[23]+tmp3[4]*pJ2j[29]+tmp3[5]*pJ2j[35];

                        ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12]+tmp4[3]*pJ2j[18]+tmp4[4]*pJ2j[24]+tmp4[5]*pJ2j[30];
                        ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13]+tmp4[3]*pJ2j[19]+tmp4[4]*pJ2j[25]+tmp4[5]*pJ2j[31];
                        ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14]+tmp4[3]*pJ2j[20]+tmp4[4]*pJ2j[26]+tmp4[5]*pJ2j[32];
                        ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15]+tmp4[3]*pJ2j[21]+tmp4[4]*pJ2j[27]+tmp4[5]*pJ2j[33];
                        ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16]+tmp4[3]*pJ2j[22]+tmp4[4]*pJ2j[28]+tmp4[5]*pJ2j[34];
                        ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17]+tmp4[3]*pJ2j[23]+tmp4[4]*pJ2j[29]+tmp4[5]*pJ2j[35];

                        ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12]+tmp5[3]*pJ2j[18]+tmp5[4]*pJ2j[24]+tmp5[5]*pJ2j[30];
                        ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13]+tmp5[3]*pJ2j[19]+tmp5[4]*pJ2j[25]+tmp5[5]*pJ2j[31];
                        ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14]+tmp5[3]*pJ2j[20]+tmp5[4]*pJ2j[26]+tmp5[5]*pJ2j[32];
                        ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15]+tmp5[3]*pJ2j[21]+tmp5[4]*pJ2j[27]+tmp5[5]*pJ2j[33];
                        ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16]+tmp5[3]*pJ2j[22]+tmp5[4]*pJ2j[28]+tmp5[5]*pJ2j[34];
                        ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17]+tmp5[3]*pJ2j[23]+tmp5[4]*pJ2j[29]+tmp5[5]*pJ2j[35];

                        ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12]+tmp6[3]*pJ2j[18]+tmp6[4]*pJ2j[24]+tmp6[5]*pJ2j[30];
                        ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13]+tmp6[3]*pJ2j[19]+tmp6[4]*pJ2j[25]+tmp6[5]*pJ2j[31];
                        ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14]+tmp6[3]*pJ2j[20]+tmp6[4]*pJ2j[26]+tmp6[5]*pJ2j[32];
                        ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15]+tmp6[3]*pJ2j[21]+tmp6[4]*pJ2j[27]+tmp6[5]*pJ2j[33];
                        ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16]+tmp6[3]*pJ2j[22]+tmp6[4]*pJ2j[28]+tmp6[5]*pJ2j[34];
                        ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17]+tmp6[3]*pJ2j[23]+tmp6[4]*pJ2j[29]+tmp6[5]*pJ2j[35];

                        ptmp = newU+m_Ui[i]*(6*6);
                        if ( m_Ui[i] <= posID )
                        {
                                // ttmp = J1(i)^T*I(j,k)*J2(j)
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];                                   
                        }
                        if ( m_Ui[i] >= posID && m_Ui[i] != m_Uj[i])
                        {
                                //if ( m_Ui[i] != m_Uj[i] )
                                //{
                                // ttmp = ttmp^T
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                        }
                }

				GMap_End.nU = n_newU;
			
                //================================================================
                // Compute newW, the number of elements in newW, as well as the ID of newW
                double* ptrW, *ptrnW, *ptrPID, *ptrV, *prtnV; // W, V 和 newW 指针

                //double* pJ1i, *pJ1j, *pJ2i, *pJ2j, *ptmp;

                int sizenW, n_newW;// n_newU is the real number of elements in newU, it should be smaller than sizenU 
                sizenW = (m_nW+n);
				GMap_End.W = (double*)malloc( sizenW*6*3 * sizeof(double) );
                double*  newW = GMap_End.W; 
                memset( newW, 0, sizenW*6*3 * sizeof(double) );
				GMap_End.feature = (int*)malloc( sizenW * sizeof(int) );
                int*  m_nfeature = GMap_End.feature;
				GMap_End.photo = (int*)malloc( sizenW * sizeof(int) );
                int*  m_nphoto = GMap_End.photo;

				GMap_End.V = (double*)malloc(n*3*3*sizeof(double));
				double*  newV = GMap_End.V;
				memset( newV, 0, n*3*3 * sizeof(double) );
                
                n_newW = 0;
                
                ptrW = m_W;
                ptrV = m_V;
                ptrnW = newW;
				prtnV = newV;

                int PreFID;
                int j=0;
                
                for ( int i = 0; i < n; i++ )
                {
                        PreFID = i;
                        ptrPID = newW + n_newW*6*3;
						GMap_End.FBlock[i] = n_newW;
                        //ptrnW = ptrPID;
                        ptrV = m_V + PreFID*3*3;
                        m_nfeature[n_newW] = PreFID;
                        m_nphoto[n_newW] = posID;
                        n_newW += 1;

                        pJ1i = J1 + m*6*6 + PreFID*3*3;
                        pJ1j = pJ1i;
                        pJ2i = J2 + m*6*6 + PreFID*3*6;
                        pJ2j = pJ2i;

					
                                //Algorithm Line 4

                                tmp1[0] = pJ2i[0]*ptrV[0]+pJ2i[6]*ptrV[3]+pJ2i[12]*ptrV[6];
                                tmp1[1] = pJ2i[0]*ptrV[1]+pJ2i[6]*ptrV[4]+pJ2i[12]*ptrV[7];
                                tmp1[2] = pJ2i[0]*ptrV[2]+pJ2i[6]*ptrV[5]+pJ2i[12]*ptrV[8];
                                

                                tmp2[0] = pJ2i[1]*ptrV[0]+pJ2i[7]*ptrV[3]+pJ2i[13]*ptrV[6];
                                tmp2[1] = pJ2i[1]*ptrV[1]+pJ2i[7]*ptrV[4]+pJ2i[13]*ptrV[7];
                                tmp2[2] = pJ2i[1]*ptrV[2]+pJ2i[7]*ptrV[5]+pJ2i[13]*ptrV[8];
                                
                                tmp3[0] = pJ2i[2]*ptrV[0]+pJ2i[8]*ptrV[3]+pJ2i[14]*ptrV[6];
                                tmp3[1] = pJ2i[2]*ptrV[1]+pJ2i[8]*ptrV[4]+pJ2i[14]*ptrV[7];
                                tmp3[2] = pJ2i[2]*ptrV[2]+pJ2i[8]*ptrV[5]+pJ2i[14]*ptrV[8];
                                
                                tmp4[0] = pJ2i[3]*ptrV[0]+pJ2i[9]*ptrV[3]+pJ2i[15]*ptrV[6];
                                tmp4[1] = pJ2i[3]*ptrV[1]+pJ2i[9]*ptrV[4]+pJ2i[15]*ptrV[7];
                                tmp4[2] = pJ2i[3]*ptrV[2]+pJ2i[9]*ptrV[5]+pJ2i[15]*ptrV[8];
                                
                                tmp5[0] = pJ2i[4]*ptrV[0]+pJ2i[10]*ptrV[3]+pJ2i[16]*ptrV[6];
                                tmp5[1] = pJ2i[4]*ptrV[1]+pJ2i[10]*ptrV[4]+pJ2i[16]*ptrV[7];
                                tmp5[2] = pJ2i[4]*ptrV[2]+pJ2i[10]*ptrV[5]+pJ2i[16]*ptrV[8];
                                
                                tmp6[0] = pJ2i[5]*ptrV[0]+pJ2i[11]*ptrV[3]+pJ2i[17]*ptrV[6];
                                tmp6[1] = pJ2i[5]*ptrV[1]+pJ2i[11]*ptrV[4]+pJ2i[17]*ptrV[7];
                                tmp6[2] = pJ2i[5]*ptrV[2]+pJ2i[11]*ptrV[5]+pJ2i[17]*ptrV[8];
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+posID*(6*6);

                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];           

                                //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID;

								ptmp[0] += ttmp1[0];
								ptmp[1] += ttmp1[1];
								ptmp[2] += ttmp1[2];
								ptmp[3] += ttmp2[0];
								ptmp[4] += ttmp2[1];
								ptmp[5] += ttmp2[2];
								ptmp[6] += ttmp3[0];
								ptmp[7] += ttmp3[1];
								ptmp[8] += ttmp3[2];
								ptmp[9] += ttmp4[0];
								ptmp[10] += ttmp4[1];
								ptmp[11] += ttmp4[2];
								ptmp[12] += ttmp5[0];
								ptmp[13] += ttmp5[1];
								ptmp[14] += ttmp5[2];
								ptmp[15] += ttmp6[0];
								ptmp[16] += ttmp6[1];
								ptmp[17] += ttmp6[2];

								//Algorithm Line 3

                                tmp1[0] = pJ1i[0]*ptrV[0]+pJ1i[3]*ptrV[3]+pJ1i[6]*ptrV[6];
                                tmp1[1] = pJ1i[0]*ptrV[1]+pJ1i[3]*ptrV[4]+pJ1i[6]*ptrV[7];
                                tmp1[2] = pJ1i[0]*ptrV[2]+pJ1i[3]*ptrV[5]+pJ1i[6]*ptrV[8];                              

                                tmp2[0] = pJ1i[1]*ptrV[0]+pJ1i[4]*ptrV[3]+pJ1i[7]*ptrV[6];
                                tmp2[1] = pJ1i[1]*ptrV[1]+pJ1i[4]*ptrV[4]+pJ1i[7]*ptrV[7];
                                tmp2[2] = pJ1i[1]*ptrV[2]+pJ1i[4]*ptrV[5]+pJ1i[7]*ptrV[8];
                                
                                tmp3[0] = pJ1i[2]*ptrV[0]+pJ1i[5]*ptrV[3]+pJ1i[8]*ptrV[6];
                                tmp3[1] = pJ1i[2]*ptrV[1]+pJ1i[5]*ptrV[4]+pJ1i[8]*ptrV[7];
                                tmp3[2] = pJ1i[2]*ptrV[2]+pJ1i[5]*ptrV[5]+pJ1i[8]*ptrV[8];

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ptmp = newV + PreFID*3*3;

                                        ptmp[0] = ttmp1[0];
                                        ptmp[1] = ttmp1[1];
                                        ptmp[2] = ttmp1[2];
                                        ptmp[3] = ttmp2[0];
                                        ptmp[4] = ttmp2[1];
                                        ptmp[5] = ttmp2[2];
                                        ptmp[6] = ttmp3[0];
                                        ptmp[7] = ttmp3[1];
                                        ptmp[8] = ttmp3[2];
																	
                        while (j<m_nW && m_feature[j]==PreFID)
                        {
                        ptrW = m_W + j*6*3;

                        pJ1i = J1 + m_photo[j]*6*6;
                        pJ1j = J1 + m*6*6 + m_feature[j]*3*3;
                        pJ2i = J2 + m_photo[j]*6*6;
                        pJ2j = J2 + m*6*6 + m_feature[j]*3*6;
                                //Algorithm Line 4

                                tmp1[0] = pJ2i[0]*ptrW[0]+pJ2i[6]*ptrW[3]+pJ2i[12]*ptrW[6]+pJ2i[18]*ptrW[9]+pJ2i[24]*ptrW[12]+pJ2i[30]*ptrW[15];
                                tmp1[1] = pJ2i[0]*ptrW[1]+pJ2i[6]*ptrW[4]+pJ2i[12]*ptrW[7]+pJ2i[18]*ptrW[10]+pJ2i[24]*ptrW[13]+pJ2i[30]*ptrW[16];
                                tmp1[2] = pJ2i[0]*ptrW[2]+pJ2i[6]*ptrW[5]+pJ2i[12]*ptrW[8]+pJ2i[18]*ptrW[11]+pJ2i[24]*ptrW[14]+pJ2i[30]*ptrW[17];
                                
                                tmp2[0] = pJ2i[1]*ptrW[0]+pJ2i[7]*ptrW[3]+pJ2i[13]*ptrW[6]+pJ2i[19]*ptrW[9]+pJ2i[25]*ptrW[12]+pJ2i[31]*ptrW[15];
                                tmp2[1] = pJ2i[1]*ptrW[1]+pJ2i[7]*ptrW[4]+pJ2i[13]*ptrW[7]+pJ2i[19]*ptrW[10]+pJ2i[25]*ptrW[13]+pJ2i[31]*ptrW[16];
                                tmp2[2] = pJ2i[1]*ptrW[2]+pJ2i[7]*ptrW[5]+pJ2i[13]*ptrW[8]+pJ2i[19]*ptrW[11]+pJ2i[25]*ptrW[14]+pJ2i[31]*ptrW[17];
                                
                                tmp3[0] = pJ2i[2]*ptrW[0]+pJ2i[8]*ptrW[3]+pJ2i[14]*ptrW[6]+pJ2i[20]*ptrW[9]+pJ2i[26]*ptrW[12]+pJ2i[32]*ptrW[15];
                                tmp3[1] = pJ2i[2]*ptrW[1]+pJ2i[8]*ptrW[4]+pJ2i[14]*ptrW[7]+pJ2i[20]*ptrW[10]+pJ2i[26]*ptrW[13]+pJ2i[32]*ptrW[16];
                                tmp3[2] = pJ2i[2]*ptrW[2]+pJ2i[8]*ptrW[5]+pJ2i[14]*ptrW[8]+pJ2i[20]*ptrW[11]+pJ2i[26]*ptrW[14]+pJ2i[32]*ptrW[17];
                                
                                tmp4[0] = pJ2i[3]*ptrW[0]+pJ2i[9]*ptrW[3]+pJ2i[15]*ptrW[6]+pJ2i[21]*ptrW[9]+pJ2i[27]*ptrW[12]+pJ2i[33]*ptrW[15];
                                tmp4[1] = pJ2i[3]*ptrW[1]+pJ2i[9]*ptrW[4]+pJ2i[15]*ptrW[7]+pJ2i[21]*ptrW[10]+pJ2i[27]*ptrW[13]+pJ2i[33]*ptrW[16];
                                tmp4[2] = pJ2i[3]*ptrW[2]+pJ2i[9]*ptrW[5]+pJ2i[15]*ptrW[8]+pJ2i[21]*ptrW[11]+pJ2i[27]*ptrW[14]+pJ2i[33]*ptrW[17];
                                
                                tmp5[0] = pJ2i[4]*ptrW[0]+pJ2i[10]*ptrW[3]+pJ2i[16]*ptrW[6]+pJ2i[22]*ptrW[9]+pJ2i[28]*ptrW[12]+pJ2i[34]*ptrW[15];
                                tmp5[1] = pJ2i[4]*ptrW[1]+pJ2i[10]*ptrW[4]+pJ2i[16]*ptrW[7]+pJ2i[22]*ptrW[10]+pJ2i[28]*ptrW[13]+pJ2i[34]*ptrW[16];
                                tmp5[2] = pJ2i[4]*ptrW[2]+pJ2i[10]*ptrW[5]+pJ2i[16]*ptrW[8]+pJ2i[22]*ptrW[11]+pJ2i[28]*ptrW[14]+pJ2i[34]*ptrW[17];
                                
                                tmp6[0] = pJ2i[5]*ptrW[0]+pJ2i[11]*ptrW[3]+pJ2i[17]*ptrW[6]+pJ2i[23]*ptrW[9]+pJ2i[29]*ptrW[12]+pJ2i[35]*ptrW[15];
                                tmp6[1] = pJ2i[5]*ptrW[1]+pJ2i[11]*ptrW[4]+pJ2i[17]*ptrW[7]+pJ2i[23]*ptrW[10]+pJ2i[29]*ptrW[13]+pJ2i[35]*ptrW[16];
                                tmp6[2] = pJ2i[5]*ptrW[2]+pJ2i[11]*ptrW[5]+pJ2i[17]*ptrW[8]+pJ2i[23]*ptrW[11]+pJ2i[29]*ptrW[14]+pJ2i[35]*ptrW[17];
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+posID*(6*6);

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                   

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];

                                        //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID;

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp2[0];
                                        ptmp[4] += ttmp2[1];
                                        ptmp[5] += ttmp2[2];
                                        ptmp[6] += ttmp3[0];
                                        ptmp[7] += ttmp3[1];
                                        ptmp[8] += ttmp3[2];
                                        ptmp[9] += ttmp4[0];
                                        ptmp[10] += ttmp4[1];
                                        ptmp[11] += ttmp4[2];
                                        ptmp[12] += ttmp5[0];
                                        ptmp[13] += ttmp5[1];
                                        ptmp[14] += ttmp5[2];
                                        ptmp[15] += ttmp6[0];
                                        ptmp[16] += ttmp6[1];
                                        ptmp[17] += ttmp6[2];

                                //Algorithm Line 3

                                tmp1[0] = pJ1i[0]*ptrW[0]+pJ1i[6]*ptrW[3]+pJ1i[12]*ptrW[6]+pJ1i[18]*ptrW[9]+pJ1i[24]*ptrW[12]+pJ1i[30]*ptrW[15];
                                tmp1[1] = pJ1i[0]*ptrW[1]+pJ1i[6]*ptrW[4]+pJ1i[12]*ptrW[7]+pJ1i[18]*ptrW[10]+pJ1i[24]*ptrW[13]+pJ1i[30]*ptrW[16];
                                tmp1[2] = pJ1i[0]*ptrW[2]+pJ1i[6]*ptrW[5]+pJ1i[12]*ptrW[8]+pJ1i[18]*ptrW[11]+pJ1i[24]*ptrW[14]+pJ1i[30]*ptrW[17];
                                
                                tmp2[0] = pJ1i[1]*ptrW[0]+pJ1i[7]*ptrW[3]+pJ1i[13]*ptrW[6]+pJ1i[19]*ptrW[9]+pJ1i[25]*ptrW[12]+pJ1i[31]*ptrW[15];
                                tmp2[1] = pJ1i[1]*ptrW[1]+pJ1i[7]*ptrW[4]+pJ1i[13]*ptrW[7]+pJ1i[19]*ptrW[10]+pJ1i[25]*ptrW[13]+pJ1i[31]*ptrW[16];
                                tmp2[2] = pJ1i[1]*ptrW[2]+pJ1i[7]*ptrW[5]+pJ1i[13]*ptrW[8]+pJ1i[19]*ptrW[11]+pJ1i[25]*ptrW[14]+pJ1i[31]*ptrW[17];
                                
                                tmp3[0] = pJ1i[2]*ptrW[0]+pJ1i[8]*ptrW[3]+pJ1i[14]*ptrW[6]+pJ1i[20]*ptrW[9]+pJ1i[26]*ptrW[12]+pJ1i[32]*ptrW[15];
                                tmp3[1] = pJ1i[2]*ptrW[1]+pJ1i[8]*ptrW[4]+pJ1i[14]*ptrW[7]+pJ1i[20]*ptrW[10]+pJ1i[26]*ptrW[13]+pJ1i[32]*ptrW[16];
                                tmp3[2] = pJ1i[2]*ptrW[2]+pJ1i[8]*ptrW[5]+pJ1i[14]*ptrW[8]+pJ1i[20]*ptrW[11]+pJ1i[26]*ptrW[14]+pJ1i[32]*ptrW[17];
                                
                                tmp4[0] = pJ1i[3]*ptrW[0]+pJ1i[9]*ptrW[3]+pJ1i[15]*ptrW[6]+pJ1i[21]*ptrW[9]+pJ1i[27]*ptrW[12]+pJ1i[33]*ptrW[15];
                                tmp4[1] = pJ1i[3]*ptrW[1]+pJ1i[9]*ptrW[4]+pJ1i[15]*ptrW[7]+pJ1i[21]*ptrW[10]+pJ1i[27]*ptrW[13]+pJ1i[33]*ptrW[16];
                                tmp4[2] = pJ1i[3]*ptrW[2]+pJ1i[9]*ptrW[5]+pJ1i[15]*ptrW[8]+pJ1i[21]*ptrW[11]+pJ1i[27]*ptrW[14]+pJ1i[33]*ptrW[17];
                                
                                tmp5[0] = pJ1i[4]*ptrW[0]+pJ1i[10]*ptrW[3]+pJ1i[16]*ptrW[6]+pJ1i[22]*ptrW[9]+pJ1i[28]*ptrW[12]+pJ1i[34]*ptrW[15];
                                tmp5[1] = pJ1i[4]*ptrW[1]+pJ1i[10]*ptrW[4]+pJ1i[16]*ptrW[7]+pJ1i[22]*ptrW[10]+pJ1i[28]*ptrW[13]+pJ1i[34]*ptrW[16];
                                tmp5[2] = pJ1i[4]*ptrW[2]+pJ1i[10]*ptrW[5]+pJ1i[16]*ptrW[8]+pJ1i[22]*ptrW[11]+pJ1i[28]*ptrW[14]+pJ1i[34]*ptrW[17];
                                
                                tmp6[0] = pJ1i[5]*ptrW[0]+pJ1i[11]*ptrW[3]+pJ1i[17]*ptrW[6]+pJ1i[23]*ptrW[9]+pJ1i[29]*ptrW[12]+pJ1i[35]*ptrW[15];
                                tmp6[1] = pJ1i[5]*ptrW[1]+pJ1i[11]*ptrW[4]+pJ1i[17]*ptrW[7]+pJ1i[23]*ptrW[10]+pJ1i[29]*ptrW[13]+pJ1i[35]*ptrW[16];
                                tmp6[2] = pJ1i[5]*ptrW[2]+pJ1i[11]*ptrW[5]+pJ1i[17]*ptrW[8]+pJ1i[23]*ptrW[11]+pJ1i[29]*ptrW[14]+pJ1i[35]*ptrW[17];
                                
                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                if ( m_photo[j] == posID )
                                {
                                        ptmp = ptrPID;
                                }
                                else
                                {
                                        ptmp = newW + n_newW*6*3;
                                        m_nfeature[n_newW] = m_feature[j];
                                        m_nphoto[n_newW] = m_photo[j];
                                        n_newW += 1;
                                }

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp2[0];
                                        ptmp[4] += ttmp2[1];
                                        ptmp[5] += ttmp2[2];
                                        ptmp[6] += ttmp3[0];
                                        ptmp[7] += ttmp3[1];
                                        ptmp[8] += ttmp3[2];
                                        ptmp[9] += ttmp4[0];
                                        ptmp[10] += ttmp4[1];
                                        ptmp[11] += ttmp4[2];
                                        ptmp[12] += ttmp5[0];
                                        ptmp[13] += ttmp5[1];
                                        ptmp[14] += ttmp5[2];
                                        ptmp[15] += ttmp6[0];
                                        ptmp[16] += ttmp6[1];
                                        ptmp[17] += ttmp6[2];

                                //Algorithm Line 6
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+m_photo[j]*(6*6);
                                if ( m_photo[j] <= posID )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                    
                                }
                                if ( m_photo[j] >= posID )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];
                                }
                                j++;
                             }
               }
			   			   
			    GMap_End.nW = n_newW;               
				//GMap_End.V = (double*)malloc(n*3*3*sizeof(double));
				//memcpy( GMap_End.V, m_V,n*3*3*sizeof(double) );

				free (J1);
				free (J2);				               
	}
}

void CLinearSFMImp::lmj_PF3D_Divide_ConquerStereo( int nLocalMapCount )
{
	double t0, t1, t2, t3;
	t0 = clock();
	int L = 0;

	while( nLocalMapCount>1 )
	{
		int N2 = nLocalMapCount%2;
		nLocalMapCount = int(nLocalMapCount/2.0 + 0.5);

		int NumLM;
		for ( int i = 0; i < nLocalMapCount; i++ )
		{
			if ( i < nLocalMapCount-1 )
				NumLM = 2;
			else
			{
				if ( N2 == 0 )
					NumLM = 2;
				else
					NumLM = 1;
			}

			for ( int j = 0; j < NumLM; j++ )
			{
				printf( "Join Level %d Local Map %d\n", L, 2*i+j+1 );

				if ( j == 0 )
				{
					m_GMapS = m_LMsetS[2*i+j];
				}
				else
				{
					LocalMapInfoStereo GMap_End;

					//t1 = clock();

					lmj_Transform_PF3DStereo( GMap_End, m_LMsetS[2*i+j].Ref );

					//t2 = clock();					
					//printf( "Tramsformation Time Use:  %lf  sec \n", (t2-t1)*0.001 );

					//t1 = clock();

					free( m_GMapS.U );			m_GMapS.U = NULL;
					free( m_GMapS.V );			m_GMapS.V = NULL;
					if(m_GMapS.nW > 0) 
					{
						free( m_GMapS.W );
						free( m_GMapS.photo );
						free( m_GMapS.feature );
						m_GMapS.W = NULL;
						m_GMapS.photo = NULL;
						m_GMapS.feature = NULL;
					}
					free( m_GMapS.Ui );			m_GMapS.Ui= NULL;
					free( m_GMapS.Uj );			m_GMapS.Uj= NULL;	
					free( m_GMapS.stno );		m_GMapS.stno = NULL;
					free( m_GMapS.stVal );		m_GMapS.stVal = NULL;
					free( m_GMapS.FBlock);		m_GMapS.FBlock = NULL;

					//t2 = clock();					
					//printf( "Free Time Use:  %lf  sec \n", (t2-t1)*0.001 );

					lmj_LinearLS_PF3DStereo(GMap_End,m_LMsetS[2*i+j]);
				}
			}

			printf( "Generate Level %d Local Map %d\n\n", L+1, i+1 );

			if ( (i+1)%2 == 0 )
			{
				//t1 = clock();

				if ( m_GMapS.Ref > m_GMapS.FRef )
				{
					LocalMapInfoStereo GMapTmp;
					lmj_Transform_PF3DStereo( GMapTmp, m_GMapS.FRef );

					//free memory
					free( m_GMapS.U );			m_GMapS.U = NULL;
					free( m_GMapS.V );			m_GMapS.V = NULL;
					if(m_GMapS.nW > 0)
					{
						free( m_GMapS.W );
						free( m_GMapS.photo );
						free( m_GMapS.feature );
						m_GMapS.W = NULL;
						m_GMapS.photo = NULL;
						m_GMapS.feature = NULL;
					}
					free( m_GMapS.Ui );			m_GMapS.Ui= NULL;
					free( m_GMapS.Uj );			m_GMapS.Uj= NULL;	
					free( m_GMapS.stno );		m_GMapS.stno = NULL;
					free( m_GMapS.stVal );		m_GMapS.stVal = NULL;
					free( m_GMapS.FBlock);		m_GMapS.FBlock = NULL;					
					
					m_GMapS = GMapTmp;
				}

				//t2 = clock();					
				//printf( "Transformation Time Use:  %lf  sec \n\n", (t2-t1)*0.001 );

			}

			m_LMsetS[i] = m_GMapS;
		}
		L++;
	}

	//t1 = clock();					
				
	if ( m_GMapS.Ref > m_GMapS.FRef )
	{
		LocalMapInfoStereo GMapTmp;
		lmj_Transform_PF3DStereo( GMapTmp, m_GMapS.FRef );

		//free memory
		free( m_GMapS.U );			m_GMapS.U = NULL;
		free( m_GMapS.V );			m_GMapS.V = NULL;
		if(m_GMapS.nW > 0)
		{
			free( m_GMapS.W );
			free( m_GMapS.photo );
			free( m_GMapS.feature );
			m_GMapS.W = NULL;
			m_GMapS.photo = NULL;
			m_GMapS.feature = NULL;
		}
		free( m_GMapS.Ui );			m_GMapS.Ui= NULL;
		free( m_GMapS.Uj );			m_GMapS.Uj= NULL;	
		free( m_GMapS.stno );		m_GMapS.stno = NULL;
		free( m_GMapS.stVal );		m_GMapS.stVal = NULL;
		free( m_GMapS.FBlock);		m_GMapS.FBlock = NULL;

		m_GMapS = GMapTmp;
	}

	//t2 = clock();					
	//printf( "Final Map Transformation Time Use:  %lf  sec \n\n", (t2-t1)*0.001 );

	t1 = clock();

	t2 = (t1-t0)*0.001;

	printf( "Total Used Time:  %lf  sec\n\n", t2 );

	//Save State Vector
	if ( m_szSt != NULL )
		lmj_SaveStateVector( m_szSt, m_GMapS.stVal, m_GMapS.stno, m_GMapS.m*6+m_GMapS.n*3 );

	if ( m_szPose != NULL && m_szFeature!= NULL )
		lmj_SavePoses_3DPF( m_szPose, m_szFeature, m_GMapS.stno, m_GMapS.stVal, m_GMapS.m*6+m_GMapS.n*3 );

	free( m_GMapS.U );			m_GMapS.U = NULL;
	free( m_GMapS.V );			m_GMapS.V = NULL;
	if(m_GMapS.nW > 0)
	{
		free( m_GMapS.W );
		free( m_GMapS.photo );
		free( m_GMapS.feature );
		m_GMapS.W = NULL;
		m_GMapS.photo = NULL;
		m_GMapS.feature = NULL;
	}			
	free( m_GMapS.Ui );			m_GMapS.Ui= NULL;
	free( m_GMapS.Uj );			m_GMapS.Uj= NULL;	
	free( m_GMapS.stno );		m_GMapS.stno = NULL;
	free( m_GMapS.stVal );		m_GMapS.stVal = NULL;
	free( m_GMapS.FBlock);		m_GMapS.FBlock = NULL;	

	system( "pause" );
}


void CLinearSFMImp::lmj_SaveStateVector( char* szSt, double* st,  int* stno,  int n )
{
	FILE *fp = fopen( szSt, "w" );
	if ( fp == NULL )
	{
		printf( "Please Input Path to Save Final State Vector!" );
		return;
	}

	for ( int i = 0; i < n; i++ )
	{
		fprintf( fp, "%d %lf\n", stno[i], st[i] );
	}

	fclose( fp );
}

void CLinearSFMImp::lmj_solveLinearSFMStereo( double* stVal, double* eb, double* ea,  double* U, double*W, double* V, int* Ui, int* Uj, int* photo, int* feature, int m, int n, int nU, int nW )
{
	int i, j, ii, jj, k, l;
	int cnp = 6, pnp = 3, Usz = 36, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[6*3], sum;

	double t0, t1, t2, t3, t4, t5, t6;

	//t0 = clock();
	char* smask = (char*)malloc( m*m*sizeof(char) );
	memset( smask, 0, m*m*sizeof(char) );
		
	int* mapPhoto = (int*)malloc( n*sizeof(int) );
	int nBase = feature[0];
	int nBaseNum = 1;
	int id = 0;	
	for ( int i = 1; i < nW; i++ )
	{
		int nCur = feature[i];
		if ( nCur == nBase )
		{
			nBaseNum++;
		}
		else
		{
			mapPhoto[id] = nBaseNum;
			id++;
			nBaseNum = 1;
			nBase = feature[i];
		}		
	}
	mapPhoto[n-1] = nBaseNum;
		
	id = 0;	
	for ( int i = 0; i < n; i++ )
	{
		int num = mapPhoto[i];
		for ( int ii = 0; ii < num; ii++ )
		{
			int P1 = photo[id+ii];
			for ( int jj = ii; jj < num; jj++ )
			{
				int P2 = photo[id+jj];
				if ( P1<P2 )
					smask[P1*m+P2] = 1;
				else
					smask[P2*m+P1] = 1;
			}

		}
		id += num;
	}
	
	for ( int i = 0; i < nU; i++ )
	{
		int a = Ui[i];
		int b = Uj[i];

		if ( a < b )
			smask[a*m+b] = 1;
		else
			smask[b*m+a] = 1;
	}

	int nuis;
	for(i=nuis=0, jj=m*m; i<jj; ++i)
		nuis+=(smask[i]!=0);

	sba_crsm Sidxij;
	sba_crsm_alloc(&Sidxij, m, m, nuis);
	for(i=k=0; i<m; ++i)
	{
		Sidxij.rowptr[i]=k;
		ii=i*m;
		for(j=0; j<m; ++j)
			if(smask[ii+j])
			{
				Sidxij.val[k]=k;
				Sidxij.colidx[k++]=j;
			}
	}
	

	Sidxij.rowptr[m]=nuis;	
	double* S = (double*)malloc( 6*6*nuis*sizeof(double) );
	double* E = (double*)malloc( 6*m*sizeof(double) );
	memset( S, 0, 6*6*nuis*sizeof(double) );
	memset( E, 0, 6*m*sizeof(double) );
	double* Vold = (double*)malloc( 3*3*n*sizeof(double) );
	memcpy( Vold, V, 3*3*n*sizeof(double));
	pba_inverseV( V, m, n );	
	//Copy U matrix to S matrix 
	for ( int i = 0; i < nU; i++ )
	{
		int a = Ui[i];
		int b = Uj[i];
		pos1 = sba_crsm_elmidx( &Sidxij, a, b );

		ptr2 = S + pos1*36;
		if ( a == b )
		{					
			ptr1 = U + i * Usz;
			for(ii=0; ii<cnp; ++ii, ptr2+=6)
			{
				ptr2[ii]= ptr1[ii*cnp+ii];
				for(jj=ii+1; jj<cnp; ++jj)
					ptr2[jj]= ptr1[ii*cnp+jj];
			}
		}
		else
		{	
			ptr1 = U + i * Usz;
			for(ii=0; ii<cnp; ++ii, ptr2+=6)
				for(jj=0; jj<cnp; ++jj)
					ptr2[jj]= ptr1[ii*cnp+jj];	
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];


	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( int i = 0; i < n; i++ )
	{
		numfea = mapPhoto[i];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = feature[(pos+j)];
			nP1 = photo[(pos+j)];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + (pos+j)*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;

			//WV
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )
						sum+=ptr3[k]*ptr2[jj*pnp+k];
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = 0; k < numfea; k++ )
			{
				nP2 = photo[pos+k];

				//W(V^-1)W^T
				ptr3 = W + (pos+k)*cnp*3;
				//ptrS = S + (nP1*m*36) + nP2*cnp;

				if ( nP1 == nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}
				if( nP1 < nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}				
			}
			//-W^tb
			ptr5 = eb + nF1*ebsz;
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj];
				ptrE[ii]-= sum;
			}					
		}

		pos += numfea;
	}	

	cholmod_common cS;	
	cholmod_start (&cS) ; 
	int *Ap  = (int*)malloc((m+1)*sizeof(int));
	int * Aii = (int*)malloc(nuis*sizeof(int));
	pba_constructAuxCSSLM( Ap, Aii, m, smask  );

	m_sparseE = cholmod_zeros( cnp*m, 1, CHOLMOD_REAL, &cS);
	double* Ex = (double*)m_sparseE->x;
	int nMaxS = (nuis-m)*36+m*21;	//maximum non-zero element in S matrix 
	
	m_sparseS = cholmod_allocate_sparse(m*6,m*6,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_sparseS->p;		//column pointer
	Si = (int*)m_sparseS->i;		//row pointer

	bool init = false;
	pba_constructCSSLM( Si, Sp, Sx, S, Sidxij, init, m, smask ); //set CSS format using S matrix
	for ( ii = 0; ii < cnp*m; ii++  )
		Ex[ii] = E[ii];	

	bool ordering = true;	
	pba_solveCholmodLM( Ap, Aii, false, ordering, cS, nuis, m);

	double* rx = (double*)m_sparseR->x;
	for ( int i = 0; i < 6*m; i++ )
		stVal[i] = rx[i];
	

	pba_solveFeatures( W, V, ea, eb, stVal, stVal+6*m, m, n, mapPhoto, photo );

	memcpy( V, Vold, 3*3*n*sizeof(double) );
	free( Vold );
	free( smask );
	free( S );
	free( E );
 	free( Ap );
 	free( Aii );
	cholmod_free_factor( &m_factorS, &m_cS );
	cholmod_free_sparse( &m_sparseS, &m_cS );
	cholmod_free_dense( &m_sparseR, &m_cS );
	cholmod_free_dense( &m_sparseE, &m_cS );
	cholmod_finish( &cS );
	sba_crsm_free(&Sidxij);
}

bool CLinearSFMImp::pba_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering, cholmod_common m_cS, int m_nS, int m )
{
	int i, j;
	VectorXi scalarPermutation, blockPermutation;

	ordering = true;
	if (!init)
	{
		if (!ordering)
			m_factorS = cholmod_analyze(m_sparseS, &m_cS); // symbolic factorization
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m);
			if (blockPermutation.size() < m) // double space if resizing
				blockPermutation.resize(2*m);

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m;
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);
			if (! amdStatus) 
				return false;


			// blow up the permutation to the scalar matrix
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_sparseS->ncol);
			if (scalarPermutation.size() < (int)m_sparseS->ncol)
				scalarPermutation.resize(2*m_sparseS->ncol);
			size_t scalarIdx = 0;

			for ( i = 0; i < m; ++i)
			{
				const int &pp = blockPermutation(i);
				int base = (pp==0) ? 0 : pp*6;
				int nCols= 6;

				for ( j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;

			}
			assert(scalarIdx == m_sparseS->ncol);

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;
			m_factorS = cholmod_analyze_p(m_sparseS, scalarPermutation.data(), NULL, 0, &m_cS);
		}
	}

	cholmod_factorize(m_sparseS, m_factorS, &m_cS); 
	m_sparseR = cholmod_solve (CHOLMOD_A, m_factorS, m_sparseE, &m_cS) ;


	return true;
}

void CLinearSFMImp::pba_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, sba_crsm& Sidxij, bool init, int m, char* smask)
{
	int ii, jj, jjj, k;
	int pos1;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;

	Sx = (double*)m_sparseS->x;
	if ( !init)
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				*Sp = nZ;
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
					}
				}
			}
		}
	}
}

void CLinearSFMImp::pba_constructAuxCSSLM( int *Ap, int *Aii, int m, char* smask )
{
	int* Cp = Ap;
	int* Ci = Aii;
	int ii, jj;
	int nZ = 0;
	for ( ii = 0; ii < m; ii++ ) 
	{
		*Cp = nZ;
		for( jj=0; jj<=ii; jj++ )
		{
			if (smask[jj*m+ii]==1)
			{
				*Ci++ = jj;
				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
}

void CLinearSFMImp::lmj_LinearLS_PF3DStereo( LocalMapInfoStereo& GMap_End, LocalMapInfoStereo& GMap_Cur )
{

	double t0, t1, t2;
	//t0 = clock();
	
        LocalMapInfoStereo GMap_Joint;   
		    
        int i;
        int m1, m2, m, n1, n2, n;
        m1 = GMap_End.m;
        m2 = GMap_Cur.m;
        n1 = GMap_End.n;
        n2 = GMap_Cur.n;
        int* p2 = GMap_Cur.stno + m2*6;
        int* iter;
        int pos;
        int* comM1 = (int*)malloc( n1*sizeof(int) );
        int* comM2 = (int*)malloc( n2*sizeof(int) );
        int FID, ncom;
        ncom = 0;
        int*  Curfeature2 = (int*)malloc( GMap_Cur.n * sizeof(int) ); 
        memset( Curfeature2, -1, GMap_Cur.n * sizeof(int) ); 

		int*  findTmp = (int*)malloc( GMap_Cur.n*sizeof(int) );
		for ( int i = 0; i < GMap_Cur.n; i++ )
		{
			findTmp[i] = p2[i*3];
		}
	   
        for ( int i = 0; i < n1; i++ )
        {
                int stno1 = GMap_End.stno[6*m1+i*3];
                //iter = find( p2, p2+n2*3, stno1 );
				iter = find( findTmp, findTmp+n2, stno1 );
                //FID = iter - p2;
				FID = iter -findTmp;

                //if ( FID < n2*3 )
				if ( FID < n2 )
                {
                        comM1[ncom] = i;
                        //comM2[ncom] = FID/3;
						comM2[ncom] = FID;
                        //Curfeature2[FID/3] = i;
						Curfeature2[FID] = i;
                        ncom++;
                }
        }

		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Find Common Features Time Use:  %lf  sec \n", t2 );
		//t0 = clock();

        GMap_Joint.n = n1 + n2 - ncom;
        GMap_Joint.m = m1 + m2;
        m = GMap_Joint.m;
        n = GMap_Joint.n;
        GMap_Joint.nU = GMap_End.nU + GMap_Cur.nU;
        GMap_Joint.nW = GMap_End.nW + GMap_Cur.nW;
        GMap_Joint.U = (double*)malloc( GMap_Joint.nU*6*6*sizeof(double) );
        GMap_Joint.W = (double*)malloc( GMap_Joint.nW*6*3*sizeof(double) );
        GMap_Joint.V = (double*)malloc( GMap_Joint.n *3*3*sizeof(double) );
        GMap_Joint.Ui= (int*)malloc( GMap_Joint.nU*sizeof(int) );
        GMap_Joint.Uj= (int*)malloc( GMap_Joint.nU*sizeof(int) );
        GMap_Joint.feature = (int*)malloc( GMap_Joint.nW*sizeof(int) );
        GMap_Joint.photo = (int*)malloc( GMap_Joint.nW*sizeof(int) );
        GMap_Joint.stno = (int*)malloc( (6*m+n*3)*sizeof(int) );
        GMap_Joint.stVal = (double*)malloc( (6*m+n*3)*sizeof(double) );
		GMap_Joint.FBlock = (int*)malloc( GMap_Joint.n * sizeof(int) );
        memset( GMap_Joint.FBlock, -1, GMap_Cur.n * sizeof(int) );

		GMap_Joint.FRef = GMap_End.FRef;  

        memcpy( GMap_Joint.stno, GMap_End.stno, 6*m1*sizeof(int) );
        memcpy( GMap_Joint.stno+6*m1, GMap_Cur.stno, 6*m2*sizeof(int) );
        memcpy( GMap_Joint.stno+m*6, GMap_End.stno+m1*6, n1*3*sizeof(int) );
		       
        int id = 0;
        for ( int i = 0; i < n2; i++ )
        {
                if (Curfeature2[i]==-1)
                {
                        GMap_Joint.stno[6*m+n1*3+id*3] = GMap_Cur.stno[6*m2+i*3];
                        GMap_Joint.stno[6*m+n1*3+id*3+1] = GMap_Cur.stno[6*m2+i*3+1];
                        GMap_Joint.stno[6*m+n1*3+id*3+2] = GMap_Cur.stno[6*m2+i*3+2];

						Curfeature2[i] = n1+id;

                        id++;
                }
        }
				
		double*  eP = (double*)malloc( m*6 * sizeof(double) ); 
		memset( eP, 0, m*6 * sizeof(double) );
		double*  eF = (double*)malloc( n*3 * sizeof(double) ); 
		memset( eF, 0, n*3 * sizeof(double) );


		
		double *ptr1, *ptr2, *ptmp, *tmp1;
		int *ptri, *ptrj;
		ptr1 = GMap_End.U;
		ptr2 = GMap_Joint.U;
		ptri = GMap_Joint.Ui;
		ptrj = GMap_Joint.Uj;
		for (int i=0; i<GMap_End.nU; i++)
		{
			for (int k=0; k < 36; k++)
				ptr2[k] = ptr1[k];

			ptri[0] = GMap_End.Ui[i];
			ptrj[0] = GMap_End.Uj[i];

			ptmp = eP + ptri[0]*6;
			tmp1 = GMap_End.stVal + GMap_End.Uj[i]*6;

			ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2]+ptr1[3]*tmp1[3]+ptr1[4]*tmp1[4]+ptr1[5]*tmp1[5];
			ptmp[1] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[10]*tmp1[4]+ptr1[11]*tmp1[5];
			ptmp[2] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[15]*tmp1[3]+ptr1[16]*tmp1[4]+ptr1[17]*tmp1[5];
			ptmp[3] += ptr1[18]*tmp1[0]+ptr1[19]*tmp1[1]+ptr1[20]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[22]*tmp1[4]+ptr1[23]*tmp1[5];
			ptmp[4] += ptr1[24]*tmp1[0]+ptr1[25]*tmp1[1]+ptr1[26]*tmp1[2]+ptr1[27]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[29]*tmp1[5];
			ptmp[5] += ptr1[30]*tmp1[0]+ptr1[31]*tmp1[1]+ptr1[32]*tmp1[2]+ptr1[33]*tmp1[3]+ptr1[34]*tmp1[4]+ptr1[35]*tmp1[5];

			if ( GMap_End.Ui[i] != GMap_End.Uj[i])
			{
				ptmp = eP + ptrj[0]*6;
				tmp1 = GMap_End.stVal + GMap_End.Ui[i]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[6]*tmp1[1]+ptr1[12]*tmp1[2]+ptr1[18]*tmp1[3]+ptr1[24]*tmp1[4]+ptr1[30]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[13]*tmp1[2]+ptr1[19]*tmp1[3]+ptr1[25]*tmp1[4]+ptr1[31]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[8]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[20]*tmp1[3]+ptr1[26]*tmp1[4]+ptr1[32]*tmp1[5];
				ptmp[3] += ptr1[3]*tmp1[0]+ptr1[9]*tmp1[1]+ptr1[15]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[27]*tmp1[4]+ptr1[33]*tmp1[5];
				ptmp[4] += ptr1[4]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[16]*tmp1[2]+ptr1[22]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[34]*tmp1[5];
				ptmp[5] += ptr1[5]*tmp1[0]+ptr1[11]*tmp1[1]+ptr1[17]*tmp1[2]+ptr1[23]*tmp1[3]+ptr1[29]*tmp1[4]+ptr1[35]*tmp1[5];

			}

			ptr1 += 36;
			ptr2 += 36;
			ptri += 1;
			ptrj += 1;
		}

		ptr1 = GMap_Cur.U;
		for (int i=0; i<GMap_Cur.nU; i++)
		{
			for (int k=0; k < 36; k++)
				ptr2[k] = ptr1[k];

			ptri[0] = GMap_Cur.Ui[i]+m1;
			ptrj[0] = GMap_Cur.Uj[i]+m1;

			
			ptmp = eP + ptri[0]*6;
			tmp1 = GMap_Cur.stVal + GMap_Cur.Uj[i]*6;

			ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2]+ptr1[3]*tmp1[3]+ptr1[4]*tmp1[4]+ptr1[5]*tmp1[5];
			ptmp[1] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[10]*tmp1[4]+ptr1[11]*tmp1[5];
			ptmp[2] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[15]*tmp1[3]+ptr1[16]*tmp1[4]+ptr1[17]*tmp1[5];
			ptmp[3] += ptr1[18]*tmp1[0]+ptr1[19]*tmp1[1]+ptr1[20]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[22]*tmp1[4]+ptr1[23]*tmp1[5];
			ptmp[4] += ptr1[24]*tmp1[0]+ptr1[25]*tmp1[1]+ptr1[26]*tmp1[2]+ptr1[27]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[29]*tmp1[5];
			ptmp[5] += ptr1[30]*tmp1[0]+ptr1[31]*tmp1[1]+ptr1[32]*tmp1[2]+ptr1[33]*tmp1[3]+ptr1[34]*tmp1[4]+ptr1[35]*tmp1[5];

			if ( GMap_Cur.Ui[i] != GMap_Cur.Uj[i])
			{
				ptmp = eP + ptrj[0]*6;
				tmp1 = GMap_Cur.stVal + GMap_Cur.Ui[i]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[6]*tmp1[1]+ptr1[12]*tmp1[2]+ptr1[18]*tmp1[3]+ptr1[24]*tmp1[4]+ptr1[30]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[13]*tmp1[2]+ptr1[19]*tmp1[3]+ptr1[25]*tmp1[4]+ptr1[31]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[8]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[20]*tmp1[3]+ptr1[26]*tmp1[4]+ptr1[32]*tmp1[5];
				ptmp[3] += ptr1[3]*tmp1[0]+ptr1[9]*tmp1[1]+ptr1[15]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[27]*tmp1[4]+ptr1[33]*tmp1[5];
				ptmp[4] += ptr1[4]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[16]*tmp1[2]+ptr1[22]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[34]*tmp1[5];
				ptmp[5] += ptr1[5]*tmp1[0]+ptr1[11]*tmp1[1]+ptr1[17]*tmp1[2]+ptr1[23]*tmp1[3]+ptr1[29]*tmp1[4]+ptr1[35]*tmp1[5];

			}
			

			ptr1 += 36;
			ptr2 += 36;
			ptri += 1;
			ptrj += 1;
		}

		int*  m_nfeature = GMap_Joint.feature;
		int*  m_nphoto = GMap_Joint.photo;
		
		double *ptr3, *ptr4, *ptr5, *ptr6;
		ptr1 = GMap_End.W;
		ptr3 = GMap_Joint.W;
		ptr4 = GMap_End.V;
		ptr5 = GMap_Joint.V;
		int j = 0, l = 0, a = 0, b, cnt;

		for ( int i = 0; i < GMap_End.n; i++ )	
		{
			for (int k=0; k < 9; k++)
				ptr5[k] = ptr4[k];

			ptmp = eF + i*3;
			tmp1 = GMap_End.stVal + 6*GMap_End.m + i*3;

			ptmp[0] += ptr4[0]*tmp1[0]+ptr4[1]*tmp1[1]+ptr4[2]*tmp1[2];
			ptmp[1] += ptr4[3]*tmp1[0]+ptr4[4]*tmp1[1]+ptr4[5]*tmp1[2];
			ptmp[2] += ptr4[6]*tmp1[0]+ptr4[7]*tmp1[1]+ptr4[8]*tmp1[2];

			
			cnt = 0;
			while (j<GMap_End.nW && GMap_End.feature[j]==i)
			{
				for (int k=0; k < 18; k++)
				{ptr3[k] = ptr1[k];}

				m_nfeature[l] = i;
				m_nphoto[l] = GMap_End.photo[j];


				ptmp = eP + m_nphoto[l]*6;
				tmp1 = GMap_End.stVal + 6*GMap_End.m + i*3;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2];
				ptmp[1] += ptr1[3]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[5]*tmp1[2];
				ptmp[2] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2];
				ptmp[3] += ptr1[9]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[11]*tmp1[2];
				ptmp[4] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2];
				ptmp[5] += ptr1[15]*tmp1[0]+ptr1[16]*tmp1[1]+ptr1[17]*tmp1[2];

				
				ptmp = eF + m_nfeature[l]*3;
				tmp1 = GMap_End.stVal + GMap_End.photo[j]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[3]*tmp1[1]+ptr1[6]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[12]*tmp1[4]+ptr1[15]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[7]*tmp1[2]+ptr1[10]*tmp1[3]+ptr1[13]*tmp1[4]+ptr1[16]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[5]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[11]*tmp1[3]+ptr1[14]*tmp1[4]+ptr1[17]*tmp1[5];
				

				ptr1 += 18;
				ptr3 += 18;
				l++;
				j++;
				cnt += 1;
			}

			if ( a<ncom && i == comM1[a] )			
			{
				ptr6 = GMap_Cur.V+comM2[a]*9;
				for (int k=0; k < 9; k++)
					ptr5[k] += ptr6[k];

				ptmp = eF + i*3;
				tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + comM2[a]*3;

				ptmp[0] += ptr6[0]*tmp1[0]+ptr6[1]*tmp1[1]+ptr6[2]*tmp1[2];
				ptmp[1] += ptr6[3]*tmp1[0]+ptr6[4]*tmp1[1]+ptr6[5]*tmp1[2];
				ptmp[2] += ptr6[6]*tmp1[0]+ptr6[7]*tmp1[1]+ptr6[8]*tmp1[2];

				b = GMap_Cur.FBlock[comM2[a]];
				if ( b !=-1 )
				{
					ptr2 = GMap_Cur.W + 18*b;
					while (b<GMap_Cur.nW && GMap_Cur.feature[b]==comM2[a])
					{
						for (int k=0; k < 18; k++)
						{ptr3[k] = ptr2[k];}

						m_nfeature[l] = i;
						m_nphoto[l] = GMap_Cur.photo[b]+GMap_End.m;


						ptmp = eP + m_nphoto[l]*6;
						tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + comM2[a]*3;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[1]*tmp1[1]+ptr2[2]*tmp1[2];
						ptmp[1] += ptr2[3]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[5]*tmp1[2];
						ptmp[2] += ptr2[6]*tmp1[0]+ptr2[7]*tmp1[1]+ptr2[8]*tmp1[2];
						ptmp[3] += ptr2[9]*tmp1[0]+ptr2[10]*tmp1[1]+ptr2[11]*tmp1[2];
						ptmp[4] += ptr2[12]*tmp1[0]+ptr2[13]*tmp1[1]+ptr2[14]*tmp1[2];
						ptmp[5] += ptr2[15]*tmp1[0]+ptr2[16]*tmp1[1]+ptr2[17]*tmp1[2];

						
						ptmp = eF + m_nfeature[l]*3;
						tmp1 = GMap_Cur.stVal + GMap_Cur.photo[b]*6;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[3]*tmp1[1]+ptr2[6]*tmp1[2]+ptr2[9]*tmp1[3]+ptr2[12]*tmp1[4]+ptr2[15]*tmp1[5];
						ptmp[1] += ptr2[1]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[7]*tmp1[2]+ptr2[10]*tmp1[3]+ptr2[13]*tmp1[4]+ptr2[16]*tmp1[5];
						ptmp[2] += ptr2[2]*tmp1[0]+ptr2[5]*tmp1[1]+ptr2[8]*tmp1[2]+ptr2[11]*tmp1[3]+ptr2[14]*tmp1[4]+ptr2[17]*tmp1[5];
						

						ptr2 += 18;
						ptr3 += 18;
						l++;
						b++;
						cnt += 1;
					}
				}
				a += 1;
			}
			
			ptr4 += 9;
			ptr5 += 9;

			
			if (cnt==0)
				GMap_Joint.FBlock[i] = -1;
			else
				GMap_Joint.FBlock[i] = l-cnt;
				
		}

		j = 0;
		ptr1 = GMap_Cur.W;
		ptr4 = GMap_Cur.V;
	
		for ( int i = 0; i < GMap_Cur.n; i++ )	
		{
			if (Curfeature2[i]>=GMap_End.n)
			{
				for (int k=0; k < 9; k++)
					ptr5[k] = ptr4[k];
				ptr5 += 9;

				ptmp = eF + Curfeature2[i]*3;
				tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + i*3;

				ptmp[0] += ptr4[0]*tmp1[0]+ptr4[1]*tmp1[1]+ptr4[2]*tmp1[2];
				ptmp[1] += ptr4[3]*tmp1[0]+ptr4[4]*tmp1[1]+ptr4[5]*tmp1[2];
				ptmp[2] += ptr4[6]*tmp1[0]+ptr4[7]*tmp1[1]+ptr4[8]*tmp1[2];		

				
				cnt = 0;
				while (j<GMap_Cur.nW && GMap_Cur.feature[j]==i)
				{				
					for (int k=0; k < 18; k++)
					{ptr3[k] = ptr1[k];}

					m_nfeature[l] = Curfeature2[i];
					m_nphoto[l] = GMap_Cur.photo[j]+GMap_End.m;

					ptmp = eP + m_nphoto[l]*6;
					tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + i*3;

					ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2];
					ptmp[1] += ptr1[3]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[5]*tmp1[2];
					ptmp[2] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2];
					ptmp[3] += ptr1[9]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[11]*tmp1[2];
					ptmp[4] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2];
					ptmp[5] += ptr1[15]*tmp1[0]+ptr1[16]*tmp1[1]+ptr1[17]*tmp1[2];

					ptmp = eF + m_nfeature[l]*3;
					tmp1 = GMap_Cur.stVal + GMap_Cur.photo[j]*6;

					ptmp[0] += ptr1[0]*tmp1[0]+ptr1[3]*tmp1[1]+ptr1[6]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[12]*tmp1[4]+ptr1[15]*tmp1[5];
					ptmp[1] += ptr1[1]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[7]*tmp1[2]+ptr1[10]*tmp1[3]+ptr1[13]*tmp1[4]+ptr1[16]*tmp1[5];
					ptmp[2] += ptr1[2]*tmp1[0]+ptr1[5]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[11]*tmp1[3]+ptr1[14]*tmp1[4]+ptr1[17]*tmp1[5];


					ptr1 += 18;
					ptr3 += 18;
					l++;
					j++;
					cnt += 1;
				}
				if (cnt==0)
					GMap_Joint.FBlock[Curfeature2[i]] = -1;
				else
					GMap_Joint.FBlock[Curfeature2[i]] = l-cnt;
			}
			else
			{
				while (j<GMap_Cur.nW && GMap_Cur.feature[j]==i)
				{
					ptr1 += 18;
					j++;
				}
			}
			
			ptr4 += 9;	
		}

		if(comM1 != NULL) free( comM1 );
		if(comM2 != NULL) free( comM2 );
		if(findTmp != NULL) free( findTmp );
		if(Curfeature2 != NULL) free( Curfeature2 );		

		if(GMap_End.U  != NULL) free( GMap_End.U );
		if(GMap_End.V != NULL) free( GMap_End.V );
		if(GMap_End.nW > 0 ) free( GMap_End.W );
		if(GMap_End.Ui != NULL) free( GMap_End.Ui );
		if(GMap_End.Uj != NULL) free( GMap_End.Uj );
		if(GMap_End.nW > 0 ) free( GMap_End.photo );
		if(GMap_End.nW > 0 ) free( GMap_End.feature );
		if( GMap_End.FBlock != NULL) free( GMap_End.FBlock );
		if(GMap_End.stno != NULL) free( GMap_End.stno );
		if(GMap_End.stVal != NULL) free( GMap_End.stVal );
		

		if(GMap_Cur.U != NULL) free( GMap_Cur.U );
		if(GMap_Cur.V != NULL) free( GMap_Cur.V );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.W );
		if(GMap_Cur.Ui != NULL) free( GMap_Cur.Ui );
		if(GMap_Cur.Uj != NULL) free( GMap_Cur.Uj );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.photo );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.feature );
		if(GMap_Cur.FBlock != NULL) free( GMap_Cur.FBlock );
		if(GMap_Cur.stno != NULL) free( GMap_Cur.stno );
		if(GMap_Cur.stVal != NULL) free( GMap_Cur.stVal );


		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Join Two Maps Time Use:  %lf  sec \n", t2 );
		//t0 = clock();

		lmj_solveLinearSFMStereo( GMap_Joint.stVal, eF, eP,  GMap_Joint.U, GMap_Joint.W, GMap_Joint.V, GMap_Joint.Ui, GMap_Joint.Uj, GMap_Joint.photo, GMap_Joint.feature, GMap_Joint.m, GMap_Joint.n, GMap_Joint.nU, GMap_Joint.nW );

		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Solve Time Use:  %lf  sec \n", t2 );

		free( eP );
		free( eF );
		GMap_Joint.Ref = GMap_Cur.Ref;
		GMap_Joint.r = 6*m + 3*n;
		m_GMapS = GMap_Joint;

}

void CLinearSFMImp::pba_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb, int m, int n, int* mapCor, int* photo )
{
	int i, j, ii, jj, pos, numfea;
	int nP1, cnp = 6, pnp = 3;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
	double sum, eb2[6];
	pos = 0;
	for ( i = 0; i < n; i++ )
	{
		ptr1 = eb  + i*3;
		ptr2 = IV   + i*3*3;
		ptr5 = dpb + i*3;
		memset( eb2, 0, sizeof(double)*cnp );		
		//numfea = mapCor[i].size();
		numfea = mapCor[i];

		for ( j = 0; j < numfea; j++ )
		{
			//nP1  = mapCor[i].at(j);
			nP1 = photo[pos];
			ptr3 = W + pos*cnp*3;
			ptr4 = dpa + nP1*cnp;
			//Wta
			for(ii=0; ii<pnp; ++ii)
			{
				for(jj=0, sum=0; jj < cnp; ++jj )
					sum += ptr3[jj*3+ii]*ptr4[jj];
				eb2[ii] += sum;
			}
			pos++;
		}

		//V*(eb-Wta)
		for(ii=0; ii<pnp; ++ii)
		{
			for(jj=0, sum = 0; jj < pnp; jj++ )
				sum += ptr2[ii*3+jj]*(ptr1[jj]-eb2[jj]);
			ptr5[ii] = sum;
		}
	}
}

void CLinearSFMImp::pba_inverseV( double* V, int m, int n )
{
	int i;	
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double *ptr1;
	Matrix3d matV, MatInv;

	//compute V inverse matrix using Eigen that has better performance than Lapack
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; // set ptr1 to point to V_i
		matV << ptr1[0],ptr1[1],ptr1[2],ptr1[3],ptr1[4],ptr1[5],ptr1[6],ptr1[7],ptr1[8];		
		MatInv = matV.inverse();
		ptr1[0] = MatInv(0,0);
		ptr1[4] = MatInv(1,1);
		ptr1[8] = MatInv(2,2);
		ptr1[1] = ptr1[3] = MatInv(0,1);
		ptr1[2] = ptr1[6] = MatInv(0,2);
		ptr1[5] = ptr1[7] = MatInv(1,2);
	}	
}

void CLinearSFMImp::lmj_readInformationStereo( LocalMapInfoStereo& GMap_End, char* szPath )
{
	int i, j;
	FILE* fpMap1 = fopen( szPath, "r" );
	int row, tmp2;
	double tmp3;
	fscanf( fpMap1, "%d", &row );
	GMap_End.Ref = row;

	GMap_End.FRef = GMap_End.Ref;

	fscanf( fpMap1, "%d", &row );
	GMap_End.r = row;
	GMap_End.stno = (int*)malloc( row*sizeof(int) );
	GMap_End.stVal = (double*)malloc( row*sizeof(double) );


	int m = 0, n = 0;
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d %lf", &tmp2, &tmp3 );
		GMap_End.stno[i] = tmp2;
		GMap_End.stVal[i] = tmp3;
	}

	fscanf( fpMap1, "%d", &row );
	GMap_End.m = row;
	fscanf( fpMap1, "%d", &row );
	GMap_End.n = row;

	fscanf( fpMap1, "%d", &row );
	GMap_End.nU = row;
	GMap_End.Ui = (int*)malloc( GMap_End.nU*sizeof(int) );
	GMap_End.Uj = (int*)malloc( GMap_End.nU*sizeof(int) );
	GMap_End.U  = (double*)malloc( GMap_End.nU*6*6*sizeof(double) );
	for ( int i = 0; i < 6*6*row; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.U[i] = tmp3;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.Ui[i] = tmp2;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.Uj[i] = tmp2;
	}

	fscanf( fpMap1, "%d", &row );
	GMap_End.nW = row;
	GMap_End.feature = (int*)malloc( GMap_End.nW*sizeof(int) );
	GMap_End.photo = (int*)malloc( GMap_End.nW*sizeof(int) );
	GMap_End.W  = (double*)malloc( GMap_End.nW*6*3*sizeof(double) );
	for ( int i = 0; i < 6*3*row; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.W[i] = tmp3;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.photo[i] = tmp2;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.feature[i] = tmp2;
	}

	GMap_End.V = (double*)malloc( GMap_End.n*3*3*sizeof(double) );
	for ( int i = 0; i < 3*3*GMap_End.n; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.V[i] = tmp3;
	}

	GMap_End.FBlock = (int*)malloc( GMap_End.n*sizeof(int) );
	for ( int i = 0; i < GMap_End.n; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.FBlock[i] = tmp2;
	}

	fclose( fpMap1 );

}

// Monocular Linear SFM

void CLinearSFMImp::runMono( char* szState, char* szPose, char* szFeature, char* szData, int nMapCount )
{
	m_szPath = szData;
	m_szSt   = szState;
	m_szPose = szPose;
	m_szFeature = szFeature;
	m_LMset = new LocalMapInfo[nMapCount];

	lmj_loadLocalMapsMono( nMapCount );
	
	lmj_PF3D_Divide_ConquerMono( nMapCount );


	free( m_LMset );
	
	system("pause");
}

void CLinearSFMImp::lmj_loadLocalMapsMono( int nMapCount )
{
	FILE* fpSt;
	int i, rows, tmp, n;

	for ( i = 0; i < nMapCount; i++)
	{
		char szSt[200];
		strcpy( szSt, m_szPath );	
		char szPathExt[20];

		sprintf( szPathExt, "/localmap_%d.txt", i+1 );
		strcat( szSt, szPathExt );

		lmj_readInformationMono( m_LMset[i], szSt );
					
	}
}

void CLinearSFMImp::lmj_Transform_PF3DMono( LocalMapInfo& GMap_End, int Ref, int ScaP, int Fix)
{

	if (m_GMap.Ref == Ref && m_GMap.ScaP == ScaP)
    {
         GMap_End = m_GMap;
     }
    else
    {	
		int pos, pos1, pos2, pos3, pos4, i, N, mFix; 
		double t[3], t2[3], ts[3], Scale,  Alpha, Beta, Gamma, R[9], R2[9], R3[9], dRA[9], dRB[9], dRG[9];
		double dA[3], dB[3], dG[3], dRA2[9], dRB2[9], dRG2[9];
		double Alpha2, Beta2, Gamma2, Ri[9]; 
		double tmp1[6], tmp2[6], tmp3[6], tmp4[6], tmp5[6], tmp6[6], ttmp1[6], ttmp2[6], ttmp3[6], ttmp4[6], ttmp5[6], ttmp6[6];
		double ddA2[3], ddB2[3], ddG2[3], ddA[3], ddB[3], ddG[3];
		double* ptr1, *ptr2;
		int* iter;
		int* stno;
		double dSdt[9], dSdA[3], dSdB[3], dSdG[3], dSdtt[9], Sign;
		double Scale2, t222[3], t22[3], dt2dt22[9], dt2dt[9], dt2dtt[9];

		N = 6*m_GMap.m + 3*m_GMap.n;
		GMap_End.r = m_GMap.r;
        
		int m = m_GMap.m;
		int n  = m_GMap.n;
		double* m_U, *m_W, *m_V;
		int* m_Ui, *m_Uj, *m_photo, *m_feature;
		int m_nU, m_nW, m_nFea;	

		m_U = m_GMap.U;
		m_V = m_GMap.V;
		m_W = m_GMap.W;
		m_Ui= m_GMap.Ui;
		m_Uj= m_GMap.Uj;
		m_photo = m_GMap.photo;
		m_feature = m_GMap.feature;
		m_nU = m_GMap.nU;
		m_nW = m_GMap.nW;
		m_nFea = m_GMap.n;
        ptr1 = m_GMap.stVal;
        stno = m_GMap.stno;

		iter = find( stno, stno+N, -Ref );
        pos1 = iter - stno;
		iter = find( stno, stno+N, -ScaP );
        pos2 = iter - stno;

		t[0] = ptr1[pos1];
        t[1] = ptr1[pos1+1];
        t[2] = ptr1[pos1+2];

        Alpha  = ptr1[pos1+3];
        Beta   = ptr1[pos1+4];
        Gamma  = ptr1[pos1+5];

        lmj_RMatrixYPR22( R, Alpha, Beta, Gamma );

		t2[0] = ptr1[pos2];
        t2[1] = ptr1[pos2+1];
        t2[2] = ptr1[pos2+2];

		ts[0] = ( R[0]*(t2[0]-t[0])+R[1]*(t2[1]-t[1])+R[2]*(t2[2]-t[2]) );
		ts[1] = ( R[3]*(t2[0]-t[0])+R[4]*(t2[1]-t[1])+R[5]*(t2[2]-t[2]) );
		ts[2] = ( R[6]*(t2[0]-t[0])+R[7]*(t2[1]-t[1])+R[8]*(t2[2]-t[2]) );

		Scale = abs(ts[Fix]);

		if ( ts[Fix]>=0 )
		{	GMap_End.Sign = 1;}
		else
		{	GMap_End.Sign = -1;}

		GMap_End.m = m_GMap.m;
		GMap_End.n = m_GMap.n;
		GMap_End.stno = (int*)malloc( (6*m+3*n)*sizeof(int) );
		GMap_End.stVal = (double*)malloc( (6*m+3*n)*sizeof(double) );

		GMap_End.FBlock = (int*)malloc(n*sizeof(int));
		memset( GMap_End.FBlock, -1, n * sizeof(int) ); 
		GMap_End.Ref = Ref;
		GMap_End.ScaP = ScaP;
		GMap_End.Fix = Fix;

		GMap_End.FRef = m_GMap.FRef;
		GMap_End.FScaP = m_GMap.FScaP;
		GMap_End.FFix = m_GMap.FFix;

		ptr2 = GMap_End.stVal;
        memcpy( GMap_End.stno, stno, (6*m+3*n)*sizeof(int) );

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================

		for ( i=0; i<N; i++ )
        {
			if ( stno[i] <=0 )
            {
				
				ptr2[i]   = ( R[0]*(ptr1[i]-t[0])+R[1]*(ptr1[i+1]-t[1])+R[2]*(ptr1[i+2]-t[2]) )/Scale;
				ptr2[i+1] = ( R[3]*(ptr1[i]-t[0])+R[4]*(ptr1[i+1]-t[1])+R[5]*(ptr1[i+2]-t[2]) )/Scale;
				ptr2[i+2] = ( R[6]*(ptr1[i]-t[0])+R[7]*(ptr1[i+1]-t[1])+R[8]*(ptr1[i+2]-t[2]) )/Scale;

				lmj_RMatrixYPR22( R2, ptr1[i+3], ptr1[i+4], ptr1[i+5] );

				lmj_TimesRRT( R3, R2, R );
				lmj_InvRotMatrixYPR22( R3, ptr2[i+3], ptr2[i+4], ptr2[i+5] );
				
				if ( i == pos1 )
                {
					ptr2[i]   = 0;
                    ptr2[i+1] = 0;
                    ptr2[i+2] = 0;
					ptr2[i+3] = 0;
                    ptr2[i+4] = 0;
                    ptr2[i+5] = 0;                                 
				}
                if ( i == pos2 )
				{
					ptr2[i+Fix] = GMap_End.Sign;
				}

				i+=5;
			}
            else
			{
				ptr2[i]   = ( R[0]*(ptr1[i]-t[0])+R[1]*(ptr1[i+1]-t[1])+R[2]*(ptr1[i+2]-t[2]) )/Scale;
				ptr2[i+1] = ( R[3]*(ptr1[i]-t[0])+R[4]*(ptr1[i+1]-t[1])+R[5]*(ptr1[i+2]-t[2]) )/Scale;
				ptr2[i+2] = ( R[6]*(ptr1[i]-t[0])+R[7]*(ptr1[i+1]-t[1])+R[8]*(ptr1[i+2]-t[2]) )/Scale;

				i = i+2;
			}
		}

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================
        iter = find( GMap_End.stno, GMap_End.stno+N, -m_GMap.Ref );
        pos3 = iter - GMap_End.stno;
		iter = find( GMap_End.stno, GMap_End.stno+N, -m_GMap.ScaP );
        pos4 = iter - GMap_End.stno;
		
		t[0] = ptr2[pos3];
        t[1] = ptr2[pos3+1];
        t[2] = ptr2[pos3+2];

        Alpha  = ptr2[pos3+3];
        Beta   = ptr2[pos3+4];
        Gamma  = ptr2[pos3+5];

        lmj_Rderivation( Alpha, Beta, Gamma, R, dRA, dRB, dRG );

        lmj_dRiTT( dA, dRA, R );
        lmj_dRiTT( dB, dRB, R );
        lmj_dRiTT( dG, dRG, R );

		t2[0] = ptr2[pos4];
        t2[1] = ptr2[pos4+1];
        t2[2] = ptr2[pos4+2];

		ts[0] = ( R[0]*(t2[0]-t[0])+R[1]*(t2[1]-t[1])+R[2]*(t2[2]-t[2]) );
		ts[1] = ( R[3]*(t2[0]-t[0])+R[4]*(t2[1]-t[1])+R[5]*(t2[2]-t[2]) );
		ts[2] = ( R[6]*(t2[0]-t[0])+R[7]*(t2[1]-t[1])+R[8]*(t2[2]-t[2]) );
		Scale = abs(ts[m_GMap.Fix]);
		Scale2 = Scale*Scale;

		if ( ts[m_GMap.Fix]>=0 )
		{	Sign = 1;}
		else
		{	Sign = -1;}
		mFix = m_GMap.Fix;


		dSdt[0] = -R[0]*Sign; dSdt[1] = -R[1]*Sign; dSdt[2] = -R[2]*Sign; dSdt[3] = -R[3]*Sign; 
		dSdt[4] = -R[4]*Sign; dSdt[5] = -R[5]*Sign; dSdt[6] = -R[6]*Sign; dSdt[7] = -R[7]*Sign;
		dSdt[8] = -R[8]*Sign; 
		
		dSdA[0] = ( dRA[0]*(t2[0]-t[0])+dRA[1]*(t2[1]-t[1])+dRA[2]*(t2[2]-t[2]) )*Sign;
		dSdA[1] = ( dRA[3]*(t2[0]-t[0])+dRA[4]*(t2[1]-t[1])+dRA[5]*(t2[2]-t[2]) )*Sign;
		dSdA[2] = ( dRA[6]*(t2[0]-t[0])+dRA[7]*(t2[1]-t[1])+dRA[8]*(t2[2]-t[2]) )*Sign;

		dSdB[0] = ( dRB[0]*(t2[0]-t[0])+dRB[1]*(t2[1]-t[1])+dRB[2]*(t2[2]-t[2]) )*Sign;
		dSdB[1] = ( dRB[3]*(t2[0]-t[0])+dRB[4]*(t2[1]-t[1])+dRB[5]*(t2[2]-t[2]) )*Sign;
		dSdB[2] = ( dRB[6]*(t2[0]-t[0])+dRB[7]*(t2[1]-t[1])+dRB[8]*(t2[2]-t[2]) )*Sign;

		dSdG[0] = ( dRG[0]*(t2[0]-t[0])+dRG[1]*(t2[1]-t[1])+dRG[2]*(t2[2]-t[2]) )*Sign;
		dSdG[1] = ( dRG[3]*(t2[0]-t[0])+dRG[4]*(t2[1]-t[1])+dRG[5]*(t2[2]-t[2]) )*Sign;
		dSdG[2] = ( dRG[6]*(t2[0]-t[0])+dRG[7]*(t2[1]-t[1])+dRG[8]*(t2[2]-t[2]) )*Sign;

		dSdtt[0] = R[0]*Sign; dSdtt[1] = R[1]*Sign; dSdtt[2] = R[2]*Sign; dSdtt[3] = R[3]*Sign; 
		dSdtt[4] = R[4]*Sign; dSdtt[5] = R[5]*Sign; dSdtt[6] = R[6]*Sign; dSdtt[7] = R[7]*Sign;
		dSdtt[8] = R[8]*Sign;

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //========================================================================================

        int size1 = m*6*6 + n*3*3;
        double *J1 = (double*)malloc( size1*sizeof(double) );
        memset( J1, 0, size1*sizeof(double) );
        int size2 = m*6*6 + n*3*6;
        double *J2 = (double*)malloc( size2*sizeof(double) );
        memset( J2, 0, size2*sizeof(double) );
		double *J3 = (double*)malloc( size2*sizeof(double) );
        memset( J3, 0, size2*sizeof(double) );

        double* pJ1, *pJ2, *pJ3;
        int Id;

        for ( i = 0; i < N; i++ )
		{
			if (GMap_End.stno[i]<=0)
			{
					t2[0] = ptr2[i];
					t2[1] = ptr2[i+1];
					t2[2] = ptr2[i+2];

					Alpha2 = ptr2[i+3];
					Beta2  = ptr2[i+4];
					Gamma2 = ptr2[i+5];

					lmj_Rderivation( Alpha2, Beta2, Gamma2, R2, dRA2, dRB2, dRG2 );

					lmj_TimesRRT( Ri, R2, R );

					double dRidA2[9], dRidB2[9], dRidG2[9];

					lmj_TimesRRT( dRidA2, dRA2, R );
					lmj_TimesRRT( dRidB2, dRB2, R );
					lmj_TimesRRT( dRidG2, dRG2, R );


					double dRidA[9], dRidB[9], dRidG[9];
					lmj_TimesRRT( dRidA, R2, dRA );
					lmj_TimesRRT( dRidB, R2, dRB );
					lmj_TimesRRT( dRidG, R2, dRG );


					lmj_dRi( ddA2, dRidA2, Ri );
					lmj_dRi( ddB2, dRidB2, Ri );
					lmj_dRi( ddG2, dRidG2, Ri );

					lmj_dRi( ddA, dRidA, Ri );
					lmj_dRi( ddB, dRidB, Ri );
					lmj_dRi( ddG, dRidG, Ri );

					t222[0] = t2[0]-t[0];
					t222[1] = t2[1]-t[1];
					t222[2] = t2[2]-t[2];

					t22[0] = R[0]*t222[0] + R[1]*t222[1] + R[2]*t222[2];
					t22[1] = R[3]*t222[0] + R[4]*t222[1] + R[5]*t222[2];
					t22[2] = R[6]*t222[0] + R[7]*t222[1] + R[8]*t222[2];

					//dt2dA = (dRdA*t22*Scale-t2*dSdA(3))/Scale^2;%%
					tmp1[0] = ((dRA[0]*t222[0] + dRA[1]*t222[1] + dRA[2]*t222[2])*Scale - t22[0]*dSdA[mFix])/Scale2;
					tmp1[1] = ((dRA[3]*t222[0] + dRA[4]*t222[1] + dRA[5]*t222[2])*Scale - t22[1]*dSdA[mFix])/Scale2;
					tmp1[2] = ((dRA[6]*t222[0] + dRA[7]*t222[1] + dRA[8]*t222[2])*Scale - t22[2]*dSdA[mFix])/Scale2;
					//dt2dB = (dRdB*t22*Scale-t2*dSdB(3))/Scale^2;%%
					tmp2[0] = ((dRB[0]*t222[0] + dRB[1]*t222[1] + dRB[2]*t222[2])*Scale - t22[0]*dSdB[mFix])/Scale2;
					tmp2[1] = ((dRB[3]*t222[0] + dRB[4]*t222[1] + dRB[5]*t222[2])*Scale - t22[1]*dSdB[mFix])/Scale2;
					tmp2[2] = ((dRB[6]*t222[0] + dRB[7]*t222[1] + dRB[8]*t222[2])*Scale - t22[2]*dSdB[mFix])/Scale2;
					//dt2dG = (dRdG*t22*Scale-t2*dSdG(3))/Scale^2;%%
					tmp3[0] = ((dRG[0]*t222[0] + dRG[1]*t222[1] + dRG[2]*t222[2])*Scale - t22[0]*dSdG[mFix])/Scale2;
					tmp3[1] = ((dRG[3]*t222[0] + dRG[4]*t222[1] + dRG[5]*t222[2])*Scale - t22[1]*dSdG[mFix])/Scale2;
					tmp3[2] = ((dRG[6]*t222[0] + dRG[7]*t222[1] + dRG[8]*t222[2])*Scale - t22[2]*dSdG[mFix])/Scale2;
					//dt2dt22 = R/Scale;%%
					dt2dt22[0] = R[0]/Scale;
					dt2dt22[1] = R[1]/Scale;
					dt2dt22[2] = R[2]/Scale;
					dt2dt22[3] = R[3]/Scale;
					dt2dt22[4] = R[4]/Scale;
					dt2dt22[5] = R[5]/Scale;
					dt2dt22[6] = R[6]/Scale;
					dt2dt22[7] = R[7]/Scale;
					dt2dt22[8] = R[8]/Scale;
					//dt2dt = (-R*Scale-t2*dSdt(3,1:3))/Scale^2;%%
					dt2dt[0] = (-R[0]*Scale-t22[0]*dSdt[3*mFix])/Scale2;
					dt2dt[1] = (-R[1]*Scale-t22[0]*dSdt[3*mFix+1])/Scale2;
					dt2dt[2] = (-R[2]*Scale-t22[0]*dSdt[3*mFix+2])/Scale2;
					dt2dt[3] = (-R[3]*Scale-t22[1]*dSdt[3*mFix])/Scale2;
					dt2dt[4] = (-R[4]*Scale-t22[1]*dSdt[3*mFix+1])/Scale2;
					dt2dt[5] = (-R[5]*Scale-t22[1]*dSdt[3*mFix+2])/Scale2;
					dt2dt[6] = (-R[6]*Scale-t22[2]*dSdt[3*mFix])/Scale2;
					dt2dt[7] = (-R[7]*Scale-t22[2]*dSdt[3*mFix+1])/Scale2;
					dt2dt[8] = (-R[8]*Scale-t22[2]*dSdt[3*mFix+2])/Scale2;
					//dt2dtt = -t2*dSdtt(3,1:3)/Scale^2;%%
					dt2dtt[0] = (-t22[0]*dSdtt[3*mFix])/Scale2;
					dt2dtt[1] = (-t22[0]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[2] = (-t22[0]*dSdtt[3*mFix+2])/Scale2;
					dt2dtt[3] = (-t22[1]*dSdtt[3*mFix])/Scale2;
					dt2dtt[4] = (-t22[1]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[5] = (-t22[1]*dSdtt[3*mFix+2])/Scale2;
					dt2dtt[6] = (-t22[2]*dSdtt[3*mFix])/Scale2;
					dt2dtt[7] = (-t22[2]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[8] = (-t22[2]*dSdtt[3*mFix+2])/Scale2;

					Id = i/6;
                    pJ1 = J1 + Id*6*6;
					pJ2 = J2 + Id*6*6;
					pJ3 = J3 + Id*6*6;

					pJ1[0] += dt2dt22[0];
					pJ1[1] += dt2dt22[1];
					pJ1[2] += dt2dt22[2];
					pJ1[6] += dt2dt22[3];
					pJ1[7] += dt2dt22[4];
					pJ1[8] += dt2dt22[5];
					pJ1[12] += dt2dt22[6];
					pJ1[13] += dt2dt22[7];
					pJ1[14] += dt2dt22[8];
					pJ1[21] += ddA2[0];
                    pJ1[27] += ddA2[1];  
                    pJ1[33] += ddA2[2];
                    pJ1[22] += ddB2[0];  
                    pJ1[28] += ddB2[1];   
                    pJ1[34] += ddB2[2];
                    pJ1[23] += ddG2[0];         
                    pJ1[29] += ddG2[1];  
                    pJ1[35] += ddG2[2]; 

					if (i==pos3)
					{
						pJ1[0] += dt2dt[0];
						pJ1[1] += dt2dt[1];
						pJ1[2] += dt2dt[2];
						pJ1[6] += dt2dt[3];
						pJ1[7] += dt2dt[4];
						pJ1[8] += dt2dt[5];
						pJ1[12] += dt2dt[6];
						pJ1[13] += dt2dt[7];
						pJ1[14] += dt2dt[8];
						pJ1[21] += ddA[0];
						pJ1[27] += ddA[1];  
				        pJ1[33] += ddA[2];
				        pJ1[22] += ddB[0];  
				        pJ1[28] += ddB[1];   
				        pJ1[34] += ddB[2];
				        pJ1[23] += ddG[0];         
				        pJ1[29] += ddG[1];  
						pJ1[35] += ddG[2];

						pJ1[3] += tmp1[0];
						pJ1[4] += tmp2[0];;
						pJ1[5] += tmp3[0];
						pJ1[9] += tmp1[1];
						pJ1[10] += tmp2[1];
						pJ1[11] += tmp3[1];
						pJ1[15] += tmp1[2];
						pJ1[16] += tmp2[2];
						pJ1[17] += tmp3[2];
					}
					else
					{
						pJ2[0] += dt2dt[0];
						pJ2[1] += dt2dt[1];
						pJ2[2] += dt2dt[2];
						pJ2[6] += dt2dt[3];
						pJ2[7] += dt2dt[4];
						pJ2[8] += dt2dt[5];
						pJ2[12] += dt2dt[6];
						pJ2[13] += dt2dt[7];
						pJ2[14] += dt2dt[8];
						pJ2[21] += ddA[0];
						pJ2[27] += ddA[1];  
				        pJ2[33] += ddA[2];
				        pJ2[22] += ddB[0];  
				        pJ2[28] += ddB[1];   
				        pJ2[34] += ddB[2];
				        pJ2[23] += ddG[0];         
				        pJ2[29] += ddG[1];  
						pJ2[35] += ddG[2];

						pJ2[3] += tmp1[0];
						pJ2[4] += tmp2[0];;
						pJ2[5] += tmp3[0];
						pJ2[9] += tmp1[1];
						pJ2[10] += tmp2[1];
						pJ2[11] += tmp3[1];
						pJ2[15] += tmp1[2];
						pJ2[16] += tmp2[2];
						pJ2[17] += tmp3[2];
					}

					if (i==pos4)
					{
						pJ1[0] += dt2dtt[0];
						pJ1[1] += dt2dtt[1];
						pJ1[2] += dt2dtt[2];
						pJ1[6] += dt2dtt[3];
						pJ1[7] += dt2dtt[4];
						pJ1[8] += dt2dtt[5];
						pJ1[12] += dt2dtt[6];
						pJ1[13] += dt2dtt[7];
						pJ1[14] += dt2dtt[8];
					}
					else
					{
						pJ3[0] += dt2dtt[0];
						pJ3[1] += dt2dtt[1];
						pJ3[2] += dt2dtt[2];
						pJ3[6] += dt2dtt[3];
						pJ3[7] += dt2dtt[4];
						pJ3[8] += dt2dtt[5];
						pJ3[12] += dt2dtt[6];
						pJ3[13] += dt2dtt[7];
						pJ3[14] += dt2dtt[8];
					}

					i += 5;
			}
            else
			{
					t2[0] = ptr2[i];
					t2[1] = ptr2[i+1];
					t2[2] = ptr2[i+2];

					t222[0] = t2[0]-t[0];
					t222[1] = t2[1]-t[1];
					t222[2] = t2[2]-t[2];

					t22[0] = R[0]*t222[0] + R[1]*t222[1] + R[2]*t222[2];
					t22[1] = R[3]*t222[0] + R[4]*t222[1] + R[5]*t222[2];
					t22[2] = R[6]*t222[0] + R[7]*t222[1] + R[8]*t222[2];

					//dt2dA = (dRdA*t22*Scale-t2*dSdA(3))/Scale^2;%%
					tmp1[0] = ((dRA[0]*t222[0] + dRA[1]*t222[1] + dRA[2]*t222[2])*Scale - t22[0]*dSdA[mFix])/Scale2;
					tmp1[1] = ((dRA[3]*t222[0] + dRA[4]*t222[1] + dRA[5]*t222[2])*Scale - t22[1]*dSdA[mFix])/Scale2;
					tmp1[2] = ((dRA[6]*t222[0] + dRA[7]*t222[1] + dRA[8]*t222[2])*Scale - t22[2]*dSdA[mFix])/Scale2;
					//dt2dB = (dRdB*t22*Scale-t2*dSdB(3))/Scale^2;%%
					tmp2[0] = ((dRB[0]*t222[0] + dRB[1]*t222[1] + dRB[2]*t222[2])*Scale - t22[0]*dSdB[mFix])/Scale2;
					tmp2[1] = ((dRB[3]*t222[0] + dRB[4]*t222[1] + dRB[5]*t222[2])*Scale - t22[1]*dSdB[mFix])/Scale2;
					tmp2[2] = ((dRB[6]*t222[0] + dRB[7]*t222[1] + dRB[8]*t222[2])*Scale - t22[2]*dSdB[mFix])/Scale2;
					//dt2dG = (dRdG*t22*Scale-t2*dSdG(3))/Scale^2;%%
					tmp3[0] = ((dRG[0]*t222[0] + dRG[1]*t222[1] + dRG[2]*t222[2])*Scale - t22[0]*dSdG[mFix])/Scale2;
					tmp3[1] = ((dRG[3]*t222[0] + dRG[4]*t222[1] + dRG[5]*t222[2])*Scale - t22[1]*dSdG[mFix])/Scale2;
					tmp3[2] = ((dRG[6]*t222[0] + dRG[7]*t222[1] + dRG[8]*t222[2])*Scale - t22[2]*dSdG[mFix])/Scale2;
					//dt2dt22 = R/Scale;%%
					dt2dt22[0] = R[0]/Scale;
					dt2dt22[1] = R[1]/Scale;
					dt2dt22[2] = R[2]/Scale;
					dt2dt22[3] = R[3]/Scale;
					dt2dt22[4] = R[4]/Scale;
					dt2dt22[5] = R[5]/Scale;
					dt2dt22[6] = R[6]/Scale;
					dt2dt22[7] = R[7]/Scale;
					dt2dt22[8] = R[8]/Scale;
					//dt2dt = (-R*Scale-t2*dSdt(3,1:3))/Scale^2;%%
					dt2dt[0] = (-R[0]*Scale-t22[0]*dSdt[3*mFix])/Scale2;
					dt2dt[1] = (-R[1]*Scale-t22[0]*dSdt[3*mFix+1])/Scale2;
					dt2dt[2] = (-R[2]*Scale-t22[0]*dSdt[3*mFix+2])/Scale2;
					dt2dt[3] = (-R[3]*Scale-t22[1]*dSdt[3*mFix])/Scale2;
					dt2dt[4] = (-R[4]*Scale-t22[1]*dSdt[3*mFix+1])/Scale2;
					dt2dt[5] = (-R[5]*Scale-t22[1]*dSdt[3*mFix+2])/Scale2;
					dt2dt[6] = (-R[6]*Scale-t22[2]*dSdt[3*mFix])/Scale2;
					dt2dt[7] = (-R[7]*Scale-t22[2]*dSdt[3*mFix+1])/Scale2;
					dt2dt[8] = (-R[8]*Scale-t22[2]*dSdt[3*mFix+2])/Scale2;
					//dt2dtt = -t2*dSdtt(3,1:3)/Scale^2;%%
					dt2dtt[0] = (-t22[0]*dSdtt[3*mFix])/Scale2;
					dt2dtt[1] = (-t22[0]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[2] = (-t22[0]*dSdtt[3*mFix+2])/Scale2;
					dt2dtt[3] = (-t22[1]*dSdtt[3*mFix])/Scale2;
					dt2dtt[4] = (-t22[1]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[5] = (-t22[1]*dSdtt[3*mFix+2])/Scale2;
					dt2dtt[6] = (-t22[2]*dSdtt[3*mFix])/Scale2;
					dt2dtt[7] = (-t22[2]*dSdtt[3*mFix+1])/Scale2;
					dt2dtt[8] = (-t22[2]*dSdtt[3*mFix+2])/Scale2;

					Id = (i - m*6)/3;
					pJ1 = J1 + m*6*6 + Id*3*3;
					pJ2 = J2 + m*6*6 + Id*3*6;
					pJ3 = J3 + m*6*6 + Id*3*6;

					pJ1[0] += dt2dt22[0];
					pJ1[1] += dt2dt22[1];
					pJ1[2] += dt2dt22[2];
					pJ1[3] += dt2dt22[3];
					pJ1[4] += dt2dt22[4];
					pJ1[5] += dt2dt22[5];
					pJ1[6] += dt2dt22[6];
					pJ1[7] += dt2dt22[7];
					pJ1[8] += dt2dt22[8];
				
						pJ2[0] += dt2dt[0];
						pJ2[1] += dt2dt[1];
						pJ2[2] += dt2dt[2];
						pJ2[6] += dt2dt[3];
						pJ2[7] += dt2dt[4];
						pJ2[8] += dt2dt[5];
						pJ2[12] += dt2dt[6];
						pJ2[13] += dt2dt[7];
						pJ2[14] += dt2dt[8];
						pJ2[3] += tmp1[0];
						pJ2[4] += tmp2[0];;
						pJ2[5] += tmp3[0];
						pJ2[9] += tmp1[1];
						pJ2[10] += tmp2[1];
						pJ2[11] += tmp3[1];
						pJ2[15] += tmp1[2];
						pJ2[16] += tmp2[2];
						pJ2[17] += tmp3[2];

						pJ3[0] += dt2dtt[0];
						pJ3[1] += dt2dtt[1];
						pJ3[2] += dt2dtt[2];
						pJ3[6] += dt2dtt[3];
						pJ3[7] += dt2dtt[4];
						pJ3[8] += dt2dtt[5];
						pJ3[12] += dt2dtt[6];
						pJ3[13] += dt2dtt[7];
						pJ3[14] += dt2dtt[8];

					i += 2;
			}
		}
		
		
		int posID3, posID4;
		posID3 = pos1/6;
		posID4 = pos2/6;

		pJ1 = J1 + 6*6*posID3;
		for(int i=0; i<36; i++)
		{pJ1[i] = 0;}

		pJ1 = J1 + 6*6*posID4;
		for(int i=0; i<6; i++)
		{pJ1[6*i+Fix] = 0;}

		if ( pos3==pos2 )
		{
			for(int i=0; i<6*m+3*n; i++)
			{J2[6*i+Fix] = 0;}
		}

		if ( pos4==pos1 )
		{memset( J3, 0, size2*sizeof(double) );}
		

                double* ptrU, *ptrnU; 

                double* pJ1i, *pJ1j, *pJ2i, *pJ2j, *pJ3i, *pJ3j, *ptmp;

                

                int sizenU, n_newU;// n_newU is the real number of elements in newU, it should be smaller than sizenU 
                sizenU = (m_nU+2*m);
                
				GMap_End.U = (double*)malloc( sizenU*6*6 * sizeof(double) ); 
				double*  newU = GMap_End.U;
                memset( newU, 0, sizenU*6*6 * sizeof(double) );   
				GMap_End.Ui = (int*)malloc( sizenU * sizeof(int) ); 
                int*  m_nUi = GMap_End.Ui;
				GMap_End.Uj = (int*)malloc( sizenU * sizeof(int) );
                int*  m_nUj =GMap_End.Uj;

                n_newU = 2*m;
                int posID, posID2;
					
				posID = pos3/6;
				posID2 = pos4/6;

                ptrU = m_U;
                ptrnU = newU+2*m*(6*6);

                for ( int i = 0; i < m; i++ )
                {
                        if ( i <= posID)
                        {
                                m_nUi[i] = i;
                                m_nUj[i] = posID;
                        }
                        else
                        {
                                m_nUi[i] = posID;
                                m_nUj[i] = i;
                        }
				}

				for ( int i = 0; i < m; i++ )
                {
                        if ( i <= posID)
                        {
                                m_nUi[i+m] = i;
                                m_nUj[i+m] = posID2;
                        }
                        else
                        {
                                m_nUi[i+m] = posID2;
                                m_nUj[i+m] = i;
                        }
				}
                
				for ( int i = 0; i < m_nU; i++ )  
                {
                        ptrU = m_U + i*6*6;

                        pJ1i = J1 + m_Ui[i]*6*6;
                        pJ1j = J1 + m_Uj[i]*6*6;
                        pJ2i = J2 + m_Ui[i]*6*6;
                        pJ2j = J2 + m_Uj[i]*6*6;
						pJ3i = J3 + m_Ui[i]*6*6;
                        pJ3j = J3 + m_Uj[i]*6*6;

                        //Algorithm Line 4
                        tmp1[0] = pJ2i[0]*ptrU[0]+pJ2i[6]*ptrU[6]+pJ2i[12]*ptrU[12]+pJ2i[18]*ptrU[18]+pJ2i[24]*ptrU[24]+pJ2i[30]*ptrU[30];
                        tmp1[1] = pJ2i[0]*ptrU[1]+pJ2i[6]*ptrU[7]+pJ2i[12]*ptrU[13]+pJ2i[18]*ptrU[19]+pJ2i[24]*ptrU[25]+pJ2i[30]*ptrU[31];
                        tmp1[2] = pJ2i[0]*ptrU[2]+pJ2i[6]*ptrU[8]+pJ2i[12]*ptrU[14]+pJ2i[18]*ptrU[20]+pJ2i[24]*ptrU[26]+pJ2i[30]*ptrU[32];
                        tmp1[3] = pJ2i[0]*ptrU[3]+pJ2i[6]*ptrU[9]+pJ2i[12]*ptrU[15]+pJ2i[18]*ptrU[21]+pJ2i[24]*ptrU[27]+pJ2i[30]*ptrU[33];
                        tmp1[4] = pJ2i[0]*ptrU[4]+pJ2i[6]*ptrU[10]+pJ2i[12]*ptrU[16]+pJ2i[18]*ptrU[22]+pJ2i[24]*ptrU[28]+pJ2i[30]*ptrU[34];
                        tmp1[5] = pJ2i[0]*ptrU[5]+pJ2i[6]*ptrU[11]+pJ2i[12]*ptrU[17]+pJ2i[18]*ptrU[23]+pJ2i[24]*ptrU[29]+pJ2i[30]*ptrU[35];

                        tmp2[0] = pJ2i[1]*ptrU[0]+pJ2i[7]*ptrU[6]+pJ2i[13]*ptrU[12]+pJ2i[19]*ptrU[18]+pJ2i[25]*ptrU[24]+pJ2i[31]*ptrU[30];
                        tmp2[1] = pJ2i[1]*ptrU[1]+pJ2i[7]*ptrU[7]+pJ2i[13]*ptrU[13]+pJ2i[19]*ptrU[19]+pJ2i[25]*ptrU[25]+pJ2i[31]*ptrU[31];
                        tmp2[2] = pJ2i[1]*ptrU[2]+pJ2i[7]*ptrU[8]+pJ2i[13]*ptrU[14]+pJ2i[19]*ptrU[20]+pJ2i[25]*ptrU[26]+pJ2i[31]*ptrU[32];
                        tmp2[3] = pJ2i[1]*ptrU[3]+pJ2i[7]*ptrU[9]+pJ2i[13]*ptrU[15]+pJ2i[19]*ptrU[21]+pJ2i[25]*ptrU[27]+pJ2i[31]*ptrU[33];
                        tmp2[4] = pJ2i[1]*ptrU[4]+pJ2i[7]*ptrU[10]+pJ2i[13]*ptrU[16]+pJ2i[19]*ptrU[22]+pJ2i[25]*ptrU[28]+pJ2i[31]*ptrU[34];
                        tmp2[5] = pJ2i[1]*ptrU[5]+pJ2i[7]*ptrU[11]+pJ2i[13]*ptrU[17]+pJ2i[19]*ptrU[23]+pJ2i[25]*ptrU[29]+pJ2i[31]*ptrU[35];

                        tmp3[0] = pJ2i[2]*ptrU[0]+pJ2i[8]*ptrU[6]+pJ2i[14]*ptrU[12]+pJ2i[20]*ptrU[18]+pJ2i[26]*ptrU[24]+pJ2i[32]*ptrU[30];
                        tmp3[1] = pJ2i[2]*ptrU[1]+pJ2i[8]*ptrU[7]+pJ2i[14]*ptrU[13]+pJ2i[20]*ptrU[19]+pJ2i[26]*ptrU[25]+pJ2i[32]*ptrU[31];
                        tmp3[2] = pJ2i[2]*ptrU[2]+pJ2i[8]*ptrU[8]+pJ2i[14]*ptrU[14]+pJ2i[20]*ptrU[20]+pJ2i[26]*ptrU[26]+pJ2i[32]*ptrU[32];
                        tmp3[3] = pJ2i[2]*ptrU[3]+pJ2i[8]*ptrU[9]+pJ2i[14]*ptrU[15]+pJ2i[20]*ptrU[21]+pJ2i[26]*ptrU[27]+pJ2i[32]*ptrU[33];
                        tmp3[4] = pJ2i[2]*ptrU[4]+pJ2i[8]*ptrU[10]+pJ2i[14]*ptrU[16]+pJ2i[20]*ptrU[22]+pJ2i[26]*ptrU[28]+pJ2i[32]*ptrU[34];
                        tmp3[5] = pJ2i[2]*ptrU[5]+pJ2i[8]*ptrU[11]+pJ2i[14]*ptrU[17]+pJ2i[20]*ptrU[23]+pJ2i[26]*ptrU[29]+pJ2i[32]*ptrU[35];

                        tmp4[0] = pJ2i[3]*ptrU[0]+pJ2i[9]*ptrU[6]+pJ2i[15]*ptrU[12]+pJ2i[21]*ptrU[18]+pJ2i[27]*ptrU[24]+pJ2i[33]*ptrU[30];
                        tmp4[1] = pJ2i[3]*ptrU[1]+pJ2i[9]*ptrU[7]+pJ2i[15]*ptrU[13]+pJ2i[21]*ptrU[19]+pJ2i[27]*ptrU[25]+pJ2i[33]*ptrU[31];
                        tmp4[2] = pJ2i[3]*ptrU[2]+pJ2i[9]*ptrU[8]+pJ2i[15]*ptrU[14]+pJ2i[21]*ptrU[20]+pJ2i[27]*ptrU[26]+pJ2i[33]*ptrU[32];
                        tmp4[3] = pJ2i[3]*ptrU[3]+pJ2i[9]*ptrU[9]+pJ2i[15]*ptrU[15]+pJ2i[21]*ptrU[21]+pJ2i[27]*ptrU[27]+pJ2i[33]*ptrU[33];
                        tmp4[4] = pJ2i[3]*ptrU[4]+pJ2i[9]*ptrU[10]+pJ2i[15]*ptrU[16]+pJ2i[21]*ptrU[22]+pJ2i[27]*ptrU[28]+pJ2i[33]*ptrU[34];
                        tmp4[5] = pJ2i[3]*ptrU[5]+pJ2i[9]*ptrU[11]+pJ2i[15]*ptrU[17]+pJ2i[21]*ptrU[23]+pJ2i[27]*ptrU[29]+pJ2i[33]*ptrU[35];

                        tmp5[0] = pJ2i[4]*ptrU[0]+pJ2i[10]*ptrU[6]+pJ2i[16]*ptrU[12]+pJ2i[22]*ptrU[18]+pJ2i[28]*ptrU[24]+pJ2i[34]*ptrU[30];
                        tmp5[1] = pJ2i[4]*ptrU[1]+pJ2i[10]*ptrU[7]+pJ2i[16]*ptrU[13]+pJ2i[22]*ptrU[19]+pJ2i[28]*ptrU[25]+pJ2i[34]*ptrU[31];
                        tmp5[2] = pJ2i[4]*ptrU[2]+pJ2i[10]*ptrU[8]+pJ2i[16]*ptrU[14]+pJ2i[22]*ptrU[20]+pJ2i[28]*ptrU[26]+pJ2i[34]*ptrU[32];
                        tmp5[3] = pJ2i[4]*ptrU[3]+pJ2i[10]*ptrU[9]+pJ2i[16]*ptrU[15]+pJ2i[22]*ptrU[21]+pJ2i[28]*ptrU[27]+pJ2i[34]*ptrU[33];
                        tmp5[4] = pJ2i[4]*ptrU[4]+pJ2i[10]*ptrU[10]+pJ2i[16]*ptrU[16]+pJ2i[22]*ptrU[22]+pJ2i[28]*ptrU[28]+pJ2i[34]*ptrU[34];
                        tmp5[5] = pJ2i[4]*ptrU[5]+pJ2i[10]*ptrU[11]+pJ2i[16]*ptrU[17]+pJ2i[22]*ptrU[23]+pJ2i[28]*ptrU[29]+pJ2i[34]*ptrU[35];

                        tmp6[0] = pJ2i[5]*ptrU[0]+pJ2i[11]*ptrU[6]+pJ2i[17]*ptrU[12]+pJ2i[23]*ptrU[18]+pJ2i[29]*ptrU[24]+pJ2i[35]*ptrU[30];
                        tmp6[1] = pJ2i[5]*ptrU[1]+pJ2i[11]*ptrU[7]+pJ2i[17]*ptrU[13]+pJ2i[23]*ptrU[19]+pJ2i[29]*ptrU[25]+pJ2i[35]*ptrU[31];
                        tmp6[2] = pJ2i[5]*ptrU[2]+pJ2i[11]*ptrU[8]+pJ2i[17]*ptrU[14]+pJ2i[23]*ptrU[20]+pJ2i[29]*ptrU[26]+pJ2i[35]*ptrU[32];
                        tmp6[3] = pJ2i[5]*ptrU[3]+pJ2i[11]*ptrU[9]+pJ2i[17]*ptrU[15]+pJ2i[23]*ptrU[21]+pJ2i[29]*ptrU[27]+pJ2i[35]*ptrU[33];
                        tmp6[4] = pJ2i[5]*ptrU[4]+pJ2i[11]*ptrU[10]+pJ2i[17]*ptrU[16]+pJ2i[23]*ptrU[22]+pJ2i[29]*ptrU[28]+pJ2i[35]*ptrU[34];
                        tmp6[5] = pJ2i[5]*ptrU[5]+pJ2i[11]*ptrU[11]+pJ2i[17]*ptrU[17]+pJ2i[23]*ptrU[23]+pJ2i[29]*ptrU[29]+pJ2i[35]*ptrU[35];

                        ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12]+tmp1[3]*pJ2j[18]+tmp1[4]*pJ2j[24]+tmp1[5]*pJ2j[30];
                        ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13]+tmp1[3]*pJ2j[19]+tmp1[4]*pJ2j[25]+tmp1[5]*pJ2j[31];
                        ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14]+tmp1[3]*pJ2j[20]+tmp1[4]*pJ2j[26]+tmp1[5]*pJ2j[32];
                        ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15]+tmp1[3]*pJ2j[21]+tmp1[4]*pJ2j[27]+tmp1[5]*pJ2j[33];
                        ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16]+tmp1[3]*pJ2j[22]+tmp1[4]*pJ2j[28]+tmp1[5]*pJ2j[34];
                        ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17]+tmp1[3]*pJ2j[23]+tmp1[4]*pJ2j[29]+tmp1[5]*pJ2j[35];

                        ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12]+tmp2[3]*pJ2j[18]+tmp2[4]*pJ2j[24]+tmp2[5]*pJ2j[30];
                        ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13]+tmp2[3]*pJ2j[19]+tmp2[4]*pJ2j[25]+tmp2[5]*pJ2j[31];
                        ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14]+tmp2[3]*pJ2j[20]+tmp2[4]*pJ2j[26]+tmp2[5]*pJ2j[32];
                        ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15]+tmp2[3]*pJ2j[21]+tmp2[4]*pJ2j[27]+tmp2[5]*pJ2j[33];
                        ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16]+tmp2[3]*pJ2j[22]+tmp2[4]*pJ2j[28]+tmp2[5]*pJ2j[34];
                        ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17]+tmp2[3]*pJ2j[23]+tmp2[4]*pJ2j[29]+tmp2[5]*pJ2j[35];

                        ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12]+tmp3[3]*pJ2j[18]+tmp3[4]*pJ2j[24]+tmp3[5]*pJ2j[30];
                        ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13]+tmp3[3]*pJ2j[19]+tmp3[4]*pJ2j[25]+tmp3[5]*pJ2j[31];
                        ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14]+tmp3[3]*pJ2j[20]+tmp3[4]*pJ2j[26]+tmp3[5]*pJ2j[32];
                        ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15]+tmp3[3]*pJ2j[21]+tmp3[4]*pJ2j[27]+tmp3[5]*pJ2j[33];
                        ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16]+tmp3[3]*pJ2j[22]+tmp3[4]*pJ2j[28]+tmp3[5]*pJ2j[34];
                        ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17]+tmp3[3]*pJ2j[23]+tmp3[4]*pJ2j[29]+tmp3[5]*pJ2j[35];

                        ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12]+tmp4[3]*pJ2j[18]+tmp4[4]*pJ2j[24]+tmp4[5]*pJ2j[30];
                        ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13]+tmp4[3]*pJ2j[19]+tmp4[4]*pJ2j[25]+tmp4[5]*pJ2j[31];
                        ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14]+tmp4[3]*pJ2j[20]+tmp4[4]*pJ2j[26]+tmp4[5]*pJ2j[32];
                        ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15]+tmp4[3]*pJ2j[21]+tmp4[4]*pJ2j[27]+tmp4[5]*pJ2j[33];
                        ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16]+tmp4[3]*pJ2j[22]+tmp4[4]*pJ2j[28]+tmp4[5]*pJ2j[34];
                        ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17]+tmp4[3]*pJ2j[23]+tmp4[4]*pJ2j[29]+tmp4[5]*pJ2j[35];

                        ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12]+tmp5[3]*pJ2j[18]+tmp5[4]*pJ2j[24]+tmp5[5]*pJ2j[30];
                        ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13]+tmp5[3]*pJ2j[19]+tmp5[4]*pJ2j[25]+tmp5[5]*pJ2j[31];
                        ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14]+tmp5[3]*pJ2j[20]+tmp5[4]*pJ2j[26]+tmp5[5]*pJ2j[32];
                        ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15]+tmp5[3]*pJ2j[21]+tmp5[4]*pJ2j[27]+tmp5[5]*pJ2j[33];
                        ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16]+tmp5[3]*pJ2j[22]+tmp5[4]*pJ2j[28]+tmp5[5]*pJ2j[34];
                        ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17]+tmp5[3]*pJ2j[23]+tmp5[4]*pJ2j[29]+tmp5[5]*pJ2j[35];

                        ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12]+tmp6[3]*pJ2j[18]+tmp6[4]*pJ2j[24]+tmp6[5]*pJ2j[30];
                        ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13]+tmp6[3]*pJ2j[19]+tmp6[4]*pJ2j[25]+tmp6[5]*pJ2j[31];
                        ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14]+tmp6[3]*pJ2j[20]+tmp6[4]*pJ2j[26]+tmp6[5]*pJ2j[32];
                        ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15]+tmp6[3]*pJ2j[21]+tmp6[4]*pJ2j[27]+tmp6[5]*pJ2j[33];
                        ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16]+tmp6[3]*pJ2j[22]+tmp6[4]*pJ2j[28]+tmp6[5]*pJ2j[34];
                        ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17]+tmp6[3]*pJ2j[23]+tmp6[4]*pJ2j[29]+tmp6[5]*pJ2j[35];

                        ptmp = newU+posID*(6*6);

                        // ttmp = J2(i)^T*I(j,k)*J2(j)
                        ptmp[0] += ttmp1[0];
                        ptmp[1] += ttmp1[1];
                        ptmp[2] += ttmp1[2];
                        ptmp[3] += ttmp1[3];
                        ptmp[4] += ttmp1[4];
                        ptmp[5] += ttmp1[5];
                        ptmp[6] += ttmp2[0];
                        ptmp[7] += ttmp2[1];
                        ptmp[8] += ttmp2[2];
                        ptmp[9] += ttmp2[3];
                        ptmp[10] += ttmp2[4];
                        ptmp[11] += ttmp2[5];
                        ptmp[12] += ttmp3[0];
                        ptmp[13] += ttmp3[1];
                        ptmp[14] += ttmp3[2];
                        ptmp[15] += ttmp3[3];
                        ptmp[16] += ttmp3[4];
                        ptmp[17] += ttmp3[5];
                        ptmp[18] += ttmp4[0];
                        ptmp[19] += ttmp4[1];
                        ptmp[20] += ttmp4[2];
                        ptmp[21] += ttmp4[3];
                        ptmp[22] += ttmp4[4];
                        ptmp[23] += ttmp4[5];
                        ptmp[24] += ttmp5[0];
                        ptmp[25] += ttmp5[1];
                        ptmp[26] += ttmp5[2];
                        ptmp[27] += ttmp5[3];
                        ptmp[28] += ttmp5[4];
                        ptmp[29] += ttmp5[5];
                        ptmp[30] += ttmp6[0];
                        ptmp[31] += ttmp6[1];
                        ptmp[32] += ttmp6[2];
                        ptmp[33] += ttmp6[3];
                        ptmp[34] += ttmp6[4];
                        ptmp[35] += ttmp6[5];                                   


                        if ( m_Ui[i] != m_Uj[i] )
                        {
                                //ttmp = ttmp^T
                                
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];           

                        }


                        //Algorithm Line 5

                        ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[6]+tmp1[2]*pJ1j[12]+tmp1[3]*pJ1j[18]+tmp1[4]*pJ1j[24]+tmp1[5]*pJ1j[30];
                        ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[7]+tmp1[2]*pJ1j[13]+tmp1[3]*pJ1j[19]+tmp1[4]*pJ1j[25]+tmp1[5]*pJ1j[31];
                        ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[8]+tmp1[2]*pJ1j[14]+tmp1[3]*pJ1j[20]+tmp1[4]*pJ1j[26]+tmp1[5]*pJ1j[32];
                        ttmp1[3] = tmp1[0]*pJ1j[3]+tmp1[1]*pJ1j[9]+tmp1[2]*pJ1j[15]+tmp1[3]*pJ1j[21]+tmp1[4]*pJ1j[27]+tmp1[5]*pJ1j[33];
                        ttmp1[4] = tmp1[0]*pJ1j[4]+tmp1[1]*pJ1j[10]+tmp1[2]*pJ1j[16]+tmp1[3]*pJ1j[22]+tmp1[4]*pJ1j[28]+tmp1[5]*pJ1j[34];
                        ttmp1[5] = tmp1[0]*pJ1j[5]+tmp1[1]*pJ1j[11]+tmp1[2]*pJ1j[17]+tmp1[3]*pJ1j[23]+tmp1[4]*pJ1j[29]+tmp1[5]*pJ1j[35];

                        ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[6]+tmp2[2]*pJ1j[12]+tmp2[3]*pJ1j[18]+tmp2[4]*pJ1j[24]+tmp2[5]*pJ1j[30];
                        ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[7]+tmp2[2]*pJ1j[13]+tmp2[3]*pJ1j[19]+tmp2[4]*pJ1j[25]+tmp2[5]*pJ1j[31];
                        ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[8]+tmp2[2]*pJ1j[14]+tmp2[3]*pJ1j[20]+tmp2[4]*pJ1j[26]+tmp2[5]*pJ1j[32];
                        ttmp2[3] = tmp2[0]*pJ1j[3]+tmp2[1]*pJ1j[9]+tmp2[2]*pJ1j[15]+tmp2[3]*pJ1j[21]+tmp2[4]*pJ1j[27]+tmp2[5]*pJ1j[33];
                        ttmp2[4] = tmp2[0]*pJ1j[4]+tmp2[1]*pJ1j[10]+tmp2[2]*pJ1j[16]+tmp2[3]*pJ1j[22]+tmp2[4]*pJ1j[28]+tmp2[5]*pJ1j[34];
                        ttmp2[5] = tmp2[0]*pJ1j[5]+tmp2[1]*pJ1j[11]+tmp2[2]*pJ1j[17]+tmp2[3]*pJ1j[23]+tmp2[4]*pJ1j[29]+tmp2[5]*pJ1j[35];

                        ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[6]+tmp3[2]*pJ1j[12]+tmp3[3]*pJ1j[18]+tmp3[4]*pJ1j[24]+tmp3[5]*pJ1j[30];
                        ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[7]+tmp3[2]*pJ1j[13]+tmp3[3]*pJ1j[19]+tmp3[4]*pJ1j[25]+tmp3[5]*pJ1j[31];
                        ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[8]+tmp3[2]*pJ1j[14]+tmp3[3]*pJ1j[20]+tmp3[4]*pJ1j[26]+tmp3[5]*pJ1j[32];
                        ttmp3[3] = tmp3[0]*pJ1j[3]+tmp3[1]*pJ1j[9]+tmp3[2]*pJ1j[15]+tmp3[3]*pJ1j[21]+tmp3[4]*pJ1j[27]+tmp3[5]*pJ1j[33];
                        ttmp3[4] = tmp3[0]*pJ1j[4]+tmp3[1]*pJ1j[10]+tmp3[2]*pJ1j[16]+tmp3[3]*pJ1j[22]+tmp3[4]*pJ1j[28]+tmp3[5]*pJ1j[34];
                        ttmp3[5] = tmp3[0]*pJ1j[5]+tmp3[1]*pJ1j[11]+tmp3[2]*pJ1j[17]+tmp3[3]*pJ1j[23]+tmp3[4]*pJ1j[29]+tmp3[5]*pJ1j[35];

                        ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[6]+tmp4[2]*pJ1j[12]+tmp4[3]*pJ1j[18]+tmp4[4]*pJ1j[24]+tmp4[5]*pJ1j[30];
                        ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[7]+tmp4[2]*pJ1j[13]+tmp4[3]*pJ1j[19]+tmp4[4]*pJ1j[25]+tmp4[5]*pJ1j[31];
                        ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[8]+tmp4[2]*pJ1j[14]+tmp4[3]*pJ1j[20]+tmp4[4]*pJ1j[26]+tmp4[5]*pJ1j[32];
                        ttmp4[3] = tmp4[0]*pJ1j[3]+tmp4[1]*pJ1j[9]+tmp4[2]*pJ1j[15]+tmp4[3]*pJ1j[21]+tmp4[4]*pJ1j[27]+tmp4[5]*pJ1j[33];
                        ttmp4[4] = tmp4[0]*pJ1j[4]+tmp4[1]*pJ1j[10]+tmp4[2]*pJ1j[16]+tmp4[3]*pJ1j[22]+tmp4[4]*pJ1j[28]+tmp4[5]*pJ1j[34];
                        ttmp4[5] = tmp4[0]*pJ1j[5]+tmp4[1]*pJ1j[11]+tmp4[2]*pJ1j[17]+tmp4[3]*pJ1j[23]+tmp4[4]*pJ1j[29]+tmp4[5]*pJ1j[35];

                        ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[6]+tmp5[2]*pJ1j[12]+tmp5[3]*pJ1j[18]+tmp5[4]*pJ1j[24]+tmp5[5]*pJ1j[30];
                        ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[7]+tmp5[2]*pJ1j[13]+tmp5[3]*pJ1j[19]+tmp5[4]*pJ1j[25]+tmp5[5]*pJ1j[31];
                        ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[8]+tmp5[2]*pJ1j[14]+tmp5[3]*pJ1j[20]+tmp5[4]*pJ1j[26]+tmp5[5]*pJ1j[32];
                        ttmp5[3] = tmp5[0]*pJ1j[3]+tmp5[1]*pJ1j[9]+tmp5[2]*pJ1j[15]+tmp5[3]*pJ1j[21]+tmp5[4]*pJ1j[27]+tmp5[5]*pJ1j[33];
                        ttmp5[4] = tmp5[0]*pJ1j[4]+tmp5[1]*pJ1j[10]+tmp5[2]*pJ1j[16]+tmp5[3]*pJ1j[22]+tmp5[4]*pJ1j[28]+tmp5[5]*pJ1j[34];
                        ttmp5[5] = tmp5[0]*pJ1j[5]+tmp5[1]*pJ1j[11]+tmp5[2]*pJ1j[17]+tmp5[3]*pJ1j[23]+tmp5[4]*pJ1j[29]+tmp5[5]*pJ1j[35];

                        ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[6]+tmp6[2]*pJ1j[12]+tmp6[3]*pJ1j[18]+tmp6[4]*pJ1j[24]+tmp6[5]*pJ1j[30];
                        ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[7]+tmp6[2]*pJ1j[13]+tmp6[3]*pJ1j[19]+tmp6[4]*pJ1j[25]+tmp6[5]*pJ1j[31];
                        ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[8]+tmp6[2]*pJ1j[14]+tmp6[3]*pJ1j[20]+tmp6[4]*pJ1j[26]+tmp6[5]*pJ1j[32];
                        ttmp6[3] = tmp6[0]*pJ1j[3]+tmp6[1]*pJ1j[9]+tmp6[2]*pJ1j[15]+tmp6[3]*pJ1j[21]+tmp6[4]*pJ1j[27]+tmp6[5]*pJ1j[33];
                        ttmp6[4] = tmp6[0]*pJ1j[4]+tmp6[1]*pJ1j[10]+tmp6[2]*pJ1j[16]+tmp6[3]*pJ1j[22]+tmp6[4]*pJ1j[28]+tmp6[5]*pJ1j[34];
                        ttmp6[5] = tmp6[0]*pJ1j[5]+tmp6[1]*pJ1j[11]+tmp6[2]*pJ1j[17]+tmp6[3]*pJ1j[23]+tmp6[4]*pJ1j[29]+tmp6[5]*pJ1j[35];

                        ptmp = newU+m_Uj[i]*(6*6);
                        if ( m_Uj[i] >= posID )
                        {
                                // ttmp = J2(i)^T*I(j,k)*J1(j)
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];
                        }
                        if ( m_Uj[i] <= posID && m_Ui[i] != m_Uj[i] )
                        {
                                //if ( m_Ui[i] != m_Uj[i] )
                                //{
                                // ttmp = ttmp^T
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                                //}


                        }

						///////////////////////////////////////////////////
						//J2i^T*I(i,j)*J3j

                        ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12]+tmp1[3]*pJ3j[18]+tmp1[4]*pJ3j[24]+tmp1[5]*pJ3j[30];
                        ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13]+tmp1[3]*pJ3j[19]+tmp1[4]*pJ3j[25]+tmp1[5]*pJ3j[31];
                        ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14]+tmp1[3]*pJ3j[20]+tmp1[4]*pJ3j[26]+tmp1[5]*pJ3j[32];
                        ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15]+tmp1[3]*pJ3j[21]+tmp1[4]*pJ3j[27]+tmp1[5]*pJ3j[33];
                        ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16]+tmp1[3]*pJ3j[22]+tmp1[4]*pJ3j[28]+tmp1[5]*pJ3j[34];
                        ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17]+tmp1[3]*pJ3j[23]+tmp1[4]*pJ3j[29]+tmp1[5]*pJ3j[35];

                        ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12]+tmp2[3]*pJ3j[18]+tmp2[4]*pJ3j[24]+tmp2[5]*pJ3j[30];
                        ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13]+tmp2[3]*pJ3j[19]+tmp2[4]*pJ3j[25]+tmp2[5]*pJ3j[31];
                        ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14]+tmp2[3]*pJ3j[20]+tmp2[4]*pJ3j[26]+tmp2[5]*pJ3j[32];
                        ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15]+tmp2[3]*pJ3j[21]+tmp2[4]*pJ3j[27]+tmp2[5]*pJ3j[33];
                        ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16]+tmp2[3]*pJ3j[22]+tmp2[4]*pJ3j[28]+tmp2[5]*pJ3j[34];
                        ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17]+tmp2[3]*pJ3j[23]+tmp2[4]*pJ3j[29]+tmp2[5]*pJ3j[35];

                        ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12]+tmp3[3]*pJ3j[18]+tmp3[4]*pJ3j[24]+tmp3[5]*pJ3j[30];
                        ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13]+tmp3[3]*pJ3j[19]+tmp3[4]*pJ3j[25]+tmp3[5]*pJ3j[31];
                        ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14]+tmp3[3]*pJ3j[20]+tmp3[4]*pJ3j[26]+tmp3[5]*pJ3j[32];
                        ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15]+tmp3[3]*pJ3j[21]+tmp3[4]*pJ3j[27]+tmp3[5]*pJ3j[33];
                        ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16]+tmp3[3]*pJ3j[22]+tmp3[4]*pJ3j[28]+tmp3[5]*pJ3j[34];
                        ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17]+tmp3[3]*pJ3j[23]+tmp3[4]*pJ3j[29]+tmp3[5]*pJ3j[35];

                        ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12]+tmp4[3]*pJ3j[18]+tmp4[4]*pJ3j[24]+tmp4[5]*pJ3j[30];
                        ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13]+tmp4[3]*pJ3j[19]+tmp4[4]*pJ3j[25]+tmp4[5]*pJ3j[31];
                        ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14]+tmp4[3]*pJ3j[20]+tmp4[4]*pJ3j[26]+tmp4[5]*pJ3j[32];
                        ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15]+tmp4[3]*pJ3j[21]+tmp4[4]*pJ3j[27]+tmp4[5]*pJ3j[33];
                        ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16]+tmp4[3]*pJ3j[22]+tmp4[4]*pJ3j[28]+tmp4[5]*pJ3j[34];
                        ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17]+tmp4[3]*pJ3j[23]+tmp4[4]*pJ3j[29]+tmp4[5]*pJ3j[35];

                        ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12]+tmp5[3]*pJ3j[18]+tmp5[4]*pJ3j[24]+tmp5[5]*pJ3j[30];
                        ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13]+tmp5[3]*pJ3j[19]+tmp5[4]*pJ3j[25]+tmp5[5]*pJ3j[31];
                        ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14]+tmp5[3]*pJ3j[20]+tmp5[4]*pJ3j[26]+tmp5[5]*pJ3j[32];
                        ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15]+tmp5[3]*pJ3j[21]+tmp5[4]*pJ3j[27]+tmp5[5]*pJ3j[33];
                        ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16]+tmp5[3]*pJ3j[22]+tmp5[4]*pJ3j[28]+tmp5[5]*pJ3j[34];
                        ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17]+tmp5[3]*pJ3j[23]+tmp5[4]*pJ3j[29]+tmp5[5]*pJ3j[35];

                        ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12]+tmp6[3]*pJ3j[18]+tmp6[4]*pJ3j[24]+tmp6[5]*pJ3j[30];
                        ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13]+tmp6[3]*pJ3j[19]+tmp6[4]*pJ3j[25]+tmp6[5]*pJ3j[31];
                        ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14]+tmp6[3]*pJ3j[20]+tmp6[4]*pJ3j[26]+tmp6[5]*pJ3j[32];
                        ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15]+tmp6[3]*pJ3j[21]+tmp6[4]*pJ3j[27]+tmp6[5]*pJ3j[33];
                        ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16]+tmp6[3]*pJ3j[22]+tmp6[4]*pJ3j[28]+tmp6[5]*pJ3j[34];
                        ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17]+tmp6[3]*pJ3j[23]+tmp6[4]*pJ3j[29]+tmp6[5]*pJ3j[35];

						ptmp = newU+posID2*(6*6);
                        if ( posID2 > posID )
                        {                        
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];
                        }
                        if ( posID2 < posID && m_Ui[i] != m_Uj[i] )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                        }
						//////////////////////////////////////////////////

                        //Algorithm Line 3

                        tmp1[0] = pJ1i[0]*ptrU[0]+pJ1i[6]*ptrU[6]+pJ1i[12]*ptrU[12]+pJ1i[18]*ptrU[18]+pJ1i[24]*ptrU[24]+pJ1i[30]*ptrU[30];
                        tmp1[1] = pJ1i[0]*ptrU[1]+pJ1i[6]*ptrU[7]+pJ1i[12]*ptrU[13]+pJ1i[18]*ptrU[19]+pJ1i[24]*ptrU[25]+pJ1i[30]*ptrU[31];
                        tmp1[2] = pJ1i[0]*ptrU[2]+pJ1i[6]*ptrU[8]+pJ1i[12]*ptrU[14]+pJ1i[18]*ptrU[20]+pJ1i[24]*ptrU[26]+pJ1i[30]*ptrU[32];
                        tmp1[3] = pJ1i[0]*ptrU[3]+pJ1i[6]*ptrU[9]+pJ1i[12]*ptrU[15]+pJ1i[18]*ptrU[21]+pJ1i[24]*ptrU[27]+pJ1i[30]*ptrU[33];
                        tmp1[4] = pJ1i[0]*ptrU[4]+pJ1i[6]*ptrU[10]+pJ1i[12]*ptrU[16]+pJ1i[18]*ptrU[22]+pJ1i[24]*ptrU[28]+pJ1i[30]*ptrU[34];
                        tmp1[5] = pJ1i[0]*ptrU[5]+pJ1i[6]*ptrU[11]+pJ1i[12]*ptrU[17]+pJ1i[18]*ptrU[23]+pJ1i[24]*ptrU[29]+pJ1i[30]*ptrU[35];

                        tmp2[0] = pJ1i[1]*ptrU[0]+pJ1i[7]*ptrU[6]+pJ1i[13]*ptrU[12]+pJ1i[19]*ptrU[18]+pJ1i[25]*ptrU[24]+pJ1i[31]*ptrU[30];
                        tmp2[1] = pJ1i[1]*ptrU[1]+pJ1i[7]*ptrU[7]+pJ1i[13]*ptrU[13]+pJ1i[19]*ptrU[19]+pJ1i[25]*ptrU[25]+pJ1i[31]*ptrU[31];
                        tmp2[2] = pJ1i[1]*ptrU[2]+pJ1i[7]*ptrU[8]+pJ1i[13]*ptrU[14]+pJ1i[19]*ptrU[20]+pJ1i[25]*ptrU[26]+pJ1i[31]*ptrU[32];
                        tmp2[3] = pJ1i[1]*ptrU[3]+pJ1i[7]*ptrU[9]+pJ1i[13]*ptrU[15]+pJ1i[19]*ptrU[21]+pJ1i[25]*ptrU[27]+pJ1i[31]*ptrU[33];
                        tmp2[4] = pJ1i[1]*ptrU[4]+pJ1i[7]*ptrU[10]+pJ1i[13]*ptrU[16]+pJ1i[19]*ptrU[22]+pJ1i[25]*ptrU[28]+pJ1i[31]*ptrU[34];
                        tmp2[5] = pJ1i[1]*ptrU[5]+pJ1i[7]*ptrU[11]+pJ1i[13]*ptrU[17]+pJ1i[19]*ptrU[23]+pJ1i[25]*ptrU[29]+pJ1i[31]*ptrU[35];

                        tmp3[0] = pJ1i[2]*ptrU[0]+pJ1i[8]*ptrU[6]+pJ1i[14]*ptrU[12]+pJ1i[20]*ptrU[18]+pJ1i[26]*ptrU[24]+pJ1i[32]*ptrU[30];
                        tmp3[1] = pJ1i[2]*ptrU[1]+pJ1i[8]*ptrU[7]+pJ1i[14]*ptrU[13]+pJ1i[20]*ptrU[19]+pJ1i[26]*ptrU[25]+pJ1i[32]*ptrU[31];
                        tmp3[2] = pJ1i[2]*ptrU[2]+pJ1i[8]*ptrU[8]+pJ1i[14]*ptrU[14]+pJ1i[20]*ptrU[20]+pJ1i[26]*ptrU[26]+pJ1i[32]*ptrU[32];
                        tmp3[3] = pJ1i[2]*ptrU[3]+pJ1i[8]*ptrU[9]+pJ1i[14]*ptrU[15]+pJ1i[20]*ptrU[21]+pJ1i[26]*ptrU[27]+pJ1i[32]*ptrU[33];
                        tmp3[4] = pJ1i[2]*ptrU[4]+pJ1i[8]*ptrU[10]+pJ1i[14]*ptrU[16]+pJ1i[20]*ptrU[22]+pJ1i[26]*ptrU[28]+pJ1i[32]*ptrU[34];
                        tmp3[5] = pJ1i[2]*ptrU[5]+pJ1i[8]*ptrU[11]+pJ1i[14]*ptrU[17]+pJ1i[20]*ptrU[23]+pJ1i[26]*ptrU[29]+pJ1i[32]*ptrU[35];

                        tmp4[0] = pJ1i[3]*ptrU[0]+pJ1i[9]*ptrU[6]+pJ1i[15]*ptrU[12]+pJ1i[21]*ptrU[18]+pJ1i[27]*ptrU[24]+pJ1i[33]*ptrU[30];
                        tmp4[1] = pJ1i[3]*ptrU[1]+pJ1i[9]*ptrU[7]+pJ1i[15]*ptrU[13]+pJ1i[21]*ptrU[19]+pJ1i[27]*ptrU[25]+pJ1i[33]*ptrU[31];
                        tmp4[2] = pJ1i[3]*ptrU[2]+pJ1i[9]*ptrU[8]+pJ1i[15]*ptrU[14]+pJ1i[21]*ptrU[20]+pJ1i[27]*ptrU[26]+pJ1i[33]*ptrU[32];
                        tmp4[3] = pJ1i[3]*ptrU[3]+pJ1i[9]*ptrU[9]+pJ1i[15]*ptrU[15]+pJ1i[21]*ptrU[21]+pJ1i[27]*ptrU[27]+pJ1i[33]*ptrU[33];
                        tmp4[4] = pJ1i[3]*ptrU[4]+pJ1i[9]*ptrU[10]+pJ1i[15]*ptrU[16]+pJ1i[21]*ptrU[22]+pJ1i[27]*ptrU[28]+pJ1i[33]*ptrU[34];
                        tmp4[5] = pJ1i[3]*ptrU[5]+pJ1i[9]*ptrU[11]+pJ1i[15]*ptrU[17]+pJ1i[21]*ptrU[23]+pJ1i[27]*ptrU[29]+pJ1i[33]*ptrU[35];

                        tmp5[0] = pJ1i[4]*ptrU[0]+pJ1i[10]*ptrU[6]+pJ1i[16]*ptrU[12]+pJ1i[22]*ptrU[18]+pJ1i[28]*ptrU[24]+pJ1i[34]*ptrU[30];
                        tmp5[1] = pJ1i[4]*ptrU[1]+pJ1i[10]*ptrU[7]+pJ1i[16]*ptrU[13]+pJ1i[22]*ptrU[19]+pJ1i[28]*ptrU[25]+pJ1i[34]*ptrU[31];
                        tmp5[2] = pJ1i[4]*ptrU[2]+pJ1i[10]*ptrU[8]+pJ1i[16]*ptrU[14]+pJ1i[22]*ptrU[20]+pJ1i[28]*ptrU[26]+pJ1i[34]*ptrU[32];
                        tmp5[3] = pJ1i[4]*ptrU[3]+pJ1i[10]*ptrU[9]+pJ1i[16]*ptrU[15]+pJ1i[22]*ptrU[21]+pJ1i[28]*ptrU[27]+pJ1i[34]*ptrU[33];
                        tmp5[4] = pJ1i[4]*ptrU[4]+pJ1i[10]*ptrU[10]+pJ1i[16]*ptrU[16]+pJ1i[22]*ptrU[22]+pJ1i[28]*ptrU[28]+pJ1i[34]*ptrU[34];
                        tmp5[5] = pJ1i[4]*ptrU[5]+pJ1i[10]*ptrU[11]+pJ1i[16]*ptrU[17]+pJ1i[22]*ptrU[23]+pJ1i[28]*ptrU[29]+pJ1i[34]*ptrU[35];

                        tmp6[0] = pJ1i[5]*ptrU[0]+pJ1i[11]*ptrU[6]+pJ1i[17]*ptrU[12]+pJ1i[23]*ptrU[18]+pJ1i[29]*ptrU[24]+pJ1i[35]*ptrU[30];
                        tmp6[1] = pJ1i[5]*ptrU[1]+pJ1i[11]*ptrU[7]+pJ1i[17]*ptrU[13]+pJ1i[23]*ptrU[19]+pJ1i[29]*ptrU[25]+pJ1i[35]*ptrU[31];
                        tmp6[2] = pJ1i[5]*ptrU[2]+pJ1i[11]*ptrU[8]+pJ1i[17]*ptrU[14]+pJ1i[23]*ptrU[20]+pJ1i[29]*ptrU[26]+pJ1i[35]*ptrU[32];
                        tmp6[3] = pJ1i[5]*ptrU[3]+pJ1i[11]*ptrU[9]+pJ1i[17]*ptrU[15]+pJ1i[23]*ptrU[21]+pJ1i[29]*ptrU[27]+pJ1i[35]*ptrU[33];
                        tmp6[4] = pJ1i[5]*ptrU[4]+pJ1i[11]*ptrU[10]+pJ1i[17]*ptrU[16]+pJ1i[23]*ptrU[22]+pJ1i[29]*ptrU[28]+pJ1i[35]*ptrU[34];
                        tmp6[5] = pJ1i[5]*ptrU[5]+pJ1i[11]*ptrU[11]+pJ1i[17]*ptrU[17]+pJ1i[23]*ptrU[23]+pJ1i[29]*ptrU[29]+pJ1i[35]*ptrU[35];

                        if ( m_Ui[i] == posID )
                        {
                                ptmp = newU+m_Uj[i]*(6*6);
                        }
                        else if ( m_Uj[i] == posID )
                        {
                                ptmp = newU+m_Ui[i]*(6*6);
                        }
						else if ( m_Ui[i] == posID2 )
                        {
                                ptmp = newU+(m+m_Uj[i])*(6*6);
                        }
                        else if ( m_Uj[i] == posID2 )
                        {
                                ptmp = newU+(m+m_Ui[i])*(6*6);
                        }
                        else
                        {
                                ptmp = newU + n_newU*6*6;
                                m_nUi[n_newU] = m_Ui[i];
                                m_nUj[n_newU] = m_Uj[i];
                                n_newU += 1;
                        }

                        // ttmp = J1(i)^T*I(j,k)*J1(j)
                        ptmp[0] += tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[6]+tmp1[2]*pJ1j[12]+tmp1[3]*pJ1j[18]+tmp1[4]*pJ1j[24]+tmp1[5]*pJ1j[30];
                        ptmp[1] += tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[7]+tmp1[2]*pJ1j[13]+tmp1[3]*pJ1j[19]+tmp1[4]*pJ1j[25]+tmp1[5]*pJ1j[31];
                        ptmp[2] += tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[8]+tmp1[2]*pJ1j[14]+tmp1[3]*pJ1j[20]+tmp1[4]*pJ1j[26]+tmp1[5]*pJ1j[32];
                        ptmp[3] += tmp1[0]*pJ1j[3]+tmp1[1]*pJ1j[9]+tmp1[2]*pJ1j[15]+tmp1[3]*pJ1j[21]+tmp1[4]*pJ1j[27]+tmp1[5]*pJ1j[33];
                        ptmp[4] += tmp1[0]*pJ1j[4]+tmp1[1]*pJ1j[10]+tmp1[2]*pJ1j[16]+tmp1[3]*pJ1j[22]+tmp1[4]*pJ1j[28]+tmp1[5]*pJ1j[34];
                        ptmp[5] += tmp1[0]*pJ1j[5]+tmp1[1]*pJ1j[11]+tmp1[2]*pJ1j[17]+tmp1[3]*pJ1j[23]+tmp1[4]*pJ1j[29]+tmp1[5]*pJ1j[35];

                        ptmp[6] += tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[6]+tmp2[2]*pJ1j[12]+tmp2[3]*pJ1j[18]+tmp2[4]*pJ1j[24]+tmp2[5]*pJ1j[30];
                        ptmp[7] += tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[7]+tmp2[2]*pJ1j[13]+tmp2[3]*pJ1j[19]+tmp2[4]*pJ1j[25]+tmp2[5]*pJ1j[31];
                        ptmp[8] += tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[8]+tmp2[2]*pJ1j[14]+tmp2[3]*pJ1j[20]+tmp2[4]*pJ1j[26]+tmp2[5]*pJ1j[32];
                        ptmp[9] += tmp2[0]*pJ1j[3]+tmp2[1]*pJ1j[9]+tmp2[2]*pJ1j[15]+tmp2[3]*pJ1j[21]+tmp2[4]*pJ1j[27]+tmp2[5]*pJ1j[33];
                        ptmp[10] += tmp2[0]*pJ1j[4]+tmp2[1]*pJ1j[10]+tmp2[2]*pJ1j[16]+tmp2[3]*pJ1j[22]+tmp2[4]*pJ1j[28]+tmp2[5]*pJ1j[34];
                        ptmp[11] += tmp2[0]*pJ1j[5]+tmp2[1]*pJ1j[11]+tmp2[2]*pJ1j[17]+tmp2[3]*pJ1j[23]+tmp2[4]*pJ1j[29]+tmp2[5]*pJ1j[35];

                        ptmp[12] += tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[6]+tmp3[2]*pJ1j[12]+tmp3[3]*pJ1j[18]+tmp3[4]*pJ1j[24]+tmp3[5]*pJ1j[30];
                        ptmp[13] += tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[7]+tmp3[2]*pJ1j[13]+tmp3[3]*pJ1j[19]+tmp3[4]*pJ1j[25]+tmp3[5]*pJ1j[31];
                        ptmp[14] += tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[8]+tmp3[2]*pJ1j[14]+tmp3[3]*pJ1j[20]+tmp3[4]*pJ1j[26]+tmp3[5]*pJ1j[32];
                        ptmp[15] += tmp3[0]*pJ1j[3]+tmp3[1]*pJ1j[9]+tmp3[2]*pJ1j[15]+tmp3[3]*pJ1j[21]+tmp3[4]*pJ1j[27]+tmp3[5]*pJ1j[33];
                        ptmp[16] += tmp3[0]*pJ1j[4]+tmp3[1]*pJ1j[10]+tmp3[2]*pJ1j[16]+tmp3[3]*pJ1j[22]+tmp3[4]*pJ1j[28]+tmp3[5]*pJ1j[34];
                        ptmp[17] += tmp3[0]*pJ1j[5]+tmp3[1]*pJ1j[11]+tmp3[2]*pJ1j[17]+tmp3[3]*pJ1j[23]+tmp3[4]*pJ1j[29]+tmp3[5]*pJ1j[35];

                        ptmp[18] += tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[6]+tmp4[2]*pJ1j[12]+tmp4[3]*pJ1j[18]+tmp4[4]*pJ1j[24]+tmp4[5]*pJ1j[30];
                        ptmp[19] += tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[7]+tmp4[2]*pJ1j[13]+tmp4[3]*pJ1j[19]+tmp4[4]*pJ1j[25]+tmp4[5]*pJ1j[31];
                        ptmp[20] += tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[8]+tmp4[2]*pJ1j[14]+tmp4[3]*pJ1j[20]+tmp4[4]*pJ1j[26]+tmp4[5]*pJ1j[32];
                        ptmp[21] += tmp4[0]*pJ1j[3]+tmp4[1]*pJ1j[9]+tmp4[2]*pJ1j[15]+tmp4[3]*pJ1j[21]+tmp4[4]*pJ1j[27]+tmp4[5]*pJ1j[33];
                        ptmp[22] += tmp4[0]*pJ1j[4]+tmp4[1]*pJ1j[10]+tmp4[2]*pJ1j[16]+tmp4[3]*pJ1j[22]+tmp4[4]*pJ1j[28]+tmp4[5]*pJ1j[34];
                        ptmp[23] += tmp4[0]*pJ1j[5]+tmp4[1]*pJ1j[11]+tmp4[2]*pJ1j[17]+tmp4[3]*pJ1j[23]+tmp4[4]*pJ1j[29]+tmp4[5]*pJ1j[35];

                        ptmp[24] += tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[6]+tmp5[2]*pJ1j[12]+tmp5[3]*pJ1j[18]+tmp5[4]*pJ1j[24]+tmp5[5]*pJ1j[30];
                        ptmp[25] += tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[7]+tmp5[2]*pJ1j[13]+tmp5[3]*pJ1j[19]+tmp5[4]*pJ1j[25]+tmp5[5]*pJ1j[31];
                        ptmp[26] += tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[8]+tmp5[2]*pJ1j[14]+tmp5[3]*pJ1j[20]+tmp5[4]*pJ1j[26]+tmp5[5]*pJ1j[32];
                        ptmp[27] += tmp5[0]*pJ1j[3]+tmp5[1]*pJ1j[9]+tmp5[2]*pJ1j[15]+tmp5[3]*pJ1j[21]+tmp5[4]*pJ1j[27]+tmp5[5]*pJ1j[33];
                        ptmp[28] += tmp5[0]*pJ1j[4]+tmp5[1]*pJ1j[10]+tmp5[2]*pJ1j[16]+tmp5[3]*pJ1j[22]+tmp5[4]*pJ1j[28]+tmp5[5]*pJ1j[34];
                        ptmp[29] += tmp5[0]*pJ1j[5]+tmp5[1]*pJ1j[11]+tmp5[2]*pJ1j[17]+tmp5[3]*pJ1j[23]+tmp5[4]*pJ1j[29]+tmp5[5]*pJ1j[35];

                        ptmp[30] += tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[6]+tmp6[2]*pJ1j[12]+tmp6[3]*pJ1j[18]+tmp6[4]*pJ1j[24]+tmp6[5]*pJ1j[30];
                        ptmp[31] += tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[7]+tmp6[2]*pJ1j[13]+tmp6[3]*pJ1j[19]+tmp6[4]*pJ1j[25]+tmp6[5]*pJ1j[31];
                        ptmp[32] += tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[8]+tmp6[2]*pJ1j[14]+tmp6[3]*pJ1j[20]+tmp6[4]*pJ1j[26]+tmp6[5]*pJ1j[32];
                        ptmp[33] += tmp6[0]*pJ1j[3]+tmp6[1]*pJ1j[9]+tmp6[2]*pJ1j[15]+tmp6[3]*pJ1j[21]+tmp6[4]*pJ1j[27]+tmp6[5]*pJ1j[33];
                        ptmp[34] += tmp6[0]*pJ1j[4]+tmp6[1]*pJ1j[10]+tmp6[2]*pJ1j[16]+tmp6[3]*pJ1j[22]+tmp6[4]*pJ1j[28]+tmp6[5]*pJ1j[34];
                        ptmp[35] += tmp6[0]*pJ1j[5]+tmp6[1]*pJ1j[11]+tmp6[2]*pJ1j[17]+tmp6[3]*pJ1j[23]+tmp6[4]*pJ1j[29]+tmp6[5]*pJ1j[35];

                        //Algorithm Line 6

                        ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12]+tmp1[3]*pJ2j[18]+tmp1[4]*pJ2j[24]+tmp1[5]*pJ2j[30];
                        ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13]+tmp1[3]*pJ2j[19]+tmp1[4]*pJ2j[25]+tmp1[5]*pJ2j[31];
                        ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14]+tmp1[3]*pJ2j[20]+tmp1[4]*pJ2j[26]+tmp1[5]*pJ2j[32];
                        ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15]+tmp1[3]*pJ2j[21]+tmp1[4]*pJ2j[27]+tmp1[5]*pJ2j[33];
                        ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16]+tmp1[3]*pJ2j[22]+tmp1[4]*pJ2j[28]+tmp1[5]*pJ2j[34];
                        ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17]+tmp1[3]*pJ2j[23]+tmp1[4]*pJ2j[29]+tmp1[5]*pJ2j[35];

                        ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12]+tmp2[3]*pJ2j[18]+tmp2[4]*pJ2j[24]+tmp2[5]*pJ2j[30];
                        ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13]+tmp2[3]*pJ2j[19]+tmp2[4]*pJ2j[25]+tmp2[5]*pJ2j[31];
                        ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14]+tmp2[3]*pJ2j[20]+tmp2[4]*pJ2j[26]+tmp2[5]*pJ2j[32];
                        ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15]+tmp2[3]*pJ2j[21]+tmp2[4]*pJ2j[27]+tmp2[5]*pJ2j[33];
                        ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16]+tmp2[3]*pJ2j[22]+tmp2[4]*pJ2j[28]+tmp2[5]*pJ2j[34];
                        ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17]+tmp2[3]*pJ2j[23]+tmp2[4]*pJ2j[29]+tmp2[5]*pJ2j[35];

                        ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12]+tmp3[3]*pJ2j[18]+tmp3[4]*pJ2j[24]+tmp3[5]*pJ2j[30];
                        ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13]+tmp3[3]*pJ2j[19]+tmp3[4]*pJ2j[25]+tmp3[5]*pJ2j[31];
                        ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14]+tmp3[3]*pJ2j[20]+tmp3[4]*pJ2j[26]+tmp3[5]*pJ2j[32];
                        ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15]+tmp3[3]*pJ2j[21]+tmp3[4]*pJ2j[27]+tmp3[5]*pJ2j[33];
                        ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16]+tmp3[3]*pJ2j[22]+tmp3[4]*pJ2j[28]+tmp3[5]*pJ2j[34];
                        ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17]+tmp3[3]*pJ2j[23]+tmp3[4]*pJ2j[29]+tmp3[5]*pJ2j[35];

                        ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12]+tmp4[3]*pJ2j[18]+tmp4[4]*pJ2j[24]+tmp4[5]*pJ2j[30];
                        ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13]+tmp4[3]*pJ2j[19]+tmp4[4]*pJ2j[25]+tmp4[5]*pJ2j[31];
                        ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14]+tmp4[3]*pJ2j[20]+tmp4[4]*pJ2j[26]+tmp4[5]*pJ2j[32];
                        ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15]+tmp4[3]*pJ2j[21]+tmp4[4]*pJ2j[27]+tmp4[5]*pJ2j[33];
                        ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16]+tmp4[3]*pJ2j[22]+tmp4[4]*pJ2j[28]+tmp4[5]*pJ2j[34];
                        ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17]+tmp4[3]*pJ2j[23]+tmp4[4]*pJ2j[29]+tmp4[5]*pJ2j[35];

                        ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12]+tmp5[3]*pJ2j[18]+tmp5[4]*pJ2j[24]+tmp5[5]*pJ2j[30];
                        ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13]+tmp5[3]*pJ2j[19]+tmp5[4]*pJ2j[25]+tmp5[5]*pJ2j[31];
                        ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14]+tmp5[3]*pJ2j[20]+tmp5[4]*pJ2j[26]+tmp5[5]*pJ2j[32];
                        ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15]+tmp5[3]*pJ2j[21]+tmp5[4]*pJ2j[27]+tmp5[5]*pJ2j[33];
                        ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16]+tmp5[3]*pJ2j[22]+tmp5[4]*pJ2j[28]+tmp5[5]*pJ2j[34];
                        ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17]+tmp5[3]*pJ2j[23]+tmp5[4]*pJ2j[29]+tmp5[5]*pJ2j[35];

                        ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12]+tmp6[3]*pJ2j[18]+tmp6[4]*pJ2j[24]+tmp6[5]*pJ2j[30];
                        ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13]+tmp6[3]*pJ2j[19]+tmp6[4]*pJ2j[25]+tmp6[5]*pJ2j[31];
                        ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14]+tmp6[3]*pJ2j[20]+tmp6[4]*pJ2j[26]+tmp6[5]*pJ2j[32];
                        ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15]+tmp6[3]*pJ2j[21]+tmp6[4]*pJ2j[27]+tmp6[5]*pJ2j[33];
                        ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16]+tmp6[3]*pJ2j[22]+tmp6[4]*pJ2j[28]+tmp6[5]*pJ2j[34];
                        ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17]+tmp6[3]*pJ2j[23]+tmp6[4]*pJ2j[29]+tmp6[5]*pJ2j[35];

                        ptmp = newU+m_Ui[i]*(6*6);
                        if ( m_Ui[i] <= posID )
                        {
                                // ttmp = J1(i)^T*I(j,k)*J2(j)
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];                                   
                        }
                        if ( m_Ui[i] >= posID && m_Ui[i] != m_Uj[i])
                        {
                                //if ( m_Ui[i] != m_Uj[i] )
                                //{
                                // ttmp = ttmp^T
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                        }

						///////////////////////////////////////////////////////////////////////
						//J1i^T*I(i,j)*J3j
                        ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12]+tmp1[3]*pJ3j[18]+tmp1[4]*pJ3j[24]+tmp1[5]*pJ3j[30];
                        ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13]+tmp1[3]*pJ3j[19]+tmp1[4]*pJ3j[25]+tmp1[5]*pJ3j[31];
                        ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14]+tmp1[3]*pJ3j[20]+tmp1[4]*pJ3j[26]+tmp1[5]*pJ3j[32];
                        ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15]+tmp1[3]*pJ3j[21]+tmp1[4]*pJ3j[27]+tmp1[5]*pJ3j[33];
                        ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16]+tmp1[3]*pJ3j[22]+tmp1[4]*pJ3j[28]+tmp1[5]*pJ3j[34];
                        ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17]+tmp1[3]*pJ3j[23]+tmp1[4]*pJ3j[29]+tmp1[5]*pJ3j[35];

                        ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12]+tmp2[3]*pJ3j[18]+tmp2[4]*pJ3j[24]+tmp2[5]*pJ3j[30];
                        ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13]+tmp2[3]*pJ3j[19]+tmp2[4]*pJ3j[25]+tmp2[5]*pJ3j[31];
                        ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14]+tmp2[3]*pJ3j[20]+tmp2[4]*pJ3j[26]+tmp2[5]*pJ3j[32];
                        ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15]+tmp2[3]*pJ3j[21]+tmp2[4]*pJ3j[27]+tmp2[5]*pJ3j[33];
                        ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16]+tmp2[3]*pJ3j[22]+tmp2[4]*pJ3j[28]+tmp2[5]*pJ3j[34];
                        ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17]+tmp2[3]*pJ3j[23]+tmp2[4]*pJ3j[29]+tmp2[5]*pJ3j[35];

                        ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12]+tmp3[3]*pJ3j[18]+tmp3[4]*pJ3j[24]+tmp3[5]*pJ3j[30];
                        ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13]+tmp3[3]*pJ3j[19]+tmp3[4]*pJ3j[25]+tmp3[5]*pJ3j[31];
                        ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14]+tmp3[3]*pJ3j[20]+tmp3[4]*pJ3j[26]+tmp3[5]*pJ3j[32];
                        ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15]+tmp3[3]*pJ3j[21]+tmp3[4]*pJ3j[27]+tmp3[5]*pJ3j[33];
                        ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16]+tmp3[3]*pJ3j[22]+tmp3[4]*pJ3j[28]+tmp3[5]*pJ3j[34];
                        ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17]+tmp3[3]*pJ3j[23]+tmp3[4]*pJ3j[29]+tmp3[5]*pJ3j[35];

                        ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12]+tmp4[3]*pJ3j[18]+tmp4[4]*pJ3j[24]+tmp4[5]*pJ3j[30];
                        ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13]+tmp4[3]*pJ3j[19]+tmp4[4]*pJ3j[25]+tmp4[5]*pJ3j[31];
                        ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14]+tmp4[3]*pJ3j[20]+tmp4[4]*pJ3j[26]+tmp4[5]*pJ3j[32];
                        ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15]+tmp4[3]*pJ3j[21]+tmp4[4]*pJ3j[27]+tmp4[5]*pJ3j[33];
                        ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16]+tmp4[3]*pJ3j[22]+tmp4[4]*pJ3j[28]+tmp4[5]*pJ3j[34];
                        ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17]+tmp4[3]*pJ3j[23]+tmp4[4]*pJ3j[29]+tmp4[5]*pJ3j[35];

                        ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12]+tmp5[3]*pJ3j[18]+tmp5[4]*pJ3j[24]+tmp5[5]*pJ3j[30];
                        ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13]+tmp5[3]*pJ3j[19]+tmp5[4]*pJ3j[25]+tmp5[5]*pJ3j[31];
                        ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14]+tmp5[3]*pJ3j[20]+tmp5[4]*pJ3j[26]+tmp5[5]*pJ3j[32];
                        ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15]+tmp5[3]*pJ3j[21]+tmp5[4]*pJ3j[27]+tmp5[5]*pJ3j[33];
                        ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16]+tmp5[3]*pJ3j[22]+tmp5[4]*pJ3j[28]+tmp5[5]*pJ3j[34];
                        ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17]+tmp5[3]*pJ3j[23]+tmp5[4]*pJ3j[29]+tmp5[5]*pJ3j[35];

                        ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12]+tmp6[3]*pJ3j[18]+tmp6[4]*pJ3j[24]+tmp6[5]*pJ3j[30];
                        ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13]+tmp6[3]*pJ3j[19]+tmp6[4]*pJ3j[25]+tmp6[5]*pJ3j[31];
                        ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14]+tmp6[3]*pJ3j[20]+tmp6[4]*pJ3j[26]+tmp6[5]*pJ3j[32];
                        ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15]+tmp6[3]*pJ3j[21]+tmp6[4]*pJ3j[27]+tmp6[5]*pJ3j[33];
                        ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16]+tmp6[3]*pJ3j[22]+tmp6[4]*pJ3j[28]+tmp6[5]*pJ3j[34];
                        ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17]+tmp6[3]*pJ3j[23]+tmp6[4]*pJ3j[29]+tmp6[5]*pJ3j[35];

                        ptmp = newU+(m+m_Ui[i])*(6*6);
                        if ( m_Ui[i] <= posID2 )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];                                   
                        }
                        if ( m_Ui[i] >= posID2 && m_Ui[i] != m_Uj[i])
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                        }
						////////////////////////////////////////////////////////////
						//J3i*I(i,j)*J3j
                        tmp1[0] = pJ3i[0]*ptrU[0]+pJ3i[6]*ptrU[6]+pJ3i[12]*ptrU[12]+pJ3i[18]*ptrU[18]+pJ3i[24]*ptrU[24]+pJ3i[30]*ptrU[30];
                        tmp1[1] = pJ3i[0]*ptrU[1]+pJ3i[6]*ptrU[7]+pJ3i[12]*ptrU[13]+pJ3i[18]*ptrU[19]+pJ3i[24]*ptrU[25]+pJ3i[30]*ptrU[31];
                        tmp1[2] = pJ3i[0]*ptrU[2]+pJ3i[6]*ptrU[8]+pJ3i[12]*ptrU[14]+pJ3i[18]*ptrU[20]+pJ3i[24]*ptrU[26]+pJ3i[30]*ptrU[32];
                        tmp1[3] = pJ3i[0]*ptrU[3]+pJ3i[6]*ptrU[9]+pJ3i[12]*ptrU[15]+pJ3i[18]*ptrU[21]+pJ3i[24]*ptrU[27]+pJ3i[30]*ptrU[33];
                        tmp1[4] = pJ3i[0]*ptrU[4]+pJ3i[6]*ptrU[10]+pJ3i[12]*ptrU[16]+pJ3i[18]*ptrU[22]+pJ3i[24]*ptrU[28]+pJ3i[30]*ptrU[34];
                        tmp1[5] = pJ3i[0]*ptrU[5]+pJ3i[6]*ptrU[11]+pJ3i[12]*ptrU[17]+pJ3i[18]*ptrU[23]+pJ3i[24]*ptrU[29]+pJ3i[30]*ptrU[35];

                        tmp2[0] = pJ3i[1]*ptrU[0]+pJ3i[7]*ptrU[6]+pJ3i[13]*ptrU[12]+pJ3i[19]*ptrU[18]+pJ3i[25]*ptrU[24]+pJ3i[31]*ptrU[30];
                        tmp2[1] = pJ3i[1]*ptrU[1]+pJ3i[7]*ptrU[7]+pJ3i[13]*ptrU[13]+pJ3i[19]*ptrU[19]+pJ3i[25]*ptrU[25]+pJ3i[31]*ptrU[31];
                        tmp2[2] = pJ3i[1]*ptrU[2]+pJ3i[7]*ptrU[8]+pJ3i[13]*ptrU[14]+pJ3i[19]*ptrU[20]+pJ3i[25]*ptrU[26]+pJ3i[31]*ptrU[32];
                        tmp2[3] = pJ3i[1]*ptrU[3]+pJ3i[7]*ptrU[9]+pJ3i[13]*ptrU[15]+pJ3i[19]*ptrU[21]+pJ3i[25]*ptrU[27]+pJ3i[31]*ptrU[33];
                        tmp2[4] = pJ3i[1]*ptrU[4]+pJ3i[7]*ptrU[10]+pJ3i[13]*ptrU[16]+pJ3i[19]*ptrU[22]+pJ3i[25]*ptrU[28]+pJ3i[31]*ptrU[34];
                        tmp2[5] = pJ3i[1]*ptrU[5]+pJ3i[7]*ptrU[11]+pJ3i[13]*ptrU[17]+pJ3i[19]*ptrU[23]+pJ3i[25]*ptrU[29]+pJ3i[31]*ptrU[35];

                        tmp3[0] = pJ3i[2]*ptrU[0]+pJ3i[8]*ptrU[6]+pJ3i[14]*ptrU[12]+pJ3i[20]*ptrU[18]+pJ3i[26]*ptrU[24]+pJ3i[32]*ptrU[30];
                        tmp3[1] = pJ3i[2]*ptrU[1]+pJ3i[8]*ptrU[7]+pJ3i[14]*ptrU[13]+pJ3i[20]*ptrU[19]+pJ3i[26]*ptrU[25]+pJ3i[32]*ptrU[31];
                        tmp3[2] = pJ3i[2]*ptrU[2]+pJ3i[8]*ptrU[8]+pJ3i[14]*ptrU[14]+pJ3i[20]*ptrU[20]+pJ3i[26]*ptrU[26]+pJ3i[32]*ptrU[32];
                        tmp3[3] = pJ3i[2]*ptrU[3]+pJ3i[8]*ptrU[9]+pJ3i[14]*ptrU[15]+pJ3i[20]*ptrU[21]+pJ3i[26]*ptrU[27]+pJ3i[32]*ptrU[33];
                        tmp3[4] = pJ3i[2]*ptrU[4]+pJ3i[8]*ptrU[10]+pJ3i[14]*ptrU[16]+pJ3i[20]*ptrU[22]+pJ3i[26]*ptrU[28]+pJ3i[32]*ptrU[34];
                        tmp3[5] = pJ3i[2]*ptrU[5]+pJ3i[8]*ptrU[11]+pJ3i[14]*ptrU[17]+pJ3i[20]*ptrU[23]+pJ3i[26]*ptrU[29]+pJ3i[32]*ptrU[35];

                        tmp4[0] = pJ3i[3]*ptrU[0]+pJ3i[9]*ptrU[6]+pJ3i[15]*ptrU[12]+pJ3i[21]*ptrU[18]+pJ3i[27]*ptrU[24]+pJ3i[33]*ptrU[30];
                        tmp4[1] = pJ3i[3]*ptrU[1]+pJ3i[9]*ptrU[7]+pJ3i[15]*ptrU[13]+pJ3i[21]*ptrU[19]+pJ3i[27]*ptrU[25]+pJ3i[33]*ptrU[31];
                        tmp4[2] = pJ3i[3]*ptrU[2]+pJ3i[9]*ptrU[8]+pJ3i[15]*ptrU[14]+pJ3i[21]*ptrU[20]+pJ3i[27]*ptrU[26]+pJ3i[33]*ptrU[32];
                        tmp4[3] = pJ3i[3]*ptrU[3]+pJ3i[9]*ptrU[9]+pJ3i[15]*ptrU[15]+pJ3i[21]*ptrU[21]+pJ3i[27]*ptrU[27]+pJ3i[33]*ptrU[33];
                        tmp4[4] = pJ3i[3]*ptrU[4]+pJ3i[9]*ptrU[10]+pJ3i[15]*ptrU[16]+pJ3i[21]*ptrU[22]+pJ3i[27]*ptrU[28]+pJ3i[33]*ptrU[34];
                        tmp4[5] = pJ3i[3]*ptrU[5]+pJ3i[9]*ptrU[11]+pJ3i[15]*ptrU[17]+pJ3i[21]*ptrU[23]+pJ3i[27]*ptrU[29]+pJ3i[33]*ptrU[35];

                        tmp5[0] = pJ3i[4]*ptrU[0]+pJ3i[10]*ptrU[6]+pJ3i[16]*ptrU[12]+pJ3i[22]*ptrU[18]+pJ3i[28]*ptrU[24]+pJ3i[34]*ptrU[30];
                        tmp5[1] = pJ3i[4]*ptrU[1]+pJ3i[10]*ptrU[7]+pJ3i[16]*ptrU[13]+pJ3i[22]*ptrU[19]+pJ3i[28]*ptrU[25]+pJ3i[34]*ptrU[31];
                        tmp5[2] = pJ3i[4]*ptrU[2]+pJ3i[10]*ptrU[8]+pJ3i[16]*ptrU[14]+pJ3i[22]*ptrU[20]+pJ3i[28]*ptrU[26]+pJ3i[34]*ptrU[32];
                        tmp5[3] = pJ3i[4]*ptrU[3]+pJ3i[10]*ptrU[9]+pJ3i[16]*ptrU[15]+pJ3i[22]*ptrU[21]+pJ3i[28]*ptrU[27]+pJ3i[34]*ptrU[33];
                        tmp5[4] = pJ3i[4]*ptrU[4]+pJ3i[10]*ptrU[10]+pJ3i[16]*ptrU[16]+pJ3i[22]*ptrU[22]+pJ3i[28]*ptrU[28]+pJ3i[34]*ptrU[34];
                        tmp5[5] = pJ3i[4]*ptrU[5]+pJ3i[10]*ptrU[11]+pJ3i[16]*ptrU[17]+pJ3i[22]*ptrU[23]+pJ3i[28]*ptrU[29]+pJ3i[34]*ptrU[35];

                        tmp6[0] = pJ3i[5]*ptrU[0]+pJ3i[11]*ptrU[6]+pJ3i[17]*ptrU[12]+pJ3i[23]*ptrU[18]+pJ3i[29]*ptrU[24]+pJ3i[35]*ptrU[30];
                        tmp6[1] = pJ3i[5]*ptrU[1]+pJ3i[11]*ptrU[7]+pJ3i[17]*ptrU[13]+pJ3i[23]*ptrU[19]+pJ3i[29]*ptrU[25]+pJ3i[35]*ptrU[31];
                        tmp6[2] = pJ3i[5]*ptrU[2]+pJ3i[11]*ptrU[8]+pJ3i[17]*ptrU[14]+pJ3i[23]*ptrU[20]+pJ3i[29]*ptrU[26]+pJ3i[35]*ptrU[32];
                        tmp6[3] = pJ3i[5]*ptrU[3]+pJ3i[11]*ptrU[9]+pJ3i[17]*ptrU[15]+pJ3i[23]*ptrU[21]+pJ3i[29]*ptrU[27]+pJ3i[35]*ptrU[33];
                        tmp6[4] = pJ3i[5]*ptrU[4]+pJ3i[11]*ptrU[10]+pJ3i[17]*ptrU[16]+pJ3i[23]*ptrU[22]+pJ3i[29]*ptrU[28]+pJ3i[35]*ptrU[34];
                        tmp6[5] = pJ3i[5]*ptrU[5]+pJ3i[11]*ptrU[11]+pJ3i[17]*ptrU[17]+pJ3i[23]*ptrU[23]+pJ3i[29]*ptrU[29]+pJ3i[35]*ptrU[35];

                        ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12]+tmp1[3]*pJ3j[18]+tmp1[4]*pJ3j[24]+tmp1[5]*pJ3j[30];
                        ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13]+tmp1[3]*pJ3j[19]+tmp1[4]*pJ3j[25]+tmp1[5]*pJ3j[31];
                        ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14]+tmp1[3]*pJ3j[20]+tmp1[4]*pJ3j[26]+tmp1[5]*pJ3j[32];
                        ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15]+tmp1[3]*pJ3j[21]+tmp1[4]*pJ3j[27]+tmp1[5]*pJ3j[33];
                        ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16]+tmp1[3]*pJ3j[22]+tmp1[4]*pJ3j[28]+tmp1[5]*pJ3j[34];
                        ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17]+tmp1[3]*pJ3j[23]+tmp1[4]*pJ3j[29]+tmp1[5]*pJ3j[35];

                        ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12]+tmp2[3]*pJ3j[18]+tmp2[4]*pJ3j[24]+tmp2[5]*pJ3j[30];
                        ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13]+tmp2[3]*pJ3j[19]+tmp2[4]*pJ3j[25]+tmp2[5]*pJ3j[31];
                        ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14]+tmp2[3]*pJ3j[20]+tmp2[4]*pJ3j[26]+tmp2[5]*pJ3j[32];
                        ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15]+tmp2[3]*pJ3j[21]+tmp2[4]*pJ3j[27]+tmp2[5]*pJ3j[33];
                        ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16]+tmp2[3]*pJ3j[22]+tmp2[4]*pJ3j[28]+tmp2[5]*pJ3j[34];
                        ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17]+tmp2[3]*pJ3j[23]+tmp2[4]*pJ3j[29]+tmp2[5]*pJ3j[35];

                        ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12]+tmp3[3]*pJ3j[18]+tmp3[4]*pJ3j[24]+tmp3[5]*pJ3j[30];
                        ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13]+tmp3[3]*pJ3j[19]+tmp3[4]*pJ3j[25]+tmp3[5]*pJ3j[31];
                        ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14]+tmp3[3]*pJ3j[20]+tmp3[4]*pJ3j[26]+tmp3[5]*pJ3j[32];
                        ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15]+tmp3[3]*pJ3j[21]+tmp3[4]*pJ3j[27]+tmp3[5]*pJ3j[33];
                        ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16]+tmp3[3]*pJ3j[22]+tmp3[4]*pJ3j[28]+tmp3[5]*pJ3j[34];
                        ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17]+tmp3[3]*pJ3j[23]+tmp3[4]*pJ3j[29]+tmp3[5]*pJ3j[35];

                        ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12]+tmp4[3]*pJ3j[18]+tmp4[4]*pJ3j[24]+tmp4[5]*pJ3j[30];
                        ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13]+tmp4[3]*pJ3j[19]+tmp4[4]*pJ3j[25]+tmp4[5]*pJ3j[31];
                        ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14]+tmp4[3]*pJ3j[20]+tmp4[4]*pJ3j[26]+tmp4[5]*pJ3j[32];
                        ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15]+tmp4[3]*pJ3j[21]+tmp4[4]*pJ3j[27]+tmp4[5]*pJ3j[33];
                        ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16]+tmp4[3]*pJ3j[22]+tmp4[4]*pJ3j[28]+tmp4[5]*pJ3j[34];
                        ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17]+tmp4[3]*pJ3j[23]+tmp4[4]*pJ3j[29]+tmp4[5]*pJ3j[35];

                        ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12]+tmp5[3]*pJ3j[18]+tmp5[4]*pJ3j[24]+tmp5[5]*pJ3j[30];
                        ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13]+tmp5[3]*pJ3j[19]+tmp5[4]*pJ3j[25]+tmp5[5]*pJ3j[31];
                        ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14]+tmp5[3]*pJ3j[20]+tmp5[4]*pJ3j[26]+tmp5[5]*pJ3j[32];
                        ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15]+tmp5[3]*pJ3j[21]+tmp5[4]*pJ3j[27]+tmp5[5]*pJ3j[33];
                        ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16]+tmp5[3]*pJ3j[22]+tmp5[4]*pJ3j[28]+tmp5[5]*pJ3j[34];
                        ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17]+tmp5[3]*pJ3j[23]+tmp5[4]*pJ3j[29]+tmp5[5]*pJ3j[35];

                        ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12]+tmp6[3]*pJ3j[18]+tmp6[4]*pJ3j[24]+tmp6[5]*pJ3j[30];
                        ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13]+tmp6[3]*pJ3j[19]+tmp6[4]*pJ3j[25]+tmp6[5]*pJ3j[31];
                        ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14]+tmp6[3]*pJ3j[20]+tmp6[4]*pJ3j[26]+tmp6[5]*pJ3j[32];
                        ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15]+tmp6[3]*pJ3j[21]+tmp6[4]*pJ3j[27]+tmp6[5]*pJ3j[33];
                        ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16]+tmp6[3]*pJ3j[22]+tmp6[4]*pJ3j[28]+tmp6[5]*pJ3j[34];
                        ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17]+tmp6[3]*pJ3j[23]+tmp6[4]*pJ3j[29]+tmp6[5]*pJ3j[35];

                        ptmp = newU+(m+posID2)*(6*6);

                        ptmp[0] += ttmp1[0];
                        ptmp[1] += ttmp1[1];
                        ptmp[2] += ttmp1[2];
                        ptmp[3] += ttmp1[3];
                        ptmp[4] += ttmp1[4];
                        ptmp[5] += ttmp1[5];
                        ptmp[6] += ttmp2[0];
                        ptmp[7] += ttmp2[1];
                        ptmp[8] += ttmp2[2];
                        ptmp[9] += ttmp2[3];
                        ptmp[10] += ttmp2[4];
                        ptmp[11] += ttmp2[5];
                        ptmp[12] += ttmp3[0];
                        ptmp[13] += ttmp3[1];
                        ptmp[14] += ttmp3[2];
                        ptmp[15] += ttmp3[3];
                        ptmp[16] += ttmp3[4];
                        ptmp[17] += ttmp3[5];
                        ptmp[18] += ttmp4[0];
                        ptmp[19] += ttmp4[1];
                        ptmp[20] += ttmp4[2];
                        ptmp[21] += ttmp4[3];
                        ptmp[22] += ttmp4[4];
                        ptmp[23] += ttmp4[5];
                        ptmp[24] += ttmp5[0];
                        ptmp[25] += ttmp5[1];
                        ptmp[26] += ttmp5[2];
                        ptmp[27] += ttmp5[3];
                        ptmp[28] += ttmp5[4];
                        ptmp[29] += ttmp5[5];
                        ptmp[30] += ttmp6[0];
                        ptmp[31] += ttmp6[1];
                        ptmp[32] += ttmp6[2];
                        ptmp[33] += ttmp6[3];
                        ptmp[34] += ttmp6[4];
                        ptmp[35] += ttmp6[5];                                   


                        if ( m_Ui[i] != m_Uj[i] )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];           

                        }

                        //J3i^T*I(i,j)*J1j

                        ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[6]+tmp1[2]*pJ1j[12]+tmp1[3]*pJ1j[18]+tmp1[4]*pJ1j[24]+tmp1[5]*pJ1j[30];
                        ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[7]+tmp1[2]*pJ1j[13]+tmp1[3]*pJ1j[19]+tmp1[4]*pJ1j[25]+tmp1[5]*pJ1j[31];
                        ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[8]+tmp1[2]*pJ1j[14]+tmp1[3]*pJ1j[20]+tmp1[4]*pJ1j[26]+tmp1[5]*pJ1j[32];
                        ttmp1[3] = tmp1[0]*pJ1j[3]+tmp1[1]*pJ1j[9]+tmp1[2]*pJ1j[15]+tmp1[3]*pJ1j[21]+tmp1[4]*pJ1j[27]+tmp1[5]*pJ1j[33];
                        ttmp1[4] = tmp1[0]*pJ1j[4]+tmp1[1]*pJ1j[10]+tmp1[2]*pJ1j[16]+tmp1[3]*pJ1j[22]+tmp1[4]*pJ1j[28]+tmp1[5]*pJ1j[34];
                        ttmp1[5] = tmp1[0]*pJ1j[5]+tmp1[1]*pJ1j[11]+tmp1[2]*pJ1j[17]+tmp1[3]*pJ1j[23]+tmp1[4]*pJ1j[29]+tmp1[5]*pJ1j[35];

                        ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[6]+tmp2[2]*pJ1j[12]+tmp2[3]*pJ1j[18]+tmp2[4]*pJ1j[24]+tmp2[5]*pJ1j[30];
                        ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[7]+tmp2[2]*pJ1j[13]+tmp2[3]*pJ1j[19]+tmp2[4]*pJ1j[25]+tmp2[5]*pJ1j[31];
                        ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[8]+tmp2[2]*pJ1j[14]+tmp2[3]*pJ1j[20]+tmp2[4]*pJ1j[26]+tmp2[5]*pJ1j[32];
                        ttmp2[3] = tmp2[0]*pJ1j[3]+tmp2[1]*pJ1j[9]+tmp2[2]*pJ1j[15]+tmp2[3]*pJ1j[21]+tmp2[4]*pJ1j[27]+tmp2[5]*pJ1j[33];
                        ttmp2[4] = tmp2[0]*pJ1j[4]+tmp2[1]*pJ1j[10]+tmp2[2]*pJ1j[16]+tmp2[3]*pJ1j[22]+tmp2[4]*pJ1j[28]+tmp2[5]*pJ1j[34];
                        ttmp2[5] = tmp2[0]*pJ1j[5]+tmp2[1]*pJ1j[11]+tmp2[2]*pJ1j[17]+tmp2[3]*pJ1j[23]+tmp2[4]*pJ1j[29]+tmp2[5]*pJ1j[35];

                        ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[6]+tmp3[2]*pJ1j[12]+tmp3[3]*pJ1j[18]+tmp3[4]*pJ1j[24]+tmp3[5]*pJ1j[30];
                        ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[7]+tmp3[2]*pJ1j[13]+tmp3[3]*pJ1j[19]+tmp3[4]*pJ1j[25]+tmp3[5]*pJ1j[31];
                        ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[8]+tmp3[2]*pJ1j[14]+tmp3[3]*pJ1j[20]+tmp3[4]*pJ1j[26]+tmp3[5]*pJ1j[32];
                        ttmp3[3] = tmp3[0]*pJ1j[3]+tmp3[1]*pJ1j[9]+tmp3[2]*pJ1j[15]+tmp3[3]*pJ1j[21]+tmp3[4]*pJ1j[27]+tmp3[5]*pJ1j[33];
                        ttmp3[4] = tmp3[0]*pJ1j[4]+tmp3[1]*pJ1j[10]+tmp3[2]*pJ1j[16]+tmp3[3]*pJ1j[22]+tmp3[4]*pJ1j[28]+tmp3[5]*pJ1j[34];
                        ttmp3[5] = tmp3[0]*pJ1j[5]+tmp3[1]*pJ1j[11]+tmp3[2]*pJ1j[17]+tmp3[3]*pJ1j[23]+tmp3[4]*pJ1j[29]+tmp3[5]*pJ1j[35];

                        ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[6]+tmp4[2]*pJ1j[12]+tmp4[3]*pJ1j[18]+tmp4[4]*pJ1j[24]+tmp4[5]*pJ1j[30];
                        ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[7]+tmp4[2]*pJ1j[13]+tmp4[3]*pJ1j[19]+tmp4[4]*pJ1j[25]+tmp4[5]*pJ1j[31];
                        ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[8]+tmp4[2]*pJ1j[14]+tmp4[3]*pJ1j[20]+tmp4[4]*pJ1j[26]+tmp4[5]*pJ1j[32];
                        ttmp4[3] = tmp4[0]*pJ1j[3]+tmp4[1]*pJ1j[9]+tmp4[2]*pJ1j[15]+tmp4[3]*pJ1j[21]+tmp4[4]*pJ1j[27]+tmp4[5]*pJ1j[33];
                        ttmp4[4] = tmp4[0]*pJ1j[4]+tmp4[1]*pJ1j[10]+tmp4[2]*pJ1j[16]+tmp4[3]*pJ1j[22]+tmp4[4]*pJ1j[28]+tmp4[5]*pJ1j[34];
                        ttmp4[5] = tmp4[0]*pJ1j[5]+tmp4[1]*pJ1j[11]+tmp4[2]*pJ1j[17]+tmp4[3]*pJ1j[23]+tmp4[4]*pJ1j[29]+tmp4[5]*pJ1j[35];

                        ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[6]+tmp5[2]*pJ1j[12]+tmp5[3]*pJ1j[18]+tmp5[4]*pJ1j[24]+tmp5[5]*pJ1j[30];
                        ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[7]+tmp5[2]*pJ1j[13]+tmp5[3]*pJ1j[19]+tmp5[4]*pJ1j[25]+tmp5[5]*pJ1j[31];
                        ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[8]+tmp5[2]*pJ1j[14]+tmp5[3]*pJ1j[20]+tmp5[4]*pJ1j[26]+tmp5[5]*pJ1j[32];
                        ttmp5[3] = tmp5[0]*pJ1j[3]+tmp5[1]*pJ1j[9]+tmp5[2]*pJ1j[15]+tmp5[3]*pJ1j[21]+tmp5[4]*pJ1j[27]+tmp5[5]*pJ1j[33];
                        ttmp5[4] = tmp5[0]*pJ1j[4]+tmp5[1]*pJ1j[10]+tmp5[2]*pJ1j[16]+tmp5[3]*pJ1j[22]+tmp5[4]*pJ1j[28]+tmp5[5]*pJ1j[34];
                        ttmp5[5] = tmp5[0]*pJ1j[5]+tmp5[1]*pJ1j[11]+tmp5[2]*pJ1j[17]+tmp5[3]*pJ1j[23]+tmp5[4]*pJ1j[29]+tmp5[5]*pJ1j[35];

                        ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[6]+tmp6[2]*pJ1j[12]+tmp6[3]*pJ1j[18]+tmp6[4]*pJ1j[24]+tmp6[5]*pJ1j[30];
                        ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[7]+tmp6[2]*pJ1j[13]+tmp6[3]*pJ1j[19]+tmp6[4]*pJ1j[25]+tmp6[5]*pJ1j[31];
                        ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[8]+tmp6[2]*pJ1j[14]+tmp6[3]*pJ1j[20]+tmp6[4]*pJ1j[26]+tmp6[5]*pJ1j[32];
                        ttmp6[3] = tmp6[0]*pJ1j[3]+tmp6[1]*pJ1j[9]+tmp6[2]*pJ1j[15]+tmp6[3]*pJ1j[21]+tmp6[4]*pJ1j[27]+tmp6[5]*pJ1j[33];
                        ttmp6[4] = tmp6[0]*pJ1j[4]+tmp6[1]*pJ1j[10]+tmp6[2]*pJ1j[16]+tmp6[3]*pJ1j[22]+tmp6[4]*pJ1j[28]+tmp6[5]*pJ1j[34];
                        ttmp6[5] = tmp6[0]*pJ1j[5]+tmp6[1]*pJ1j[11]+tmp6[2]*pJ1j[17]+tmp6[3]*pJ1j[23]+tmp6[4]*pJ1j[29]+tmp6[5]*pJ1j[35];

                        ptmp = newU+(m+m_Uj[i])*(6*6);
                        if ( m_Uj[i] >= posID2 )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];
                        }
                        if ( m_Uj[i] <= posID2 && m_Ui[i] != m_Uj[i] )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];

                        }
						////////////////////////////////////////////////////////////
						///////////////////////////////////////////////////
						//J3i^T*I(i,j)*J2j

                        ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12]+tmp1[3]*pJ2j[18]+tmp1[4]*pJ2j[24]+tmp1[5]*pJ2j[30];
                        ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13]+tmp1[3]*pJ2j[19]+tmp1[4]*pJ2j[25]+tmp1[5]*pJ2j[31];
                        ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14]+tmp1[3]*pJ2j[20]+tmp1[4]*pJ2j[26]+tmp1[5]*pJ2j[32];
                        ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15]+tmp1[3]*pJ2j[21]+tmp1[4]*pJ2j[27]+tmp1[5]*pJ2j[33];
                        ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16]+tmp1[3]*pJ2j[22]+tmp1[4]*pJ2j[28]+tmp1[5]*pJ2j[34];
                        ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17]+tmp1[3]*pJ2j[23]+tmp1[4]*pJ2j[29]+tmp1[5]*pJ2j[35];

                        ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12]+tmp2[3]*pJ2j[18]+tmp2[4]*pJ2j[24]+tmp2[5]*pJ2j[30];
                        ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13]+tmp2[3]*pJ2j[19]+tmp2[4]*pJ2j[25]+tmp2[5]*pJ2j[31];
                        ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14]+tmp2[3]*pJ2j[20]+tmp2[4]*pJ2j[26]+tmp2[5]*pJ2j[32];
                        ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15]+tmp2[3]*pJ2j[21]+tmp2[4]*pJ2j[27]+tmp2[5]*pJ2j[33];
                        ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16]+tmp2[3]*pJ2j[22]+tmp2[4]*pJ2j[28]+tmp2[5]*pJ2j[34];
                        ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17]+tmp2[3]*pJ2j[23]+tmp2[4]*pJ2j[29]+tmp2[5]*pJ2j[35];

                        ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12]+tmp3[3]*pJ2j[18]+tmp3[4]*pJ2j[24]+tmp3[5]*pJ2j[30];
                        ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13]+tmp3[3]*pJ2j[19]+tmp3[4]*pJ2j[25]+tmp3[5]*pJ2j[31];
                        ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14]+tmp3[3]*pJ2j[20]+tmp3[4]*pJ2j[26]+tmp3[5]*pJ2j[32];
                        ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15]+tmp3[3]*pJ2j[21]+tmp3[4]*pJ2j[27]+tmp3[5]*pJ2j[33];
                        ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16]+tmp3[3]*pJ2j[22]+tmp3[4]*pJ2j[28]+tmp3[5]*pJ2j[34];
                        ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17]+tmp3[3]*pJ2j[23]+tmp3[4]*pJ2j[29]+tmp3[5]*pJ2j[35];

                        ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12]+tmp4[3]*pJ2j[18]+tmp4[4]*pJ2j[24]+tmp4[5]*pJ2j[30];
                        ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13]+tmp4[3]*pJ2j[19]+tmp4[4]*pJ2j[25]+tmp4[5]*pJ2j[31];
                        ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14]+tmp4[3]*pJ2j[20]+tmp4[4]*pJ2j[26]+tmp4[5]*pJ2j[32];
                        ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15]+tmp4[3]*pJ2j[21]+tmp4[4]*pJ2j[27]+tmp4[5]*pJ2j[33];
                        ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16]+tmp4[3]*pJ2j[22]+tmp4[4]*pJ2j[28]+tmp4[5]*pJ2j[34];
                        ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17]+tmp4[3]*pJ2j[23]+tmp4[4]*pJ2j[29]+tmp4[5]*pJ2j[35];

                        ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12]+tmp5[3]*pJ2j[18]+tmp5[4]*pJ2j[24]+tmp5[5]*pJ2j[30];
                        ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13]+tmp5[3]*pJ2j[19]+tmp5[4]*pJ2j[25]+tmp5[5]*pJ2j[31];
                        ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14]+tmp5[3]*pJ2j[20]+tmp5[4]*pJ2j[26]+tmp5[5]*pJ2j[32];
                        ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15]+tmp5[3]*pJ2j[21]+tmp5[4]*pJ2j[27]+tmp5[5]*pJ2j[33];
                        ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16]+tmp5[3]*pJ2j[22]+tmp5[4]*pJ2j[28]+tmp5[5]*pJ2j[34];
                        ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17]+tmp5[3]*pJ2j[23]+tmp5[4]*pJ2j[29]+tmp5[5]*pJ2j[35];

                        ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12]+tmp6[3]*pJ2j[18]+tmp6[4]*pJ2j[24]+tmp6[5]*pJ2j[30];
                        ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13]+tmp6[3]*pJ2j[19]+tmp6[4]*pJ2j[25]+tmp6[5]*pJ2j[31];
                        ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14]+tmp6[3]*pJ2j[20]+tmp6[4]*pJ2j[26]+tmp6[5]*pJ2j[32];
                        ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15]+tmp6[3]*pJ2j[21]+tmp6[4]*pJ2j[27]+tmp6[5]*pJ2j[33];
                        ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16]+tmp6[3]*pJ2j[22]+tmp6[4]*pJ2j[28]+tmp6[5]*pJ2j[34];
                        ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17]+tmp6[3]*pJ2j[23]+tmp6[4]*pJ2j[29]+tmp6[5]*pJ2j[35];

						ptmp = newU+posID2*(6*6);
                        if ( posID2 < posID )
                        {                        
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];
                        }
                        if ( posID2 > posID && m_Ui[i] != m_Uj[i] )
                        {
                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
                        }
						//////////////////////////////////////////////////
                }

				GMap_End.nU = n_newU;

				
				//================================================================
                // Compute newW, the number of elements in newW, as well as the ID of newW

                double* ptrW, *ptrnW, *ptrPID, *ptrPID2, *ptrV, *ptrV2; 

                //double* pJ1i, *pJ1j, *pJ2i, *pJ2j, *ptmp;

                int sizenW, n_newW;// n_newU is the real number of elements in newU, it should be smaller than sizenU 
                sizenW = (m_nW+n*2);
				GMap_End.W = (double*)malloc( sizenW*6*3 * sizeof(double) );
                double*  newW = GMap_End.W; 
                memset( newW, 0, sizenW*6*3 * sizeof(double) );
				GMap_End.feature = (int*)malloc( sizenW * sizeof(int) );
                int*  m_nfeature = GMap_End.feature; 
				GMap_End.photo = (int*)malloc( sizenW * sizeof(int) );
                int*  m_nphoto = GMap_End.photo;
				GMap_End.V = (double*)malloc(n*3*3*sizeof(double));

                n_newW = 0;

                ptrW = m_W;
                ptrV = m_V;
				ptrV2 = GMap_End.V;
                ptrnW = newW;

                int PreFID;
                int j=0;
                
                for ( int i = 0; i < n; i++ ) 
                {
                        PreFID = i;
                        ptrPID = newW + n_newW*6*3;
						GMap_End.FBlock[i] = n_newW;
                        ptrV = m_V + PreFID*3*3;
						ptrV2 = GMap_End.V + PreFID*3*3;
                        m_nfeature[n_newW] = PreFID;
                        m_nphoto[n_newW] = posID;
                        n_newW += 1;
						ptrPID2 = newW + n_newW*6*3;
						m_nfeature[n_newW] = PreFID;
                        m_nphoto[n_newW] = posID2;
						n_newW += 1;


                        pJ1i = J1 + m*6*6 + PreFID*3*3;
                        pJ1j = pJ1i;
                        pJ2i = J2 + m*6*6 + PreFID*3*6;
                        pJ2j = pJ2i;
                        pJ3i = J3 + m*6*6 + PreFID*3*6;
                        pJ3j = pJ3i;

					
                                //Algorithm Line 4

                                tmp1[0] = pJ2i[0]*ptrV[0]+pJ2i[6]*ptrV[3]+pJ2i[12]*ptrV[6];
                                tmp1[1] = pJ2i[0]*ptrV[1]+pJ2i[6]*ptrV[4]+pJ2i[12]*ptrV[7];
                                tmp1[2] = pJ2i[0]*ptrV[2]+pJ2i[6]*ptrV[5]+pJ2i[12]*ptrV[8];
                                

                                tmp2[0] = pJ2i[1]*ptrV[0]+pJ2i[7]*ptrV[3]+pJ2i[13]*ptrV[6];
                                tmp2[1] = pJ2i[1]*ptrV[1]+pJ2i[7]*ptrV[4]+pJ2i[13]*ptrV[7];
                                tmp2[2] = pJ2i[1]*ptrV[2]+pJ2i[7]*ptrV[5]+pJ2i[13]*ptrV[8];
                                
                                tmp3[0] = pJ2i[2]*ptrV[0]+pJ2i[8]*ptrV[3]+pJ2i[14]*ptrV[6];
                                tmp3[1] = pJ2i[2]*ptrV[1]+pJ2i[8]*ptrV[4]+pJ2i[14]*ptrV[7];
                                tmp3[2] = pJ2i[2]*ptrV[2]+pJ2i[8]*ptrV[5]+pJ2i[14]*ptrV[8];
                                
                                tmp4[0] = pJ2i[3]*ptrV[0]+pJ2i[9]*ptrV[3]+pJ2i[15]*ptrV[6];
                                tmp4[1] = pJ2i[3]*ptrV[1]+pJ2i[9]*ptrV[4]+pJ2i[15]*ptrV[7];
                                tmp4[2] = pJ2i[3]*ptrV[2]+pJ2i[9]*ptrV[5]+pJ2i[15]*ptrV[8];
                                
                                tmp5[0] = pJ2i[4]*ptrV[0]+pJ2i[10]*ptrV[3]+pJ2i[16]*ptrV[6];
                                tmp5[1] = pJ2i[4]*ptrV[1]+pJ2i[10]*ptrV[4]+pJ2i[16]*ptrV[7];
                                tmp5[2] = pJ2i[4]*ptrV[2]+pJ2i[10]*ptrV[5]+pJ2i[16]*ptrV[8];
                                
                                tmp6[0] = pJ2i[5]*ptrV[0]+pJ2i[11]*ptrV[3]+pJ2i[17]*ptrV[6];
                                tmp6[1] = pJ2i[5]*ptrV[1]+pJ2i[11]*ptrV[4]+pJ2i[17]*ptrV[7];
                                tmp6[2] = pJ2i[5]*ptrV[2]+pJ2i[11]*ptrV[5]+pJ2i[17]*ptrV[8];
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+posID*(6*6);

                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];           

                                //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID;

								ptmp[0] += ttmp1[0];
								ptmp[1] += ttmp1[1];
								ptmp[2] += ttmp1[2];
								ptmp[3] += ttmp2[0];
								ptmp[4] += ttmp2[1];
								ptmp[5] += ttmp2[2];
								ptmp[6] += ttmp3[0];
								ptmp[7] += ttmp3[1];
								ptmp[8] += ttmp3[2];
								ptmp[9] += ttmp4[0];
								ptmp[10] += ttmp4[1];
								ptmp[11] += ttmp4[2];
								ptmp[12] += ttmp5[0];
								ptmp[13] += ttmp5[1];
								ptmp[14] += ttmp5[2];
								ptmp[15] += ttmp6[0];
								ptmp[16] += ttmp6[1];
								ptmp[17] += ttmp6[2];

								///////////////////////////////////////////////////////////
								//J2 I J3

                                ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12];
                                ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13];
                                ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14];
                                ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15];
                                ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16];
                                ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17];

                                ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12];
                                ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13];
                                ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14];
                                ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15];
                                ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16];
                                ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17];

                                ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12];
                                ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13];
                                ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14];
                                ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15];
                                ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16];
                                ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17];

                                ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12];
                                ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13];
                                ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14];
                                ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15];
                                ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16];
                                ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17];

                                ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12];
                                ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13];
                                ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14];
                                ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15];
                                ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16];
                                ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17];

                                ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12];
                                ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13];
                                ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14];
                                ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15];
                                ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16];
                                ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17];

								ptmp = newU+posID2*(6*6);

								if ( posID < posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5]; 
								}
								if( posID > posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
								}
								///////////////////////////////////////////////////////////

								//Algorithm Line 3

                                tmp1[0] = pJ1i[0]*ptrV[0]+pJ1i[3]*ptrV[3]+pJ1i[6]*ptrV[6];
                                tmp1[1] = pJ1i[0]*ptrV[1]+pJ1i[3]*ptrV[4]+pJ1i[6]*ptrV[7];
                                tmp1[2] = pJ1i[0]*ptrV[2]+pJ1i[3]*ptrV[5]+pJ1i[6]*ptrV[8];                              

                                tmp2[0] = pJ1i[1]*ptrV[0]+pJ1i[4]*ptrV[3]+pJ1i[7]*ptrV[6];
                                tmp2[1] = pJ1i[1]*ptrV[1]+pJ1i[4]*ptrV[4]+pJ1i[7]*ptrV[7];
                                tmp2[2] = pJ1i[1]*ptrV[2]+pJ1i[4]*ptrV[5]+pJ1i[7]*ptrV[8];
                                
                                tmp3[0] = pJ1i[2]*ptrV[0]+pJ1i[5]*ptrV[3]+pJ1i[8]*ptrV[6];
                                tmp3[1] = pJ1i[2]*ptrV[1]+pJ1i[5]*ptrV[4]+pJ1i[8]*ptrV[7];
                                tmp3[2] = pJ1i[2]*ptrV[2]+pJ1i[5]*ptrV[5]+pJ1i[8]*ptrV[8];

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ptmp = ptrV2;

                                        ptmp[0] = ttmp1[0];
                                        ptmp[1] = ttmp1[1];
                                        ptmp[2] = ttmp1[2];
                                        ptmp[3] = ttmp2[0];
                                        ptmp[4] = ttmp2[1];
                                        ptmp[5] = ttmp2[2];
                                        ptmp[6] = ttmp3[0];
                                        ptmp[7] = ttmp3[1];
                                        ptmp[8] = ttmp3[2];

								////////////////////////////////////////////////////////////////
								//Algorithm Line 4

                                tmp1[0] = pJ3i[0]*ptrV[0]+pJ3i[6]*ptrV[3]+pJ3i[12]*ptrV[6];
                                tmp1[1] = pJ3i[0]*ptrV[1]+pJ3i[6]*ptrV[4]+pJ3i[12]*ptrV[7];
                                tmp1[2] = pJ3i[0]*ptrV[2]+pJ3i[6]*ptrV[5]+pJ3i[12]*ptrV[8];
                                

                                tmp2[0] = pJ3i[1]*ptrV[0]+pJ3i[7]*ptrV[3]+pJ3i[13]*ptrV[6];
                                tmp2[1] = pJ3i[1]*ptrV[1]+pJ3i[7]*ptrV[4]+pJ3i[13]*ptrV[7];
                                tmp2[2] = pJ3i[1]*ptrV[2]+pJ3i[7]*ptrV[5]+pJ3i[13]*ptrV[8];
                                
                                tmp3[0] = pJ3i[2]*ptrV[0]+pJ3i[8]*ptrV[3]+pJ3i[14]*ptrV[6];
                                tmp3[1] = pJ3i[2]*ptrV[1]+pJ3i[8]*ptrV[4]+pJ3i[14]*ptrV[7];
                                tmp3[2] = pJ3i[2]*ptrV[2]+pJ3i[8]*ptrV[5]+pJ3i[14]*ptrV[8];
                                
                                tmp4[0] = pJ3i[3]*ptrV[0]+pJ3i[9]*ptrV[3]+pJ3i[15]*ptrV[6];
                                tmp4[1] = pJ3i[3]*ptrV[1]+pJ3i[9]*ptrV[4]+pJ3i[15]*ptrV[7];
                                tmp4[2] = pJ3i[3]*ptrV[2]+pJ3i[9]*ptrV[5]+pJ3i[15]*ptrV[8];
                                
                                tmp5[0] = pJ3i[4]*ptrV[0]+pJ3i[10]*ptrV[3]+pJ3i[16]*ptrV[6];
                                tmp5[1] = pJ3i[4]*ptrV[1]+pJ3i[10]*ptrV[4]+pJ3i[16]*ptrV[7];
                                tmp5[2] = pJ3i[4]*ptrV[2]+pJ3i[10]*ptrV[5]+pJ3i[16]*ptrV[8];
                                
                                tmp6[0] = pJ3i[5]*ptrV[0]+pJ3i[11]*ptrV[3]+pJ3i[17]*ptrV[6];
                                tmp6[1] = pJ3i[5]*ptrV[1]+pJ3i[11]*ptrV[4]+pJ3i[17]*ptrV[7];
                                tmp6[2] = pJ3i[5]*ptrV[2]+pJ3i[11]*ptrV[5]+pJ3i[17]*ptrV[8];
                                
                                ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12];
                                ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13];
                                ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14];
                                ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15];
                                ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16];
                                ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17];

                                ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12];
                                ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13];
                                ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14];
                                ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15];
                                ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16];
                                ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17];

                                ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12];
                                ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13];
                                ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14];
                                ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15];
                                ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16];
                                ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17];

                                ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12];
                                ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13];
                                ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14];
                                ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15];
                                ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16];
                                ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17];

                                ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12];
                                ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13];
                                ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14];
                                ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15];
                                ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16];
                                ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17];

                                ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12];
                                ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13];
                                ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14];
                                ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15];
                                ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16];
                                ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17];

                                ptmp = newU+(m+posID2)*(6*6);

                                ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5];           

                                //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID2;

								ptmp[0] += ttmp1[0];
								ptmp[1] += ttmp1[1];
								ptmp[2] += ttmp1[2];
								ptmp[3] += ttmp2[0];
								ptmp[4] += ttmp2[1];
								ptmp[5] += ttmp2[2];
								ptmp[6] += ttmp3[0];
								ptmp[7] += ttmp3[1];
								ptmp[8] += ttmp3[2];
								ptmp[9] += ttmp4[0];
								ptmp[10] += ttmp4[1];
								ptmp[11] += ttmp4[2];
								ptmp[12] += ttmp5[0];
								ptmp[13] += ttmp5[1];
								ptmp[14] += ttmp5[2];
								ptmp[15] += ttmp6[0];
								ptmp[16] += ttmp6[1];
								ptmp[17] += ttmp6[2];
								////////////////////////////////////////////////////////////////
																	

                        while (j<m_nW && m_feature[j]==PreFID)
                        {
                        ptrW = m_W + j*6*3;

                        pJ1i = J1 + m_photo[j]*6*6;
                        pJ1j = J1 + m*6*6 + m_feature[j]*3*3;
                        pJ2i = J2 + m_photo[j]*6*6;
                        pJ2j = J2 + m*6*6 + m_feature[j]*3*6;
						pJ3i = J3 + m_photo[j]*6*6;
                        pJ3j = J3 + m*6*6 + m_feature[j]*3*6;

                                //Algorithm Line 4

                                tmp1[0] = pJ2i[0]*ptrW[0]+pJ2i[6]*ptrW[3]+pJ2i[12]*ptrW[6]+pJ2i[18]*ptrW[9]+pJ2i[24]*ptrW[12]+pJ2i[30]*ptrW[15];
                                tmp1[1] = pJ2i[0]*ptrW[1]+pJ2i[6]*ptrW[4]+pJ2i[12]*ptrW[7]+pJ2i[18]*ptrW[10]+pJ2i[24]*ptrW[13]+pJ2i[30]*ptrW[16];
                                tmp1[2] = pJ2i[0]*ptrW[2]+pJ2i[6]*ptrW[5]+pJ2i[12]*ptrW[8]+pJ2i[18]*ptrW[11]+pJ2i[24]*ptrW[14]+pJ2i[30]*ptrW[17];
                                
                                tmp2[0] = pJ2i[1]*ptrW[0]+pJ2i[7]*ptrW[3]+pJ2i[13]*ptrW[6]+pJ2i[19]*ptrW[9]+pJ2i[25]*ptrW[12]+pJ2i[31]*ptrW[15];
                                tmp2[1] = pJ2i[1]*ptrW[1]+pJ2i[7]*ptrW[4]+pJ2i[13]*ptrW[7]+pJ2i[19]*ptrW[10]+pJ2i[25]*ptrW[13]+pJ2i[31]*ptrW[16];
                                tmp2[2] = pJ2i[1]*ptrW[2]+pJ2i[7]*ptrW[5]+pJ2i[13]*ptrW[8]+pJ2i[19]*ptrW[11]+pJ2i[25]*ptrW[14]+pJ2i[31]*ptrW[17];
                                
                                tmp3[0] = pJ2i[2]*ptrW[0]+pJ2i[8]*ptrW[3]+pJ2i[14]*ptrW[6]+pJ2i[20]*ptrW[9]+pJ2i[26]*ptrW[12]+pJ2i[32]*ptrW[15];
                                tmp3[1] = pJ2i[2]*ptrW[1]+pJ2i[8]*ptrW[4]+pJ2i[14]*ptrW[7]+pJ2i[20]*ptrW[10]+pJ2i[26]*ptrW[13]+pJ2i[32]*ptrW[16];
                                tmp3[2] = pJ2i[2]*ptrW[2]+pJ2i[8]*ptrW[5]+pJ2i[14]*ptrW[8]+pJ2i[20]*ptrW[11]+pJ2i[26]*ptrW[14]+pJ2i[32]*ptrW[17];
                                
                                tmp4[0] = pJ2i[3]*ptrW[0]+pJ2i[9]*ptrW[3]+pJ2i[15]*ptrW[6]+pJ2i[21]*ptrW[9]+pJ2i[27]*ptrW[12]+pJ2i[33]*ptrW[15];
                                tmp4[1] = pJ2i[3]*ptrW[1]+pJ2i[9]*ptrW[4]+pJ2i[15]*ptrW[7]+pJ2i[21]*ptrW[10]+pJ2i[27]*ptrW[13]+pJ2i[33]*ptrW[16];
                                tmp4[2] = pJ2i[3]*ptrW[2]+pJ2i[9]*ptrW[5]+pJ2i[15]*ptrW[8]+pJ2i[21]*ptrW[11]+pJ2i[27]*ptrW[14]+pJ2i[33]*ptrW[17];
                                
                                tmp5[0] = pJ2i[4]*ptrW[0]+pJ2i[10]*ptrW[3]+pJ2i[16]*ptrW[6]+pJ2i[22]*ptrW[9]+pJ2i[28]*ptrW[12]+pJ2i[34]*ptrW[15];
                                tmp5[1] = pJ2i[4]*ptrW[1]+pJ2i[10]*ptrW[4]+pJ2i[16]*ptrW[7]+pJ2i[22]*ptrW[10]+pJ2i[28]*ptrW[13]+pJ2i[34]*ptrW[16];
                                tmp5[2] = pJ2i[4]*ptrW[2]+pJ2i[10]*ptrW[5]+pJ2i[16]*ptrW[8]+pJ2i[22]*ptrW[11]+pJ2i[28]*ptrW[14]+pJ2i[34]*ptrW[17];
                                
                                tmp6[0] = pJ2i[5]*ptrW[0]+pJ2i[11]*ptrW[3]+pJ2i[17]*ptrW[6]+pJ2i[23]*ptrW[9]+pJ2i[29]*ptrW[12]+pJ2i[35]*ptrW[15];
                                tmp6[1] = pJ2i[5]*ptrW[1]+pJ2i[11]*ptrW[4]+pJ2i[17]*ptrW[7]+pJ2i[23]*ptrW[10]+pJ2i[29]*ptrW[13]+pJ2i[35]*ptrW[16];
                                tmp6[2] = pJ2i[5]*ptrW[2]+pJ2i[11]*ptrW[5]+pJ2i[17]*ptrW[8]+pJ2i[23]*ptrW[11]+pJ2i[29]*ptrW[14]+pJ2i[35]*ptrW[17];
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+posID*(6*6);

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                   

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];

								//////////////////////////////////////////////////////////////////////////////
								//J2 I J3
								ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12];
                                ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13];
                                ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14];
                                ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15];
                                ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16];
                                ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17];

                                ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12];
                                ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13];
                                ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14];
                                ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15];
                                ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16];
                                ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17];

                                ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12];
                                ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13];
                                ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14];
                                ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15];
                                ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16];
                                ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17];

                                ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12];
                                ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13];
                                ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14];
                                ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15];
                                ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16];
                                ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17];

                                ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12];
                                ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13];
                                ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14];
                                ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15];
                                ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16];
                                ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17];

                                ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12];
                                ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13];
                                ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14];
                                ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15];
                                ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16];
                                ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17];
								
								ptmp = newU+posID2*(6*6);

								if ( posID < posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5]; 
								}
								if( posID > posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
								}

								////////////////////////////////////////////////////////////////////

                                        //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID;

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp2[0];
                                        ptmp[4] += ttmp2[1];
                                        ptmp[5] += ttmp2[2];
                                        ptmp[6] += ttmp3[0];
                                        ptmp[7] += ttmp3[1];
                                        ptmp[8] += ttmp3[2];
                                        ptmp[9] += ttmp4[0];
                                        ptmp[10] += ttmp4[1];
                                        ptmp[11] += ttmp4[2];
                                        ptmp[12] += ttmp5[0];
                                        ptmp[13] += ttmp5[1];
                                        ptmp[14] += ttmp5[2];
                                        ptmp[15] += ttmp6[0];
                                        ptmp[16] += ttmp6[1];
                                        ptmp[17] += ttmp6[2];

                                //Algorithm Line 3

                                tmp1[0] = pJ1i[0]*ptrW[0]+pJ1i[6]*ptrW[3]+pJ1i[12]*ptrW[6]+pJ1i[18]*ptrW[9]+pJ1i[24]*ptrW[12]+pJ1i[30]*ptrW[15];
                                tmp1[1] = pJ1i[0]*ptrW[1]+pJ1i[6]*ptrW[4]+pJ1i[12]*ptrW[7]+pJ1i[18]*ptrW[10]+pJ1i[24]*ptrW[13]+pJ1i[30]*ptrW[16];
                                tmp1[2] = pJ1i[0]*ptrW[2]+pJ1i[6]*ptrW[5]+pJ1i[12]*ptrW[8]+pJ1i[18]*ptrW[11]+pJ1i[24]*ptrW[14]+pJ1i[30]*ptrW[17];
                                
                                tmp2[0] = pJ1i[1]*ptrW[0]+pJ1i[7]*ptrW[3]+pJ1i[13]*ptrW[6]+pJ1i[19]*ptrW[9]+pJ1i[25]*ptrW[12]+pJ1i[31]*ptrW[15];
                                tmp2[1] = pJ1i[1]*ptrW[1]+pJ1i[7]*ptrW[4]+pJ1i[13]*ptrW[7]+pJ1i[19]*ptrW[10]+pJ1i[25]*ptrW[13]+pJ1i[31]*ptrW[16];
                                tmp2[2] = pJ1i[1]*ptrW[2]+pJ1i[7]*ptrW[5]+pJ1i[13]*ptrW[8]+pJ1i[19]*ptrW[11]+pJ1i[25]*ptrW[14]+pJ1i[31]*ptrW[17];
                                
                                tmp3[0] = pJ1i[2]*ptrW[0]+pJ1i[8]*ptrW[3]+pJ1i[14]*ptrW[6]+pJ1i[20]*ptrW[9]+pJ1i[26]*ptrW[12]+pJ1i[32]*ptrW[15];
                                tmp3[1] = pJ1i[2]*ptrW[1]+pJ1i[8]*ptrW[4]+pJ1i[14]*ptrW[7]+pJ1i[20]*ptrW[10]+pJ1i[26]*ptrW[13]+pJ1i[32]*ptrW[16];
                                tmp3[2] = pJ1i[2]*ptrW[2]+pJ1i[8]*ptrW[5]+pJ1i[14]*ptrW[8]+pJ1i[20]*ptrW[11]+pJ1i[26]*ptrW[14]+pJ1i[32]*ptrW[17];
                                
                                tmp4[0] = pJ1i[3]*ptrW[0]+pJ1i[9]*ptrW[3]+pJ1i[15]*ptrW[6]+pJ1i[21]*ptrW[9]+pJ1i[27]*ptrW[12]+pJ1i[33]*ptrW[15];
                                tmp4[1] = pJ1i[3]*ptrW[1]+pJ1i[9]*ptrW[4]+pJ1i[15]*ptrW[7]+pJ1i[21]*ptrW[10]+pJ1i[27]*ptrW[13]+pJ1i[33]*ptrW[16];
                                tmp4[2] = pJ1i[3]*ptrW[2]+pJ1i[9]*ptrW[5]+pJ1i[15]*ptrW[8]+pJ1i[21]*ptrW[11]+pJ1i[27]*ptrW[14]+pJ1i[33]*ptrW[17];
                                
                                tmp5[0] = pJ1i[4]*ptrW[0]+pJ1i[10]*ptrW[3]+pJ1i[16]*ptrW[6]+pJ1i[22]*ptrW[9]+pJ1i[28]*ptrW[12]+pJ1i[34]*ptrW[15];
                                tmp5[1] = pJ1i[4]*ptrW[1]+pJ1i[10]*ptrW[4]+pJ1i[16]*ptrW[7]+pJ1i[22]*ptrW[10]+pJ1i[28]*ptrW[13]+pJ1i[34]*ptrW[16];
                                tmp5[2] = pJ1i[4]*ptrW[2]+pJ1i[10]*ptrW[5]+pJ1i[16]*ptrW[8]+pJ1i[22]*ptrW[11]+pJ1i[28]*ptrW[14]+pJ1i[34]*ptrW[17];
                                
                                tmp6[0] = pJ1i[5]*ptrW[0]+pJ1i[11]*ptrW[3]+pJ1i[17]*ptrW[6]+pJ1i[23]*ptrW[9]+pJ1i[29]*ptrW[12]+pJ1i[35]*ptrW[15];
                                tmp6[1] = pJ1i[5]*ptrW[1]+pJ1i[11]*ptrW[4]+pJ1i[17]*ptrW[7]+pJ1i[23]*ptrW[10]+pJ1i[29]*ptrW[13]+pJ1i[35]*ptrW[16];
                                tmp6[2] = pJ1i[5]*ptrW[2]+pJ1i[11]*ptrW[5]+pJ1i[17]*ptrW[8]+pJ1i[23]*ptrW[11]+pJ1i[29]*ptrW[14]+pJ1i[35]*ptrW[17];
                                
                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                if ( m_photo[j] == posID )
                                {
                                        ptmp = ptrPID;
                                }
								else if ( m_photo[j] == posID2 )
								{
                                        ptmp = ptrPID2;
                                }
                                else
                                {
                                        ptmp = newW + n_newW*6*3;
                                        m_nfeature[n_newW] = m_feature[j];
                                        m_nphoto[n_newW] = m_photo[j];
                                        n_newW += 1;
                                }

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp2[0];
                                        ptmp[4] += ttmp2[1];
                                        ptmp[5] += ttmp2[2];
                                        ptmp[6] += ttmp3[0];
                                        ptmp[7] += ttmp3[1];
                                        ptmp[8] += ttmp3[2];
                                        ptmp[9] += ttmp4[0];
                                        ptmp[10] += ttmp4[1];
                                        ptmp[11] += ttmp4[2];
                                        ptmp[12] += ttmp5[0];
                                        ptmp[13] += ttmp5[1];
                                        ptmp[14] += ttmp5[2];
                                        ptmp[15] += ttmp6[0];
                                        ptmp[16] += ttmp6[1];
                                        ptmp[17] += ttmp6[2];

                                //Algorithm Line 6
                                
                                ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

                                ptmp = newU+m_photo[j]*(6*6);
                                if ( m_photo[j] <= posID )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                    
                                }
                                if ( m_photo[j] >= posID )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];
                                }

								ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12];
                                ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13];
                                ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14];
                                ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15];
                                ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16];
                                ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17];

                                ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12];
                                ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13];
                                ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14];
                                ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15];
                                ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16];
                                ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17];

                                ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12];
                                ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13];
                                ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14];
                                ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15];
                                ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16];
                                ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17];

                                ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12];
                                ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13];
                                ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14];
                                ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15];
                                ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16];
                                ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17];

                                ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12];
                                ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13];
                                ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14];
                                ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15];
                                ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16];
                                ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17];

                                ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12];
                                ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13];
                                ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14];
                                ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15];
                                ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16];
                                ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17];

                                ptmp = newU+(m+m_photo[j])*(6*6);
                                if ( m_photo[j] <= posID2 )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                    
                                }
                                if ( m_photo[j] >= posID2 )
                                {
                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];
                                }

								/////////////////////////////////////////////////////////////////////
								 //Algorithm Line 4

                                tmp1[0] = pJ3i[0]*ptrW[0]+pJ3i[6]*ptrW[3]+pJ3i[12]*ptrW[6]+pJ3i[18]*ptrW[9]+pJ3i[24]*ptrW[12]+pJ3i[30]*ptrW[15];
                                tmp1[1] = pJ3i[0]*ptrW[1]+pJ3i[6]*ptrW[4]+pJ3i[12]*ptrW[7]+pJ3i[18]*ptrW[10]+pJ3i[24]*ptrW[13]+pJ3i[30]*ptrW[16];
                                tmp1[2] = pJ3i[0]*ptrW[2]+pJ3i[6]*ptrW[5]+pJ3i[12]*ptrW[8]+pJ3i[18]*ptrW[11]+pJ3i[24]*ptrW[14]+pJ3i[30]*ptrW[17];
                                
                                tmp2[0] = pJ3i[1]*ptrW[0]+pJ3i[7]*ptrW[3]+pJ3i[13]*ptrW[6]+pJ3i[19]*ptrW[9]+pJ3i[25]*ptrW[12]+pJ3i[31]*ptrW[15];
                                tmp2[1] = pJ3i[1]*ptrW[1]+pJ3i[7]*ptrW[4]+pJ3i[13]*ptrW[7]+pJ3i[19]*ptrW[10]+pJ3i[25]*ptrW[13]+pJ3i[31]*ptrW[16];
                                tmp2[2] = pJ3i[1]*ptrW[2]+pJ3i[7]*ptrW[5]+pJ3i[13]*ptrW[8]+pJ3i[19]*ptrW[11]+pJ3i[25]*ptrW[14]+pJ3i[31]*ptrW[17];
                                
                                tmp3[0] = pJ3i[2]*ptrW[0]+pJ3i[8]*ptrW[3]+pJ3i[14]*ptrW[6]+pJ3i[20]*ptrW[9]+pJ3i[26]*ptrW[12]+pJ3i[32]*ptrW[15];
                                tmp3[1] = pJ3i[2]*ptrW[1]+pJ3i[8]*ptrW[4]+pJ3i[14]*ptrW[7]+pJ3i[20]*ptrW[10]+pJ3i[26]*ptrW[13]+pJ3i[32]*ptrW[16];
                                tmp3[2] = pJ3i[2]*ptrW[2]+pJ3i[8]*ptrW[5]+pJ3i[14]*ptrW[8]+pJ3i[20]*ptrW[11]+pJ3i[26]*ptrW[14]+pJ3i[32]*ptrW[17];
                                
                                tmp4[0] = pJ3i[3]*ptrW[0]+pJ3i[9]*ptrW[3]+pJ3i[15]*ptrW[6]+pJ3i[21]*ptrW[9]+pJ3i[27]*ptrW[12]+pJ3i[33]*ptrW[15];
                                tmp4[1] = pJ3i[3]*ptrW[1]+pJ3i[9]*ptrW[4]+pJ3i[15]*ptrW[7]+pJ3i[21]*ptrW[10]+pJ3i[27]*ptrW[13]+pJ3i[33]*ptrW[16];
                                tmp4[2] = pJ3i[3]*ptrW[2]+pJ3i[9]*ptrW[5]+pJ3i[15]*ptrW[8]+pJ3i[21]*ptrW[11]+pJ3i[27]*ptrW[14]+pJ3i[33]*ptrW[17];
                                
                                tmp5[0] = pJ3i[4]*ptrW[0]+pJ3i[10]*ptrW[3]+pJ3i[16]*ptrW[6]+pJ3i[22]*ptrW[9]+pJ3i[28]*ptrW[12]+pJ3i[34]*ptrW[15];
                                tmp5[1] = pJ3i[4]*ptrW[1]+pJ3i[10]*ptrW[4]+pJ3i[16]*ptrW[7]+pJ3i[22]*ptrW[10]+pJ3i[28]*ptrW[13]+pJ3i[34]*ptrW[16];
                                tmp5[2] = pJ3i[4]*ptrW[2]+pJ3i[10]*ptrW[5]+pJ3i[16]*ptrW[8]+pJ3i[22]*ptrW[11]+pJ3i[28]*ptrW[14]+pJ3i[34]*ptrW[17];
                                
                                tmp6[0] = pJ3i[5]*ptrW[0]+pJ3i[11]*ptrW[3]+pJ3i[17]*ptrW[6]+pJ3i[23]*ptrW[9]+pJ3i[29]*ptrW[12]+pJ3i[35]*ptrW[15];
                                tmp6[1] = pJ3i[5]*ptrW[1]+pJ3i[11]*ptrW[4]+pJ3i[17]*ptrW[7]+pJ3i[23]*ptrW[10]+pJ3i[29]*ptrW[13]+pJ3i[35]*ptrW[16];
                                tmp6[2] = pJ3i[5]*ptrW[2]+pJ3i[11]*ptrW[5]+pJ3i[17]*ptrW[8]+pJ3i[23]*ptrW[11]+pJ3i[29]*ptrW[14]+pJ3i[35]*ptrW[17];
                                
                                ttmp1[0] = tmp1[0]*pJ3j[0]+tmp1[1]*pJ3j[6]+tmp1[2]*pJ3j[12];
                                ttmp1[1] = tmp1[0]*pJ3j[1]+tmp1[1]*pJ3j[7]+tmp1[2]*pJ3j[13];
                                ttmp1[2] = tmp1[0]*pJ3j[2]+tmp1[1]*pJ3j[8]+tmp1[2]*pJ3j[14];
                                ttmp1[3] = tmp1[0]*pJ3j[3]+tmp1[1]*pJ3j[9]+tmp1[2]*pJ3j[15];
                                ttmp1[4] = tmp1[0]*pJ3j[4]+tmp1[1]*pJ3j[10]+tmp1[2]*pJ3j[16];
                                ttmp1[5] = tmp1[0]*pJ3j[5]+tmp1[1]*pJ3j[11]+tmp1[2]*pJ3j[17];

                                ttmp2[0] = tmp2[0]*pJ3j[0]+tmp2[1]*pJ3j[6]+tmp2[2]*pJ3j[12];
                                ttmp2[1] = tmp2[0]*pJ3j[1]+tmp2[1]*pJ3j[7]+tmp2[2]*pJ3j[13];
                                ttmp2[2] = tmp2[0]*pJ3j[2]+tmp2[1]*pJ3j[8]+tmp2[2]*pJ3j[14];
                                ttmp2[3] = tmp2[0]*pJ3j[3]+tmp2[1]*pJ3j[9]+tmp2[2]*pJ3j[15];
                                ttmp2[4] = tmp2[0]*pJ3j[4]+tmp2[1]*pJ3j[10]+tmp2[2]*pJ3j[16];
                                ttmp2[5] = tmp2[0]*pJ3j[5]+tmp2[1]*pJ3j[11]+tmp2[2]*pJ3j[17];

                                ttmp3[0] = tmp3[0]*pJ3j[0]+tmp3[1]*pJ3j[6]+tmp3[2]*pJ3j[12];
                                ttmp3[1] = tmp3[0]*pJ3j[1]+tmp3[1]*pJ3j[7]+tmp3[2]*pJ3j[13];
                                ttmp3[2] = tmp3[0]*pJ3j[2]+tmp3[1]*pJ3j[8]+tmp3[2]*pJ3j[14];
                                ttmp3[3] = tmp3[0]*pJ3j[3]+tmp3[1]*pJ3j[9]+tmp3[2]*pJ3j[15];
                                ttmp3[4] = tmp3[0]*pJ3j[4]+tmp3[1]*pJ3j[10]+tmp3[2]*pJ3j[16];
                                ttmp3[5] = tmp3[0]*pJ3j[5]+tmp3[1]*pJ3j[11]+tmp3[2]*pJ3j[17];

                                ttmp4[0] = tmp4[0]*pJ3j[0]+tmp4[1]*pJ3j[6]+tmp4[2]*pJ3j[12];
                                ttmp4[1] = tmp4[0]*pJ3j[1]+tmp4[1]*pJ3j[7]+tmp4[2]*pJ3j[13];
                                ttmp4[2] = tmp4[0]*pJ3j[2]+tmp4[1]*pJ3j[8]+tmp4[2]*pJ3j[14];
                                ttmp4[3] = tmp4[0]*pJ3j[3]+tmp4[1]*pJ3j[9]+tmp4[2]*pJ3j[15];
                                ttmp4[4] = tmp4[0]*pJ3j[4]+tmp4[1]*pJ3j[10]+tmp4[2]*pJ3j[16];
                                ttmp4[5] = tmp4[0]*pJ3j[5]+tmp4[1]*pJ3j[11]+tmp4[2]*pJ3j[17];

                                ttmp5[0] = tmp5[0]*pJ3j[0]+tmp5[1]*pJ3j[6]+tmp5[2]*pJ3j[12];
                                ttmp5[1] = tmp5[0]*pJ3j[1]+tmp5[1]*pJ3j[7]+tmp5[2]*pJ3j[13];
                                ttmp5[2] = tmp5[0]*pJ3j[2]+tmp5[1]*pJ3j[8]+tmp5[2]*pJ3j[14];
                                ttmp5[3] = tmp5[0]*pJ3j[3]+tmp5[1]*pJ3j[9]+tmp5[2]*pJ3j[15];
                                ttmp5[4] = tmp5[0]*pJ3j[4]+tmp5[1]*pJ3j[10]+tmp5[2]*pJ3j[16];
                                ttmp5[5] = tmp5[0]*pJ3j[5]+tmp5[1]*pJ3j[11]+tmp5[2]*pJ3j[17];

                                ttmp6[0] = tmp6[0]*pJ3j[0]+tmp6[1]*pJ3j[6]+tmp6[2]*pJ3j[12];
                                ttmp6[1] = tmp6[0]*pJ3j[1]+tmp6[1]*pJ3j[7]+tmp6[2]*pJ3j[13];
                                ttmp6[2] = tmp6[0]*pJ3j[2]+tmp6[1]*pJ3j[8]+tmp6[2]*pJ3j[14];
                                ttmp6[3] = tmp6[0]*pJ3j[3]+tmp6[1]*pJ3j[9]+tmp6[2]*pJ3j[15];
                                ttmp6[4] = tmp6[0]*pJ3j[4]+tmp6[1]*pJ3j[10]+tmp6[2]*pJ3j[16];
                                ttmp6[5] = tmp6[0]*pJ3j[5]+tmp6[1]*pJ3j[11]+tmp6[2]*pJ3j[17];

                                ptmp = newU+(m+posID2)*(6*6);

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp1[3];
                                        ptmp[4] += ttmp1[4];
                                        ptmp[5] += ttmp1[5];
                                        ptmp[6] += ttmp2[0];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp2[2];
                                        ptmp[9] += ttmp2[3];
                                        ptmp[10] += ttmp2[4];
                                        ptmp[11] += ttmp2[5];
                                        ptmp[12] += ttmp3[0];
                                        ptmp[13] += ttmp3[1];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp3[3];
                                        ptmp[16] += ttmp3[4];
                                        ptmp[17] += ttmp3[5];
                                        ptmp[18] += ttmp4[0];
                                        ptmp[19] += ttmp4[1];
                                        ptmp[20] += ttmp4[2];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp4[4];
                                        ptmp[23] += ttmp4[5];
                                        ptmp[24] += ttmp5[0];
                                        ptmp[25] += ttmp5[1];
                                        ptmp[26] += ttmp5[2];
                                        ptmp[27] += ttmp5[3];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp5[5];
                                        ptmp[30] += ttmp6[0];
                                        ptmp[31] += ttmp6[1];
                                        ptmp[32] += ttmp6[2];
                                        ptmp[33] += ttmp6[3];
                                        ptmp[34] += ttmp6[4];
                                        ptmp[35] += ttmp6[5];                                   

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp2[0];
                                        ptmp[2] += ttmp3[0];
                                        ptmp[3] += ttmp4[0];
                                        ptmp[4] += ttmp5[0];
                                        ptmp[5] += ttmp6[0];
                                        ptmp[6] += ttmp1[1];
                                        ptmp[7] += ttmp2[1];
                                        ptmp[8] += ttmp3[1];
                                        ptmp[9] += ttmp4[1];
                                        ptmp[10] += ttmp5[1];
                                        ptmp[11] += ttmp6[1];
                                        ptmp[12] += ttmp1[2];
                                        ptmp[13] += ttmp2[2];
                                        ptmp[14] += ttmp3[2];
                                        ptmp[15] += ttmp4[2];
                                        ptmp[16] += ttmp5[2];
                                        ptmp[17] += ttmp6[2];
                                        ptmp[18] += ttmp1[3];
                                        ptmp[19] += ttmp2[3];
                                        ptmp[20] += ttmp3[3];
                                        ptmp[21] += ttmp4[3];
                                        ptmp[22] += ttmp5[3];
                                        ptmp[23] += ttmp6[3];
                                        ptmp[24] += ttmp1[4];
                                        ptmp[25] += ttmp2[4];
                                        ptmp[26] += ttmp3[4];
                                        ptmp[27] += ttmp4[4];
                                        ptmp[28] += ttmp5[4];
                                        ptmp[29] += ttmp6[4];
                                        ptmp[30] += ttmp1[5];
                                        ptmp[31] += ttmp2[5];
                                        ptmp[32] += ttmp3[5];
                                        ptmp[33] += ttmp4[5];
                                        ptmp[34] += ttmp5[5];
                                        ptmp[35] += ttmp6[5];
								
								ttmp1[0] = tmp1[0]*pJ2j[0]+tmp1[1]*pJ2j[6]+tmp1[2]*pJ2j[12];
                                ttmp1[1] = tmp1[0]*pJ2j[1]+tmp1[1]*pJ2j[7]+tmp1[2]*pJ2j[13];
                                ttmp1[2] = tmp1[0]*pJ2j[2]+tmp1[1]*pJ2j[8]+tmp1[2]*pJ2j[14];
                                ttmp1[3] = tmp1[0]*pJ2j[3]+tmp1[1]*pJ2j[9]+tmp1[2]*pJ2j[15];
                                ttmp1[4] = tmp1[0]*pJ2j[4]+tmp1[1]*pJ2j[10]+tmp1[2]*pJ2j[16];
                                ttmp1[5] = tmp1[0]*pJ2j[5]+tmp1[1]*pJ2j[11]+tmp1[2]*pJ2j[17];

                                ttmp2[0] = tmp2[0]*pJ2j[0]+tmp2[1]*pJ2j[6]+tmp2[2]*pJ2j[12];
                                ttmp2[1] = tmp2[0]*pJ2j[1]+tmp2[1]*pJ2j[7]+tmp2[2]*pJ2j[13];
                                ttmp2[2] = tmp2[0]*pJ2j[2]+tmp2[1]*pJ2j[8]+tmp2[2]*pJ2j[14];
                                ttmp2[3] = tmp2[0]*pJ2j[3]+tmp2[1]*pJ2j[9]+tmp2[2]*pJ2j[15];
                                ttmp2[4] = tmp2[0]*pJ2j[4]+tmp2[1]*pJ2j[10]+tmp2[2]*pJ2j[16];
                                ttmp2[5] = tmp2[0]*pJ2j[5]+tmp2[1]*pJ2j[11]+tmp2[2]*pJ2j[17];

                                ttmp3[0] = tmp3[0]*pJ2j[0]+tmp3[1]*pJ2j[6]+tmp3[2]*pJ2j[12];
                                ttmp3[1] = tmp3[0]*pJ2j[1]+tmp3[1]*pJ2j[7]+tmp3[2]*pJ2j[13];
                                ttmp3[2] = tmp3[0]*pJ2j[2]+tmp3[1]*pJ2j[8]+tmp3[2]*pJ2j[14];
                                ttmp3[3] = tmp3[0]*pJ2j[3]+tmp3[1]*pJ2j[9]+tmp3[2]*pJ2j[15];
                                ttmp3[4] = tmp3[0]*pJ2j[4]+tmp3[1]*pJ2j[10]+tmp3[2]*pJ2j[16];
                                ttmp3[5] = tmp3[0]*pJ2j[5]+tmp3[1]*pJ2j[11]+tmp3[2]*pJ2j[17];

                                ttmp4[0] = tmp4[0]*pJ2j[0]+tmp4[1]*pJ2j[6]+tmp4[2]*pJ2j[12];
                                ttmp4[1] = tmp4[0]*pJ2j[1]+tmp4[1]*pJ2j[7]+tmp4[2]*pJ2j[13];
                                ttmp4[2] = tmp4[0]*pJ2j[2]+tmp4[1]*pJ2j[8]+tmp4[2]*pJ2j[14];
                                ttmp4[3] = tmp4[0]*pJ2j[3]+tmp4[1]*pJ2j[9]+tmp4[2]*pJ2j[15];
                                ttmp4[4] = tmp4[0]*pJ2j[4]+tmp4[1]*pJ2j[10]+tmp4[2]*pJ2j[16];
                                ttmp4[5] = tmp4[0]*pJ2j[5]+tmp4[1]*pJ2j[11]+tmp4[2]*pJ2j[17];

                                ttmp5[0] = tmp5[0]*pJ2j[0]+tmp5[1]*pJ2j[6]+tmp5[2]*pJ2j[12];
                                ttmp5[1] = tmp5[0]*pJ2j[1]+tmp5[1]*pJ2j[7]+tmp5[2]*pJ2j[13];
                                ttmp5[2] = tmp5[0]*pJ2j[2]+tmp5[1]*pJ2j[8]+tmp5[2]*pJ2j[14];
                                ttmp5[3] = tmp5[0]*pJ2j[3]+tmp5[1]*pJ2j[9]+tmp5[2]*pJ2j[15];
                                ttmp5[4] = tmp5[0]*pJ2j[4]+tmp5[1]*pJ2j[10]+tmp5[2]*pJ2j[16];
                                ttmp5[5] = tmp5[0]*pJ2j[5]+tmp5[1]*pJ2j[11]+tmp5[2]*pJ2j[17];

                                ttmp6[0] = tmp6[0]*pJ2j[0]+tmp6[1]*pJ2j[6]+tmp6[2]*pJ2j[12];
                                ttmp6[1] = tmp6[0]*pJ2j[1]+tmp6[1]*pJ2j[7]+tmp6[2]*pJ2j[13];
                                ttmp6[2] = tmp6[0]*pJ2j[2]+tmp6[1]*pJ2j[8]+tmp6[2]*pJ2j[14];
                                ttmp6[3] = tmp6[0]*pJ2j[3]+tmp6[1]*pJ2j[9]+tmp6[2]*pJ2j[15];
                                ttmp6[4] = tmp6[0]*pJ2j[4]+tmp6[1]*pJ2j[10]+tmp6[2]*pJ2j[16];
                                ttmp6[5] = tmp6[0]*pJ2j[5]+tmp6[1]*pJ2j[11]+tmp6[2]*pJ2j[17];

								ptmp = newU+posID2*(6*6);

								if ( posID > posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp1[1];
                                ptmp[2] += ttmp1[2];
                                ptmp[3] += ttmp1[3];
                                ptmp[4] += ttmp1[4];
                                ptmp[5] += ttmp1[5];
                                ptmp[6] += ttmp2[0];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp2[2];
                                ptmp[9] += ttmp2[3];
                                ptmp[10] += ttmp2[4];
                                ptmp[11] += ttmp2[5];
                                ptmp[12] += ttmp3[0];
                                ptmp[13] += ttmp3[1];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp3[3];
                                ptmp[16] += ttmp3[4];
                                ptmp[17] += ttmp3[5];
                                ptmp[18] += ttmp4[0];
                                ptmp[19] += ttmp4[1];
                                ptmp[20] += ttmp4[2];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp4[4];
                                ptmp[23] += ttmp4[5];
                                ptmp[24] += ttmp5[0];
                                ptmp[25] += ttmp5[1];
                                ptmp[26] += ttmp5[2];
                                ptmp[27] += ttmp5[3];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp5[5];
                                ptmp[30] += ttmp6[0];
                                ptmp[31] += ttmp6[1];
                                ptmp[32] += ttmp6[2];
                                ptmp[33] += ttmp6[3];
                                ptmp[34] += ttmp6[4];
                                ptmp[35] += ttmp6[5]; 
								}
								if( posID < posID2 )
								{
								ptmp[0] += ttmp1[0];
                                ptmp[1] += ttmp2[0];
                                ptmp[2] += ttmp3[0];
                                ptmp[3] += ttmp4[0];
                                ptmp[4] += ttmp5[0];
                                ptmp[5] += ttmp6[0];
                                ptmp[6] += ttmp1[1];
                                ptmp[7] += ttmp2[1];
                                ptmp[8] += ttmp3[1];
                                ptmp[9] += ttmp4[1];
                                ptmp[10] += ttmp5[1];
                                ptmp[11] += ttmp6[1];
                                ptmp[12] += ttmp1[2];
                                ptmp[13] += ttmp2[2];
                                ptmp[14] += ttmp3[2];
                                ptmp[15] += ttmp4[2];
                                ptmp[16] += ttmp5[2];
                                ptmp[17] += ttmp6[2];
                                ptmp[18] += ttmp1[3];
                                ptmp[19] += ttmp2[3];
                                ptmp[20] += ttmp3[3];
                                ptmp[21] += ttmp4[3];
                                ptmp[22] += ttmp5[3];
                                ptmp[23] += ttmp6[3];
                                ptmp[24] += ttmp1[4];
                                ptmp[25] += ttmp2[4];
                                ptmp[26] += ttmp3[4];
                                ptmp[27] += ttmp4[4];
                                ptmp[28] += ttmp5[4];
                                ptmp[29] += ttmp6[4];
                                ptmp[30] += ttmp1[5];
                                ptmp[31] += ttmp2[5];
                                ptmp[32] += ttmp3[5];
                                ptmp[33] += ttmp4[5];
                                ptmp[34] += ttmp5[5];
                                ptmp[35] += ttmp6[5];
								}

                                        //Algorithm Line 5

                                ttmp1[0] = tmp1[0]*pJ1j[0]+tmp1[1]*pJ1j[3]+tmp1[2]*pJ1j[6];
                                ttmp1[1] = tmp1[0]*pJ1j[1]+tmp1[1]*pJ1j[4]+tmp1[2]*pJ1j[7];
                                ttmp1[2] = tmp1[0]*pJ1j[2]+tmp1[1]*pJ1j[5]+tmp1[2]*pJ1j[8];

                                ttmp2[0] = tmp2[0]*pJ1j[0]+tmp2[1]*pJ1j[3]+tmp2[2]*pJ1j[6];
                                ttmp2[1] = tmp2[0]*pJ1j[1]+tmp2[1]*pJ1j[4]+tmp2[2]*pJ1j[7];
                                ttmp2[2] = tmp2[0]*pJ1j[2]+tmp2[1]*pJ1j[5]+tmp2[2]*pJ1j[8];

                                ttmp3[0] = tmp3[0]*pJ1j[0]+tmp3[1]*pJ1j[3]+tmp3[2]*pJ1j[6];
                                ttmp3[1] = tmp3[0]*pJ1j[1]+tmp3[1]*pJ1j[4]+tmp3[2]*pJ1j[7];
                                ttmp3[2] = tmp3[0]*pJ1j[2]+tmp3[1]*pJ1j[5]+tmp3[2]*pJ1j[8];

                                ttmp4[0] = tmp4[0]*pJ1j[0]+tmp4[1]*pJ1j[3]+tmp4[2]*pJ1j[6];
                                ttmp4[1] = tmp4[0]*pJ1j[1]+tmp4[1]*pJ1j[4]+tmp4[2]*pJ1j[7];
                                ttmp4[2] = tmp4[0]*pJ1j[2]+tmp4[1]*pJ1j[5]+tmp4[2]*pJ1j[8];

                                ttmp5[0] = tmp5[0]*pJ1j[0]+tmp5[1]*pJ1j[3]+tmp5[2]*pJ1j[6];
                                ttmp5[1] = tmp5[0]*pJ1j[1]+tmp5[1]*pJ1j[4]+tmp5[2]*pJ1j[7];
                                ttmp5[2] = tmp5[0]*pJ1j[2]+tmp5[1]*pJ1j[5]+tmp5[2]*pJ1j[8];

                                ttmp6[0] = tmp6[0]*pJ1j[0]+tmp6[1]*pJ1j[3]+tmp6[2]*pJ1j[6];
                                ttmp6[1] = tmp6[0]*pJ1j[1]+tmp6[1]*pJ1j[4]+tmp6[2]*pJ1j[7];
                                ttmp6[2] = tmp6[0]*pJ1j[2]+tmp6[1]*pJ1j[5]+tmp6[2]*pJ1j[8];

                                ptmp = ptrPID2;

                                        ptmp[0] += ttmp1[0];
                                        ptmp[1] += ttmp1[1];
                                        ptmp[2] += ttmp1[2];
                                        ptmp[3] += ttmp2[0];
                                        ptmp[4] += ttmp2[1];
                                        ptmp[5] += ttmp2[2];
                                        ptmp[6] += ttmp3[0];
                                        ptmp[7] += ttmp3[1];
                                        ptmp[8] += ttmp3[2];
                                        ptmp[9] += ttmp4[0];
                                        ptmp[10] += ttmp4[1];
                                        ptmp[11] += ttmp4[2];
                                        ptmp[12] += ttmp5[0];
                                        ptmp[13] += ttmp5[1];
                                        ptmp[14] += ttmp5[2];
                                        ptmp[15] += ttmp6[0];
                                        ptmp[16] += ttmp6[1];
                                        ptmp[17] += ttmp6[2];
                                
                                
								/////////////////////////////////////////////////////////////////////


                                j++;
                             }
							 
               }
			   			   
			   GMap_End.nW = n_newW;               
			   
			   free (J1);
			   free (J2);
			   free (J3);
	}
}

void CLinearSFMImp::lmj_PF3D_Divide_ConquerMono( int nLocalMapCount )
{
	double t0, t1, t2, t3;
	t0 = clock();
	int L = 0;

	while( nLocalMapCount>1 )
	{
		int N2 = nLocalMapCount%2;
		nLocalMapCount = int(nLocalMapCount/2.0 + 0.5);

		int NumLM;
		for ( int i = 0; i < nLocalMapCount; i++ )
		{
			if ( i < nLocalMapCount-1 )
				NumLM = 2;
			else
			{
				if ( N2 == 0 )
					NumLM = 2;
				else
					NumLM = 1;
			}

			for ( int j = 0; j < NumLM; j++ )
			{
				printf( "Join Level %d Local Map %d\n", L, 2*i+j+1 );

				if ( j == 0 )
				{
					m_GMap = m_LMset[2*i+j];
				}
				else
				{
					LocalMapInfo GMap_End;

					//t1 = clock();

					lmj_Transform_PF3DMono( GMap_End, m_LMset[2*i+j].Ref, m_LMset[2*i+j].ScaP, m_LMset[2*i+j].Fix );

					//t2 = clock();					
					//printf( "Tramsformation Time Use:  %lf  sec \n", (t2-t1)*0.001 );

					//t1 = clock();
					
					free( m_GMap.U );			m_GMap.U = NULL;
					free( m_GMap.V );			m_GMap.V = NULL;
					if(m_GMap.nW > 0) free( m_GMap.W );			
					free( m_GMap.Ui );			m_GMap.Ui= NULL;
					free( m_GMap.Uj );			m_GMap.Uj= NULL;
					if(m_GMap.nW > 0) free( m_GMap.photo );		
					if(m_GMap.nW > 0) free( m_GMap.feature );	
					free( m_GMap.stno );		m_GMap.stno = NULL;
					free( m_GMap.stVal );		m_GMap.stVal = NULL;
					free( m_GMap.FBlock);		m_GMap.FBlock = NULL;
					
					//t2 = clock();					
					//printf( "Free Time Use:  %lf  sec \n", (t2-t1)*0.001 );

					lmj_LinearLS_PF3DMono(GMap_End,m_LMset[2*i+j]);
				}
			}

			printf( "Generate Level %d Local Map %d\n\n", L+1, i+1 );

			if ( (i+1)%2 == 0 )
			{
				//t1 = clock();

				if ( m_GMap.Ref > m_GMap.FRef )
				{
					LocalMapInfo GMapTmp;
					lmj_Transform_PF3DMono( GMapTmp, m_GMap.FRef, m_GMap.FScaP, m_GMap.FFix );

					//free memory
					free( m_GMap.U );			m_GMap.U = NULL;
					free( m_GMap.V );			m_GMap.V = NULL;
					if(m_GMap.nW > 0) free( m_GMap.W );			
					free( m_GMap.Ui );			m_GMap.Ui= NULL;
					free( m_GMap.Uj );			m_GMap.Uj= NULL;
					if(m_GMap.nW > 0) free( m_GMap.photo );		
					if(m_GMap.nW > 0) free( m_GMap.feature );	
					free( m_GMap.stno );		m_GMap.stno = NULL;
					free( m_GMap.stVal );		m_GMap.stVal = NULL;
					free( m_GMap.FBlock);		m_GMap.FBlock = NULL;					
					m_GMap = GMapTmp;
				}

				//t2 = clock();					
				//printf( "Transformation Time Use:  %lf  sec \n\n", (t2-t1)*0.001 );

			}

			m_LMset[i] = m_GMap;
		}
		L++;
	}


	
	//t1 = clock();					
				
	if ( m_GMap.Ref > m_GMap.FRef )
	{
		LocalMapInfo GMapTmp;
		lmj_Transform_PF3DMono( GMapTmp, m_GMap.FRef, m_GMap.FScaP, m_GMap.FFix );

		//free memory
		free( m_GMap.U );			m_GMap.U = NULL;
		free( m_GMap.V );			m_GMap.V = NULL;
		if(m_GMap.nW > 0) free( m_GMap.W );			
		free( m_GMap.Ui );			m_GMap.Ui= NULL;
		free( m_GMap.Uj );			m_GMap.Uj= NULL;
		if(m_GMap.nW > 0) free( m_GMap.photo );		
		if(m_GMap.nW > 0) free( m_GMap.feature );	
		free( m_GMap.stno );		m_GMap.stno = NULL;
		free( m_GMap.stVal );		m_GMap.stVal = NULL;
		free( m_GMap.FBlock);		m_GMap.FBlock = NULL;		
		m_GMap = GMapTmp;
	}

	//t2 = clock();					
	//printf( "Final Map Transformation Time Use:  %lf  sec \n\n", (t2-t1)*0.001 );

	t1 = clock();

	t2 = (t1-t0)*0.001;

	printf( "Total Used Time:  %lf  sec\n\n", t2 );
	
	if ( m_szSt != NULL )
		lmj_SaveStateVector( m_szSt, m_GMap.stVal, m_GMap.stno, m_GMap.m*6+m_GMap.n*3 );
	if ( m_szPose != NULL && m_szFeature!= NULL )
		lmj_SavePoses_3DPF( m_szPose, m_szFeature, m_GMap.stno, m_GMap.stVal, m_GMap.m*6+m_GMap.n*3 );

	free( m_GMap.U );			m_GMap.U = NULL;
	free( m_GMap.V );			m_GMap.V = NULL;
	if(m_GMap.nW > 0) free( m_GMap.W );			
	free( m_GMap.Ui );			m_GMap.Ui= NULL;
	free( m_GMap.Uj );			m_GMap.Uj= NULL;
	if(m_GMap.nW > 0) free( m_GMap.photo );		
	if(m_GMap.nW > 0) free( m_GMap.feature );	
	free( m_GMap.stno );		m_GMap.stno = NULL;
	free( m_GMap.stVal );		m_GMap.stVal = NULL;
	free( m_GMap.FBlock);		m_GMap.FBlock = NULL;	

	system( "pause" );
}

void CLinearSFMImp::lmj_readInformationMono( LocalMapInfo& GMap_End, char* szPath )
{
	int i, j;
	FILE* fpMap1 = fopen( szPath, "r" );
	int row, tmp2;
	double tmp3;
	fscanf( fpMap1, "%d", &row );
	GMap_End.Ref = row;
	GMap_End.FRef = GMap_End.Ref;
	fscanf( fpMap1, "%d", &row );
	GMap_End.ScaP = row;
	GMap_End.FScaP = GMap_End.ScaP;
	fscanf( fpMap1, "%d", &row );
	GMap_End.Fix = row;
	GMap_End.FFix = GMap_End.Fix;
	fscanf( fpMap1, "%d", &row );
	GMap_End.Sign = row;

	fscanf( fpMap1, "%d", &row );
	GMap_End.r = row;
	GMap_End.stno = (int*)malloc( row*sizeof(int) );
	GMap_End.stVal = (double*)malloc( row*sizeof(double) );


	int m = 0, n = 0;
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d %lf", &tmp2, &tmp3 );
		GMap_End.stno[i] = tmp2;
		GMap_End.stVal[i] = tmp3;
	}

	fscanf( fpMap1, "%d", &row );
	GMap_End.m = row;
	fscanf( fpMap1, "%d", &row );
	GMap_End.n = row;

	fscanf( fpMap1, "%d", &row );
	GMap_End.nU = row;
	GMap_End.Ui = (int*)malloc( GMap_End.nU*sizeof(int) );
	GMap_End.Uj = (int*)malloc( GMap_End.nU*sizeof(int) );
	GMap_End.U  = (double*)malloc( GMap_End.nU*6*6*sizeof(double) );
	for ( int i = 0; i < 6*6*row; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.U[i] = tmp3;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.Ui[i] = tmp2;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.Uj[i] = tmp2;
	}

	fscanf( fpMap1, "%d", &row );
	GMap_End.nW = row;
	GMap_End.feature = (int*)malloc( GMap_End.nW*sizeof(int) );
	GMap_End.photo = (int*)malloc( GMap_End.nW*sizeof(int) );
	GMap_End.W  = (double*)malloc( GMap_End.nW*6*3*sizeof(double) );
	for ( int i = 0; i < 6*3*row; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.W[i] = tmp3;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.photo[i] = tmp2;
	}
	for ( int i = 0; i < row; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.feature[i] = tmp2;
	}

	GMap_End.V = (double*)malloc( GMap_End.n*3*3*sizeof(double) );
	for ( int i = 0; i < 3*3*GMap_End.n; i++ )
	{
		fscanf( fpMap1, "%lf", &tmp3 );
		GMap_End.V[i] = tmp3;
	}

	GMap_End.FBlock = (int*)malloc( GMap_End.n*sizeof(int) );
	for ( int i = 0; i < GMap_End.n; i++ )
	{
		fscanf( fpMap1, "%d", &tmp2 );
		GMap_End.FBlock[i] = tmp2;
	}

	fclose( fpMap1 );
}

void CLinearSFMImp::lmj_solveLinearSFMMono( double* stVal, double* eb, double* ea,  double* U, double*W, double* V, int* Ui, int* Uj, int* photo, int* feature, int m, int n, int nU, int nW, int Ref, int ScaP, int Fix, int Sign, int FixBlk )
{
	int i, j, ii, jj, k, l;
	int cnp = 6, pnp = 3, Usz = 36, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[6*3], sum;

	double t0, t1, t2, t3, t4, t5, t6;

	//t0 = clock();

	char* smask = (char*)malloc( m*m*sizeof(char) );
	memset( smask, 0, m*m*sizeof(char) );
		
	int* mapPhoto = (int*)malloc( n*sizeof(int) );
	int nBase = feature[0];
	int nBaseNum = 1;
	int id = 0;	

	for ( int i = 1; i < nW; i++ )
	{
		int nCur = feature[i];
		if ( nCur == nBase )
		{
			nBaseNum++;
		}
		else
		{
			mapPhoto[id] = nBaseNum;
			id++;
			nBaseNum = 1;
			nBase = feature[i];
		}		
	}
	mapPhoto[n-1] = nBaseNum;
		
	id = 0;	
	for ( int i = 0; i < n; i++ )
	{
		int num = mapPhoto[i];
		for ( int ii = 0; ii < num; ii++ )
		{
			int P1 = photo[id+ii];
			for ( int jj = ii; jj < num; jj++ )
			{
				int P2 = photo[id+jj];
				if ( P1<P2 )
					smask[P1*m+P2] = 1;
				else
					smask[P2*m+P1] = 1;
			}

		}
		id += num;
	}
	
	for ( int i = 0; i < nU; i++ )
	{
		int a = Ui[i];
		int b = Uj[i];

		if ( a < b )
			smask[a*m+b] = 1;
		else
			smask[b*m+a] = 1;
	}

	int nuis;
	for(i=nuis=0, jj=m*m; i<jj; ++i)
		nuis+=(smask[i]!=0);

	sba_crsm Sidxij;
	sba_crsm_alloc(&Sidxij, m, m, nuis);
	for(i=k=0; i<m; ++i)
	{
		Sidxij.rowptr[i]=k;
		ii=i*m;
		for(j=0; j<m; ++j)
			if(smask[ii+j])
			{
				Sidxij.val[k]=k;
				Sidxij.colidx[k++]=j;
			}
	}
	

	Sidxij.rowptr[m]=nuis;	

	double* S = (double*)malloc( 6*6*nuis*sizeof(double) );
	double* E = (double*)malloc( 6*m*sizeof(double) );
	memset( S, 0, 6*6*nuis*sizeof(double) );
	memset( E, 0, 6*m*sizeof(double) );
	double* Vold = (double*)malloc( 3*3*n*sizeof(double) );
	memcpy( Vold, V, 3*3*n*sizeof(double));
	pba_inverseV( V, m, n );	

	for ( int i = 0; i < nU; i++ )
	{
		int a = Ui[i];
		int b = Uj[i];
		pos1 = sba_crsm_elmidx( &Sidxij, a, b );

		ptr2 = S + pos1*36;
		if ( a == b )
		{					
			ptr1 = U + i * Usz;
			for(ii=0; ii<cnp; ++ii, ptr2+=6)
			{
				ptr2[ii] += ptr1[ii*cnp+ii];//
				for(jj=ii+1; jj<cnp; ++jj)
					ptr2[jj] += ptr1[ii*cnp+jj];//
			}
		}
		else
		{	
			ptr1 = U + i * Usz;
			for(ii=0; ii<cnp; ++ii, ptr2+=6)
				for(jj=0; jj<cnp; ++jj)
					ptr2[jj] += ptr1[ii*cnp+jj];//	
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];


	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( int i = 0; i < n; i++ )
	{
		numfea = mapPhoto[i];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = feature[(pos+j)];
			nP1 = photo[(pos+j)];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + (pos+j)*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;

			//WV
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )
						sum+=ptr3[k]*ptr2[jj*pnp+k];
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = 0; k < numfea; k++ )
			{
				nP2 = photo[pos+k];

				//W(V^-1)W^T
				ptr3 = W + (pos+k)*cnp*3;
				//ptrS = S + (nP1*m*36) + nP2*cnp;

				if ( nP1 == nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}
				if( nP1 < nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}				
			}
			//-W^tb
			ptr5 = eb + nF1*ebsz;
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj];
				ptrE[ii]-= sum;
			}
			//pos++;					
		}

		pos += numfea;
	}	

	cholmod_common cS;	
	cholmod_start (&cS) ; 
	int *Ap  = (int*)malloc((m+1)*sizeof(int));
	int * Aii = (int*)malloc(nuis*sizeof(int));
	pba_constructAuxCSSGN( Ap, Aii, m, smask, Ref );

	m_sparseE = cholmod_zeros( cnp*m-7, 1, CHOLMOD_REAL, &cS);
	double* Ex = (double*)m_sparseE->x;
	int nMaxS = (nuis-m)*36+m*21;	//maximum non-zero element in S matrix 
	
	m_sparseS = cholmod_allocate_sparse(m*6-7,m*6-7,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_sparseS->p;		//column pointer
	Si = (int*)m_sparseS->i;		//row pointer

	bool init = false;
	pba_constructCSSGN( Si, Sp, Sx, S, Sidxij, init, m, smask, Ref, ScaP, Fix ); //set CSS format using S matrix
	int nEx = 0;
	for ( ii = 0; ii < cnp*m; ii++ )
	{
		if ( (ii<ScaP || ii>ScaP+5) && ii!=Fix )
		{
			Ex[nEx] = E[ii];
			nEx += 1;
		}
	}

	//bool ordering = true;
	bool ordering = false;

	pba_solveCholmodGN( Ap, Aii, false, ordering, cS, nuis, m, FixBlk);
	
	double* rx = (double*)m_sparseR->x;
	nEx = 0;
	for ( int i = 0; i < 6*m; i++ )
	{
		if ((i>=ScaP && i<=ScaP+5) || i==Fix)
		{
			stVal[i] = 0;
		}
		else
		{
			stVal[i] = rx[nEx];
			nEx += 1;
		}
	}
	

	pba_solveFeatures( W, V, ea, eb, stVal, stVal+6*m, m, n, mapPhoto, photo );

	stVal[Fix] = Sign;

	memcpy( V, Vold, 3*3*n*sizeof(double) );
	free( Vold );
	free( smask );
	free( S );
	free( E );
 	free( Ap );
 	free( Aii );
	cholmod_free_factor( &m_factorS, &m_cS );
	cholmod_free_sparse( &m_sparseS, &m_cS );
	cholmod_free_dense( &m_sparseR, &m_cS );
	cholmod_free_dense( &m_sparseE, &m_cS );
	cholmod_finish( &cS );
	sba_crsm_free(&Sidxij);
}

bool CLinearSFMImp::pba_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering, cholmod_common m_cS, int m_nS, int m, int FixBlk )
{
	int i, j;
	VectorXi scalarPermutation, blockPermutation;

	//ordering = true;

	if (!init)
	{
		if (!ordering)
		{
			m_cS.nmethods = 1;
			m_cS.method[0].ordering = CHOLMOD_AMD; //CHOLMOD_COLAMD
			m_factorS = cholmod_analyze(m_sparseS, &m_cS); // symbolic factorization
		}
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m);
			if (blockPermutation.size() < m) // double space if resizing
				blockPermutation.resize(2*m);

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m;
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);
			if (! amdStatus) 
				return false;


			// blow up the permutation to the scalar matrix
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_sparseS->ncol);
			if (scalarPermutation.size() < (int)m_sparseS->ncol)
				scalarPermutation.resize(2*m_sparseS->ncol);
			size_t scalarIdx = 0;


			for ( i = 0; i < m-1; ++i)
			{
				const int &pp = blockPermutation(i);

				int base = (pp <= FixBlk-1) ? pp*6 : pp*6-1;
				int nCols= (pp==FixBlk-1) ? 5 : 6;

				for ( j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;

			}


			assert(scalarIdx == m_sparseS->ncol);

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;
			m_factorS = cholmod_analyze_p(m_sparseS, scalarPermutation.data(), NULL, 0, &m_cS);
		}
	}

	cholmod_factorize(m_sparseS, m_factorS, &m_cS); 
	m_sparseR = cholmod_solve (CHOLMOD_A, m_factorS, m_sparseE, &m_cS) ;


	return true;
}

void CLinearSFMImp::pba_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, sba_crsm& Sidxij, bool init, int m, char* smask, int Ref, int ScaP, int Fix)
{
	int ii, jj, jjj, k;
	int pos1;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;

	Sx = (double*)m_sparseS->x;
	if ( !init)
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			if (ii==Ref)
				continue;

			for ( k = 0; k < 6; k++ )
			{
				*Sp = nZ;

				if (ii*6+k==Fix)
				continue;

				//*Sp = nZ;

				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (jj==Ref)
						continue;

					if (smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if (jj*6+jjj==Fix)
									continue;

								if (jj*6+jjj<ScaP)
									*Si++ = jj*6 + jjj;
								else if (jj*6+jjj>ScaP+5 && jj*6+jjj<Fix)
									*Si++ = jj*6 + jjj - 6;
								else if (jj*6+jjj>Fix)
									*Si++ = jj*6 + jjj - 7;

								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								if (jj*6+jjj==Fix)
									continue;

								if (jj*6+jjj<ScaP)
									*Si++ = jj*6 + jjj;
								else if (jj*6+jjj>ScaP+5 && jj*6+jjj<Fix)
									*Si++ = jj*6 + jjj - 6;
								else if (jj*6+jjj>Fix)
									*Si++ = jj*6 + jjj - 7;

								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			if (ii==Ref)
				continue;

			for ( k = 0; k < 6; k++ )
			{
				if (ii*6+k==Fix)
					continue;

				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (jj==Ref)
						continue;

					if (smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if (jj*6+jjj==Fix)
									continue;
								*Sx++ = ptr5[jjj*6+k];
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								if (jj*6+jjj==Fix)
									continue;
								*Sx++ = ptr5[jjj*6+k];
							}
						}
					}
				}
			}
		}
	}
}

void CLinearSFMImp::pba_constructAuxCSSGN( int *Ap, int *Aii, int m, char* smask, int Ref)
{
	int* Cp = Ap;
	int* Ci = Aii;
	int ii, jj;
	int nZ = 0;
	for ( ii = 0; ii < m; ii++ ) 
	{
		
		if (ii==Ref)
			continue;

		*Cp = nZ;

		for( jj=0; jj<=ii; jj++ )
		{
			if (jj==Ref)
			continue;

			if (smask[jj*m+ii]==1)
			{
				if (jj<Ref)
					*Ci++ = jj;
				else if (jj>Ref)
					*Ci++ = jj-1;

				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
}

void CLinearSFMImp::lmj_LinearLS_PF3DMono( LocalMapInfo& GMap_End, LocalMapInfo& GMap_Cur )
{

	double t0, t1, t2;
	//t0 = clock();
	
        LocalMapInfo GMap_Joint;   
		    
        int i;
        int m1, m2, m, n1, n2, n;
        m1 = GMap_End.m;
        m2 = GMap_Cur.m;
        n1 = GMap_End.n;
        n2 = GMap_Cur.n;
        int* p2 = GMap_Cur.stno + m2*6;
        int* iter;
        int pos;
		int pos1, pos2, posID1, posID2;
        int* comM1 = (int*)malloc( n1*sizeof(int) );
        int* comM2 = (int*)malloc( n2*sizeof(int) );
        int FID, ncom;
        ncom = 0;
        int*  Curfeature2 = (int*)malloc( GMap_Cur.n * sizeof(int) ); 
        memset( Curfeature2, -1, GMap_Cur.n * sizeof(int) );   
		int*  CurPose2 = (int*)malloc( GMap_Cur.m * sizeof(int) ); 
        memset( CurPose2, -1, GMap_Cur.m * sizeof(int) );  

		iter = find( GMap_End.stno, GMap_End.stno+GMap_End.r, -GMap_End.Ref );
        pos1 = iter - GMap_End.stno;
		posID1 = pos1/6;
		iter = find( GMap_End.stno, GMap_End.stno+GMap_End.r, -GMap_End.ScaP );
        pos2 = iter - GMap_End.stno;
		posID2 = pos2/6;

		int*  findTmp = (int*)malloc( GMap_Cur.n*sizeof(int) );
		for ( int i = 0; i < GMap_Cur.n; i++ )
		{
			findTmp[i] = p2[i*3];
		}
	   
        for ( int i = 0; i < n1; i++ )
        {
                int stno1 = GMap_End.stno[6*m1+i*3];
                //iter = find( p2, p2+n2*3, stno1 );
				iter = find( findTmp, findTmp+n2, stno1 );
                //FID = iter - p2;
				FID = iter -findTmp;

                //if ( FID < n2*3 )
				if ( FID < n2 )
                {
                        comM1[ncom] = i;
                        //comM2[ncom] = FID/3;
						comM2[ncom] = FID;
                        //Curfeature2[FID/3] = i;
						Curfeature2[FID] = i;
                        ncom++;
                }
        }

		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Find Common Features Time Use:  %lf  sec \n", t2 );
		//t0 = clock();

        GMap_Joint.n = n1 + n2 - ncom;
        GMap_Joint.m = m1 + m2 -2;
        m = GMap_Joint.m;
        n = GMap_Joint.n;
        GMap_Joint.nU = GMap_End.nU + GMap_Cur.nU;
        GMap_Joint.nW = GMap_End.nW + GMap_Cur.nW;
        GMap_Joint.U = (double*)malloc( GMap_Joint.nU*6*6*sizeof(double) );
        GMap_Joint.W = (double*)malloc( GMap_Joint.nW*6*3*sizeof(double) );
        GMap_Joint.V = (double*)malloc( GMap_Joint.n *3*3*sizeof(double) );
        GMap_Joint.Ui= (int*)malloc( GMap_Joint.nU*sizeof(int) );
        GMap_Joint.Uj= (int*)malloc( GMap_Joint.nU*sizeof(int) );
        GMap_Joint.feature = (int*)malloc( GMap_Joint.nW*sizeof(int) );
        GMap_Joint.photo = (int*)malloc( GMap_Joint.nW*sizeof(int) );
        GMap_Joint.stno = (int*)malloc( (6*m+n*3)*sizeof(int) );
        GMap_Joint.stVal = (double*)malloc( (6*m+n*3)*sizeof(double) );
		GMap_Joint.FBlock = (int*)malloc( GMap_Joint.n * sizeof(int) );
        memset( GMap_Joint.FBlock, -1, GMap_Cur.n * sizeof(int) );

		GMap_Joint.Ref = GMap_Cur.Ref;
		GMap_Joint.r = 6*m + 3*n;
		GMap_Joint.ScaP = GMap_Cur.ScaP;
		GMap_Joint.Fix = GMap_Cur.Fix;
		GMap_Joint.Sign = GMap_Cur.Sign;

		GMap_Joint.FRef = GMap_End.FRef;
		GMap_Joint.FScaP = GMap_End.FScaP;
		GMap_Joint.FFix = GMap_End.FFix; 

		double*  eP = (double*)malloc( m*6 * sizeof(double) ); 
		memset( eP, 0, m*6 * sizeof(double) );
		double*  eF = (double*)malloc( n*3 * sizeof(double) ); 
		memset( eF, 0, n*3 * sizeof(double) ); 

		memcpy( GMap_Joint.stno, GMap_End.stno, 6*m1*sizeof(int) );

		int id = 0, pos12, pos22;
        for ( int i = 0; i < m2; i++ )
        {
                if (GMap_Cur.stno[6*i]==-GMap_Cur.Ref)
				{
					CurPose2[i] = posID1;
					pos12 = 6*i;
				}
				else if (GMap_Cur.stno[6*i]==-GMap_Cur.ScaP)
				{
					CurPose2[i] = posID2;
					pos22 = 6*i;
				}
				else
                {
                    
					GMap_Joint.stno[6*m1+id*6] = GMap_Cur.stno[i*6];
					GMap_Joint.stno[6*m1+id*6+1] = GMap_Cur.stno[i*6+1];
					GMap_Joint.stno[6*m1+id*6+2] = GMap_Cur.stno[i*6+2];
					GMap_Joint.stno[6*m1+id*6+3] = GMap_Cur.stno[i*6+3];
					GMap_Joint.stno[6*m1+id*6+4] = GMap_Cur.stno[i*6+4];
					GMap_Joint.stno[6*m1+id*6+5] = GMap_Cur.stno[i*6+5];

					CurPose2[i] = m1+id;

                    id++;
                }
        }
        
        memcpy( GMap_Joint.stno+m*6, GMap_End.stno+m1*6, n1*3*sizeof(int) );
		       
        id = 0;
        for ( int i = 0; i < n2; i++ )
        {
                if (Curfeature2[i]==-1)
                {
                        GMap_Joint.stno[6*m+n1*3+id*3] = GMap_Cur.stno[6*m2+i*3];
                        GMap_Joint.stno[6*m+n1*3+id*3+1] = GMap_Cur.stno[6*m2+i*3+1];
                        GMap_Joint.stno[6*m+n1*3+id*3+2] = GMap_Cur.stno[6*m2+i*3+2];

						Curfeature2[i] = n1+id;
                        id++;
                }
        }

		// wraparound
		double *wrp1, *wrp2, Errwrp;
		int tmpwrp;
		wrp1 = GMap_End.stVal+pos2+3;
		wrp2 = GMap_Cur.stVal+pos22+3;

		for (int i=0; i<3; i++)
		{
			if ( wrp1[i] > PI )
			{
				tmpwrp = (int)(wrp1[i]/(2*PI));
				wrp1[i] -= (tmpwrp+1)*(2*PI);
			}

			if ( wrp1[i] < -PI )
			{
				tmpwrp = (int)(wrp1[i]/(2*PI));
				wrp1[i] -= (tmpwrp-1)*(2*PI);
			}

			if ( wrp2[i] > PI )
			{
				tmpwrp = (int)(wrp2[i]/(2*PI));
				wrp2[i] -= (tmpwrp+1)*(2*PI);
			}

			if ( wrp2[i] < -PI )
			{
				tmpwrp = (int)(wrp2[i]/(2*PI));
				wrp2[i] -= (tmpwrp-1)*(2*PI);
			}

			Errwrp = wrp2[i]-wrp1[i];
			if ( Errwrp>PI )
				wrp2[i] -= 2*PI;
			else if ( Errwrp<-PI )
				wrp2[i] += 2*PI;

		}



		
		int Fl = 0;
		double* FlA;

		double *ptr1, *ptr2, *ptmp, *tmp1;
		int *ptri, *ptrj;
		ptr1 = GMap_End.U;
		ptr2 = GMap_Joint.U;
		ptri = GMap_Joint.Ui;
		ptrj = GMap_Joint.Uj;
		int ul = 0;
		for (int i=0; i<GMap_End.nU; i++)
		{
		if ( GMap_End.Ui[i]!=posID1 && GMap_End.Uj[i]!=posID1 )
		{
			if ( GMap_End.Ui[i]==posID2 && GMap_End.Uj[i]==posID2 )
			{
				Fl = 1;
				FlA = ptr2;
			}

			for (int k=0; k < 36; k++)
				ptr2[k] = ptr1[k];

			ptri[0] = GMap_End.Ui[i];
			ptrj[0] = GMap_End.Uj[i];

			ptmp = eP + ptri[0]*6;
			tmp1 = GMap_End.stVal + GMap_End.Uj[i]*6;

			ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2]+ptr1[3]*tmp1[3]+ptr1[4]*tmp1[4]+ptr1[5]*tmp1[5];
			ptmp[1] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[10]*tmp1[4]+ptr1[11]*tmp1[5];
			ptmp[2] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[15]*tmp1[3]+ptr1[16]*tmp1[4]+ptr1[17]*tmp1[5];
			ptmp[3] += ptr1[18]*tmp1[0]+ptr1[19]*tmp1[1]+ptr1[20]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[22]*tmp1[4]+ptr1[23]*tmp1[5];
			ptmp[4] += ptr1[24]*tmp1[0]+ptr1[25]*tmp1[1]+ptr1[26]*tmp1[2]+ptr1[27]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[29]*tmp1[5];
			ptmp[5] += ptr1[30]*tmp1[0]+ptr1[31]*tmp1[1]+ptr1[32]*tmp1[2]+ptr1[33]*tmp1[3]+ptr1[34]*tmp1[4]+ptr1[35]*tmp1[5];

			if ( GMap_End.Ui[i] != GMap_End.Uj[i])
			{
				ptmp = eP + ptrj[0]*6;
				tmp1 = GMap_End.stVal + GMap_End.Ui[i]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[6]*tmp1[1]+ptr1[12]*tmp1[2]+ptr1[18]*tmp1[3]+ptr1[24]*tmp1[4]+ptr1[30]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[13]*tmp1[2]+ptr1[19]*tmp1[3]+ptr1[25]*tmp1[4]+ptr1[31]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[8]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[20]*tmp1[3]+ptr1[26]*tmp1[4]+ptr1[32]*tmp1[5];
				ptmp[3] += ptr1[3]*tmp1[0]+ptr1[9]*tmp1[1]+ptr1[15]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[27]*tmp1[4]+ptr1[33]*tmp1[5];
				ptmp[4] += ptr1[4]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[16]*tmp1[2]+ptr1[22]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[34]*tmp1[5];
				ptmp[5] += ptr1[5]*tmp1[0]+ptr1[11]*tmp1[1]+ptr1[17]*tmp1[2]+ptr1[23]*tmp1[3]+ptr1[29]*tmp1[4]+ptr1[35]*tmp1[5];

			}
			
			ptr2 += 36;
			ptri += 1;
			ptrj += 1;
			ul += 1;
		}
			ptr1 += 36;
		}

		ptr1 = GMap_Cur.U;
		for (int i=0; i<GMap_Cur.nU; i++)
		{
		if ( CurPose2[GMap_Cur.Ui[i]]!=posID1 && CurPose2[GMap_Cur.Uj[i]]!=posID1 )
		{
			if ( CurPose2[GMap_Cur.Ui[i]]==posID2 && CurPose2[GMap_Cur.Uj[i]]==posID2 && Fl==1)
			{
				for (int k=0; k < 36; k++)
					FlA[k] += ptr1[k];

				ptmp = eP + posID2*6;
				tmp1 = GMap_Cur.stVal + GMap_Cur.Uj[i]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2]+ptr1[3]*tmp1[3]+ptr1[4]*tmp1[4]+ptr1[5]*tmp1[5];
				ptmp[1] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[10]*tmp1[4]+ptr1[11]*tmp1[5];
				ptmp[2] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[15]*tmp1[3]+ptr1[16]*tmp1[4]+ptr1[17]*tmp1[5];
				ptmp[3] += ptr1[18]*tmp1[0]+ptr1[19]*tmp1[1]+ptr1[20]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[22]*tmp1[4]+ptr1[23]*tmp1[5];
				ptmp[4] += ptr1[24]*tmp1[0]+ptr1[25]*tmp1[1]+ptr1[26]*tmp1[2]+ptr1[27]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[29]*tmp1[5];
				ptmp[5] += ptr1[30]*tmp1[0]+ptr1[31]*tmp1[1]+ptr1[32]*tmp1[2]+ptr1[33]*tmp1[3]+ptr1[34]*tmp1[4]+ptr1[35]*tmp1[5];

			}
			else
			{
				for (int k=0; k < 36; k++)
					ptr2[k] = ptr1[k];

				ptri[0] = CurPose2[GMap_Cur.Ui[i]];
				ptrj[0] = CurPose2[GMap_Cur.Uj[i]];
						
			ptmp = eP + ptri[0]*6;
			tmp1 = GMap_Cur.stVal + GMap_Cur.Uj[i]*6;

			ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2]+ptr1[3]*tmp1[3]+ptr1[4]*tmp1[4]+ptr1[5]*tmp1[5];
			ptmp[1] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[10]*tmp1[4]+ptr1[11]*tmp1[5];
			ptmp[2] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[15]*tmp1[3]+ptr1[16]*tmp1[4]+ptr1[17]*tmp1[5];
			ptmp[3] += ptr1[18]*tmp1[0]+ptr1[19]*tmp1[1]+ptr1[20]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[22]*tmp1[4]+ptr1[23]*tmp1[5];
			ptmp[4] += ptr1[24]*tmp1[0]+ptr1[25]*tmp1[1]+ptr1[26]*tmp1[2]+ptr1[27]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[29]*tmp1[5];
			ptmp[5] += ptr1[30]*tmp1[0]+ptr1[31]*tmp1[1]+ptr1[32]*tmp1[2]+ptr1[33]*tmp1[3]+ptr1[34]*tmp1[4]+ptr1[35]*tmp1[5];

			if ( GMap_Cur.Ui[i] != GMap_Cur.Uj[i])
			{
				ptmp = eP + ptrj[0]*6;
				tmp1 = GMap_Cur.stVal + GMap_Cur.Ui[i]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[6]*tmp1[1]+ptr1[12]*tmp1[2]+ptr1[18]*tmp1[3]+ptr1[24]*tmp1[4]+ptr1[30]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[13]*tmp1[2]+ptr1[19]*tmp1[3]+ptr1[25]*tmp1[4]+ptr1[31]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[8]*tmp1[1]+ptr1[14]*tmp1[2]+ptr1[20]*tmp1[3]+ptr1[26]*tmp1[4]+ptr1[32]*tmp1[5];
				ptmp[3] += ptr1[3]*tmp1[0]+ptr1[9]*tmp1[1]+ptr1[15]*tmp1[2]+ptr1[21]*tmp1[3]+ptr1[27]*tmp1[4]+ptr1[33]*tmp1[5];
				ptmp[4] += ptr1[4]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[16]*tmp1[2]+ptr1[22]*tmp1[3]+ptr1[28]*tmp1[4]+ptr1[34]*tmp1[5];
				ptmp[5] += ptr1[5]*tmp1[0]+ptr1[11]*tmp1[1]+ptr1[17]*tmp1[2]+ptr1[23]*tmp1[3]+ptr1[29]*tmp1[4]+ptr1[35]*tmp1[5];

			}

			ptr2 += 36;
			ptri += 1;
			ptrj += 1;
			ul += 1;
			}
		}
		ptr1 += 36;
		}

		GMap_Joint.nU = ul;

		int*  m_nfeature = GMap_Joint.feature;
		int*  m_nphoto = GMap_Joint.photo;
		
		double *ptr3, *ptr4, *ptr5, *ptr6;
		ptr1 = GMap_End.W;
		ptr3 = GMap_Joint.W;
		ptr4 = GMap_End.V;
		ptr5 = GMap_Joint.V;
		int j = 0, l = 0, a = 0, b, cnt;

		for ( int i = 0; i < GMap_End.n; i++ )	
		{
			for (int k=0; k < 9; k++)
				ptr5[k] = ptr4[k];

			ptmp = eF + i*3;
			tmp1 = GMap_End.stVal + 6*GMap_End.m + i*3;

			ptmp[0] += ptr4[0]*tmp1[0]+ptr4[1]*tmp1[1]+ptr4[2]*tmp1[2];
			ptmp[1] += ptr4[3]*tmp1[0]+ptr4[4]*tmp1[1]+ptr4[5]*tmp1[2];
			ptmp[2] += ptr4[6]*tmp1[0]+ptr4[7]*tmp1[1]+ptr4[8]*tmp1[2];

			
			cnt = 0;
			Fl = 0;
			while (j<GMap_End.nW && GMap_End.feature[j]==i)
			{
			if ( GMap_End.photo[j]!=posID1 )
			{
				if ( GMap_End.photo[j]==posID2 )
				{
					Fl = 1;
					FlA = ptr3;
				}
				
				for (int k=0; k < 18; k++)
				{ptr3[k] = ptr1[k];}

				m_nfeature[l] = i;
				m_nphoto[l] = GMap_End.photo[j];

				ptmp = eP + m_nphoto[l]*6;
				tmp1 = GMap_End.stVal + 6*GMap_End.m + i*3;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2];
				ptmp[1] += ptr1[3]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[5]*tmp1[2];
				ptmp[2] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2];
				ptmp[3] += ptr1[9]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[11]*tmp1[2];
				ptmp[4] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2];
				ptmp[5] += ptr1[15]*tmp1[0]+ptr1[16]*tmp1[1]+ptr1[17]*tmp1[2];
				
				ptmp = eF + m_nfeature[l]*3;
				tmp1 = GMap_End.stVal + GMap_End.photo[j]*6;

				ptmp[0] += ptr1[0]*tmp1[0]+ptr1[3]*tmp1[1]+ptr1[6]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[12]*tmp1[4]+ptr1[15]*tmp1[5];
				ptmp[1] += ptr1[1]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[7]*tmp1[2]+ptr1[10]*tmp1[3]+ptr1[13]*tmp1[4]+ptr1[16]*tmp1[5];
				ptmp[2] += ptr1[2]*tmp1[0]+ptr1[5]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[11]*tmp1[3]+ptr1[14]*tmp1[4]+ptr1[17]*tmp1[5];				
				
				ptr3 += 18;
				l++;
				cnt += 1;
			}
			j++;
			ptr1 += 18;
			}

			if ( a<ncom && i == comM1[a] )
			{
				ptr6 = GMap_Cur.V+comM2[a]*9;
				for (int k=0; k < 9; k++)
					ptr5[k] += ptr6[k];

				ptmp = eF + i*3;
				tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + comM2[a]*3;

				ptmp[0] += ptr6[0]*tmp1[0]+ptr6[1]*tmp1[1]+ptr6[2]*tmp1[2];
				ptmp[1] += ptr6[3]*tmp1[0]+ptr6[4]*tmp1[1]+ptr6[5]*tmp1[2];
				ptmp[2] += ptr6[6]*tmp1[0]+ptr6[7]*tmp1[1]+ptr6[8]*tmp1[2];

				b = GMap_Cur.FBlock[comM2[a]];
				if ( b !=-1 )
				{
					ptr2 = GMap_Cur.W + 18*b;

					while (b<GMap_Cur.nW && GMap_Cur.feature[b]==comM2[a])
					{
					if ( CurPose2[GMap_Cur.photo[b]]!=posID1 )
					{
					if ( CurPose2[GMap_Cur.photo[b]]==posID2 && Fl==1 )
					{
						for (int k=0; k < 18; k++)
						{FlA[k] += ptr2[k];}
						
						ptmp = eP + posID2*6;
						tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + comM2[a]*3;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[1]*tmp1[1]+ptr2[2]*tmp1[2];
						ptmp[1] += ptr2[3]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[5]*tmp1[2];
						ptmp[2] += ptr2[6]*tmp1[0]+ptr2[7]*tmp1[1]+ptr2[8]*tmp1[2];
						ptmp[3] += ptr2[9]*tmp1[0]+ptr2[10]*tmp1[1]+ptr2[11]*tmp1[2];
						ptmp[4] += ptr2[12]*tmp1[0]+ptr2[13]*tmp1[1]+ptr2[14]*tmp1[2];
						ptmp[5] += ptr2[15]*tmp1[0]+ptr2[16]*tmp1[1]+ptr2[17]*tmp1[2];
						
						ptmp = eF + i*3;
						tmp1 = GMap_Cur.stVal + GMap_Cur.photo[b]*6;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[3]*tmp1[1]+ptr2[6]*tmp1[2]+ptr2[9]*tmp1[3]+ptr2[12]*tmp1[4]+ptr2[15]*tmp1[5];
						ptmp[1] += ptr2[1]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[7]*tmp1[2]+ptr2[10]*tmp1[3]+ptr2[13]*tmp1[4]+ptr2[16]*tmp1[5];
						ptmp[2] += ptr2[2]*tmp1[0]+ptr2[5]*tmp1[1]+ptr2[8]*tmp1[2]+ptr2[11]*tmp1[3]+ptr2[14]*tmp1[4]+ptr2[17]*tmp1[5];								

					}
					else
					{
						for (int k=0; k < 18; k++)
						{ptr3[k] = ptr2[k];}

						m_nfeature[l] = i;
						m_nphoto[l] = CurPose2[GMap_Cur.photo[b]];

						ptmp = eP + m_nphoto[l]*6;
						tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + comM2[a]*3;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[1]*tmp1[1]+ptr2[2]*tmp1[2];
						ptmp[1] += ptr2[3]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[5]*tmp1[2];
						ptmp[2] += ptr2[6]*tmp1[0]+ptr2[7]*tmp1[1]+ptr2[8]*tmp1[2];
						ptmp[3] += ptr2[9]*tmp1[0]+ptr2[10]*tmp1[1]+ptr2[11]*tmp1[2];
						ptmp[4] += ptr2[12]*tmp1[0]+ptr2[13]*tmp1[1]+ptr2[14]*tmp1[2];
						ptmp[5] += ptr2[15]*tmp1[0]+ptr2[16]*tmp1[1]+ptr2[17]*tmp1[2];
						
						ptmp = eF + m_nfeature[l]*3;
						tmp1 = GMap_Cur.stVal + GMap_Cur.photo[b]*6;

						ptmp[0] += ptr2[0]*tmp1[0]+ptr2[3]*tmp1[1]+ptr2[6]*tmp1[2]+ptr2[9]*tmp1[3]+ptr2[12]*tmp1[4]+ptr2[15]*tmp1[5];
						ptmp[1] += ptr2[1]*tmp1[0]+ptr2[4]*tmp1[1]+ptr2[7]*tmp1[2]+ptr2[10]*tmp1[3]+ptr2[13]*tmp1[4]+ptr2[16]*tmp1[5];
						ptmp[2] += ptr2[2]*tmp1[0]+ptr2[5]*tmp1[1]+ptr2[8]*tmp1[2]+ptr2[11]*tmp1[3]+ptr2[14]*tmp1[4]+ptr2[17]*tmp1[5];						
						
						ptr3 += 18;
						l++;						
						cnt += 1;
					}
					}
					ptr2 += 18;
					b++;
					}
				}
				a += 1;
			}
			
			ptr4 += 9;
			ptr5 += 9;
			
			if (cnt==0)
				GMap_Joint.FBlock[i] = -1;
			else
				GMap_Joint.FBlock[i] = l-cnt;
				
		}

		j = 0;
		ptr1 = GMap_Cur.W;
		ptr4 = GMap_Cur.V;
	
		for ( int i = 0; i < GMap_Cur.n; i++ )	
		{
			if (Curfeature2[i]>=GMap_End.n)
			{

				for (int k=0; k < 9; k++)
					ptr5[k] = ptr4[k];
				ptr5 += 9;

				ptmp = eF + Curfeature2[i]*3;
				tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + i*3;

				ptmp[0] += ptr4[0]*tmp1[0]+ptr4[1]*tmp1[1]+ptr4[2]*tmp1[2];
				ptmp[1] += ptr4[3]*tmp1[0]+ptr4[4]*tmp1[1]+ptr4[5]*tmp1[2];
				ptmp[2] += ptr4[6]*tmp1[0]+ptr4[7]*tmp1[1]+ptr4[8]*tmp1[2];		

				
				cnt = 0;
				while (j<GMap_Cur.nW && GMap_Cur.feature[j]==i)
				{
				if ( CurPose2[GMap_Cur.photo[j]]!=posID1 )
				{
					for (int k=0; k < 18; k++)
					{ptr3[k] = ptr1[k];}

					m_nfeature[l] = Curfeature2[i];
					m_nphoto[l] = CurPose2[GMap_Cur.photo[j]];

					ptmp = eP + m_nphoto[l]*6;
					tmp1 = GMap_Cur.stVal + 6*GMap_Cur.m + i*3;

					ptmp[0] += ptr1[0]*tmp1[0]+ptr1[1]*tmp1[1]+ptr1[2]*tmp1[2];
					ptmp[1] += ptr1[3]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[5]*tmp1[2];
					ptmp[2] += ptr1[6]*tmp1[0]+ptr1[7]*tmp1[1]+ptr1[8]*tmp1[2];
					ptmp[3] += ptr1[9]*tmp1[0]+ptr1[10]*tmp1[1]+ptr1[11]*tmp1[2];
					ptmp[4] += ptr1[12]*tmp1[0]+ptr1[13]*tmp1[1]+ptr1[14]*tmp1[2];
					ptmp[5] += ptr1[15]*tmp1[0]+ptr1[16]*tmp1[1]+ptr1[17]*tmp1[2];

					ptmp = eF + m_nfeature[l]*3;
					tmp1 = GMap_Cur.stVal + GMap_Cur.photo[j]*6;

					ptmp[0] += ptr1[0]*tmp1[0]+ptr1[3]*tmp1[1]+ptr1[6]*tmp1[2]+ptr1[9]*tmp1[3]+ptr1[12]*tmp1[4]+ptr1[15]*tmp1[5];
					ptmp[1] += ptr1[1]*tmp1[0]+ptr1[4]*tmp1[1]+ptr1[7]*tmp1[2]+ptr1[10]*tmp1[3]+ptr1[13]*tmp1[4]+ptr1[16]*tmp1[5];
					ptmp[2] += ptr1[2]*tmp1[0]+ptr1[5]*tmp1[1]+ptr1[8]*tmp1[2]+ptr1[11]*tmp1[3]+ptr1[14]*tmp1[4]+ptr1[17]*tmp1[5];
					
					ptr3 += 18;
					l++;					
					cnt += 1;
				}
				ptr1 += 18;
				j++;
				}
				if (cnt==0)
					GMap_Joint.FBlock[Curfeature2[i]] = -1;
				else
					GMap_Joint.FBlock[Curfeature2[i]] = l-cnt;
			}
			else
			{
				while (j<GMap_Cur.nW && GMap_Cur.feature[j]==i)
				{
					ptr1 += 18;
					j++;
				}
			}
			
			ptr4 += 9;	
		}

		GMap_Joint.nW = l;


		if(comM1 != NULL) free( comM1 );
		if(comM2 != NULL) free( comM2 );
		if(CurPose2 != NULL) free( CurPose2 );
		if(Curfeature2 != NULL) free( Curfeature2 );
		if(findTmp != NULL) free( findTmp );

		if(GMap_End.U  != NULL) free( GMap_End.U );
		if( GMap_End.V != NULL) free( GMap_End.V );
		if(GMap_End.nW > 0 ) free( GMap_End.W );
		if(GMap_End.Ui != NULL) free( GMap_End.Ui );
		if(GMap_End.Uj != NULL) free( GMap_End.Uj );
		if(GMap_End.nW > 0 ) free( GMap_End.photo );
		if(GMap_End.nW > 0 ) free( GMap_End.feature );
		if( GMap_End.FBlock != NULL) free( GMap_End.FBlock );
		if(GMap_End.stno != NULL) free( GMap_End.stno );
		if(GMap_End.stVal != NULL) free( GMap_End.stVal );

		if(GMap_Cur.U != NULL) free( GMap_Cur.U );
		if(GMap_Cur.V != NULL) free( GMap_Cur.V );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.W );
		if(GMap_Cur.Ui != NULL) free( GMap_Cur.Ui );
		if(GMap_Cur.Uj != NULL) free( GMap_Cur.Uj );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.photo );
		if(GMap_Cur.nW > 0 ) free( GMap_Cur.feature );
		if(GMap_Cur.FBlock != NULL) free( GMap_Cur.FBlock );
		if(GMap_Cur.stno != NULL) free( GMap_Cur.stno );
		if(GMap_Cur.stVal != NULL) free( GMap_Cur.stVal );
		

		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Join Two Maps Time Use:  %lf  sec \n", t2 );
		//t0 = clock();

		int posFix, FixBlk;
		posFix = pos2+GMap_End.Fix;
		FixBlk = posID2-1;

		lmj_solveLinearSFMMono( GMap_Joint.stVal, eF, eP,  GMap_Joint.U, GMap_Joint.W, GMap_Joint.V, GMap_Joint.Ui, GMap_Joint.Uj, GMap_Joint.photo, GMap_Joint.feature, GMap_Joint.m, GMap_Joint.n, GMap_Joint.nU, GMap_Joint.nW, posID1, pos1, posFix, GMap_End.Sign, FixBlk );

		//t1 = clock();
		//t2 = (t1-t0)*0.001;
		//printf( "Solve Time Use:  %lf  sec \n", t2 );

		free( eP );
		free( eF );
		m_GMap = GMap_Joint;

}

void CLinearSFMImp::lmj_SavePoses_3DPF( char* szPose, char* szFea, int* stno, double* st, int n )
{
	bool bSavePose = false, bSaveFea = false;
	FILE *fpPose = NULL, *fpFea = NULL;

	set<int> poseID, featureID;
	vector<double> pose, feature;
	map<int, int>  mapIndexPose, mapIndexFea;

	if (szPose)
	{
		bSavePose = true;
		fpPose    = fopen( szPose, "w" );
	}

	if (szFea)
	{
		bSaveFea  = true;
		fpFea     = fopen( szFea, "w" );
	}

	if (!szPose && !szFea)
		return;

	int nIndexPose = 0, nIndexFea = 0;
	for ( int i = 0; i < n; i++ )
	{
		if (stno[i] <= 0 )
		{
			if ( bSavePose )
			{
				mapIndexPose[-stno[i]] = nIndexPose;
				poseID.insert( -stno[i]);
				pose.push_back( st[i] );
				pose.push_back( st[i+1] );
				pose.push_back( st[i+2] );
				pose.push_back( st[i+3] );
				pose.push_back( st[i+4] );
				pose.push_back( st[i+5] );

				nIndexPose++;
			}

			i += 5;
		}
		else
		{
			if ( bSaveFea )
			{
				mapIndexFea[stno[i]] = nIndexFea;	
				featureID.insert(stno[i]);
				feature.push_back( st[i] );
				feature.push_back( st[i+1] );
				feature.push_back( st[i+2] );
				nIndexFea++;
			}	
			i += 2;
		}
	}

	int m1 = 0, n1 = 0;
	map<int, int>::iterator find1, find2;
	if (bSavePose)
	{
		for (set<int>::iterator it1 = poseID.begin(); it1 != poseID.end(); it1++ )
		{
			m1 = *it1;
			find1 = mapIndexPose.find(m1);
			n1 = find1->second;

			fprintf( fpPose, "%d  %lf  %lf  %lf %lf  %lf  %lf\n", m1, pose[n1*6], pose[n1*6+1], pose[n1*6+2],
								pose[n1*6+3], pose[n1*6+4], pose[n1*6+5]);
		}
	}

	if (bSaveFea)
	{
		for (set<int>::iterator it2 = featureID.begin(); it2 != featureID.end(); it2++ )
		{
			m1 = *it2;
			find1 = mapIndexFea.find(m1);
			n1 = find1->second;

			fprintf( fpFea, "%d  %lf  %lf %lf\n", m1, feature[n1*3], feature[n1*3+1], feature[n1*3+2] );
		}
	}

	if( fpPose )
		fclose( fpPose );
	if( fpFea)
		fclose( fpFea );
}
