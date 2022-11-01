/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/
#include "PCACCAAnalysis.h"
#include "kaimingCommon.h"

namespace KML{

float pearsonCorrelation(int i, MatrixShen *U, int j, MatrixShen *V)
{
	float *A;
	float *B;

	A=Falloc1d(U->height);
	B=Falloc1d(V->height);

	for(int k=0;k<U->height;k++)
		A[k]=U->data[k][i];

	for(int l=0;l<V->height;l++)
		B[l]=V->data[l][j];

	////calculate standard deviation of A and B////////////
	float meanA=0.0;
	float meanB=0.0;
	float tmp=0.0;

	for(int k=0;k<U->height;k++)
		tmp+=A[k];
	meanA=tmp/U->height;

	tmp=0.0;
	for(int k=0;k<V->height;k++)
		tmp+=B[k];
	meanB=tmp/V->height;

	float SDA=0.0;
	float SDB=0.0;

	tmp=0.0;
	for(int k=0;k<U->height;k++)
		tmp+=(A[k]-meanA)*(A[k]-meanA);
	SDA=sqrt(tmp/U->height);

	tmp=0.0;
	for(int k=0;k<V->height;k++)
		tmp+=(B[k]-meanB)*(B[k]-meanB);
	SDB=sqrt(tmp/V->height);

	/////calucalte cov(A, B)/////////////////////////////
	tmp=0.0;
	float covAB=0.0;
	for(int k=0;k<U->height;k++)
		tmp+=(A[k]-meanA)*(B[k]-meanB);
	covAB=tmp/U->height;

	////calculate correlation coefficient
	float corr=0.0;
	corr=covAB/(SDA*SDB);

	return(corr);
	free(A);
	free(B);
}
//////////



void principleComponentAnalysis(MatrixShen *var, float *eigenValue, MatrixShen *eigenVector)
{
	////////////////////////////////////////////////////////////////////////////////////
	//this program is to find the principle comonents of data MatrixShen var. 
	//In var, each column represents the intensity-time series of a voxel in a specified ROI. 
	//The goal is to find the principle components. The underlying mathmatics of PCAis online everywhere

	//initialize var by substract its average value
	int m=var->height; //the length of the time-series
	int n=var->width;  //the number of voxels in the ROI
	//printf("intilization data\n");

	float *averVar;
	averVar=Falloc1d(n);
	float sumTmp;
	for(int i=0; i<n; i++)
	{
		sumTmp=0.0;
		for(int j=0; j<m;j++)
			sumTmp+=var->data[j][i];
		averVar[i]=sumTmp/m;
	}

	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			var->data[i][j]=var->data[i][j]-averVar[j];
	//printf("initialization done!\n");
	free(averVar);
	//initializaiton done!

	//calculation COV MatrixShen
	//printf("calculate COV MatrixShen!\n");
	MatrixShen *covMatrixShen;
	CreateMatrix(&covMatrixShen, n, n);
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		{
			sumTmp=0.0;
			for(int k=0; k<m; k++)
				sumTmp+=var->data[k][i]*var->data[k][j];
			covMatrixShen->data[i][j]=sumTmp/(n-1);
		}

	//calculate enginVectors and enginValues
	float *eigenValueTemp;
	MatrixShen *eigenVectorTemp;
	eigenValueTemp = vectorSHEN(0, n-1);
	CreateMatrix(&eigenVectorTemp, n, n);

	Mat_Calculate_EigenVectors_EigenValues(covMatrixShen, eigenValueTemp, eigenVectorTemp, 0);
	FreeMatrix(covMatrixShen);

	/////////////////////////////////////////////////////////////////////////////////////////////
	//compotent (eigenVector) = var*eigenvectorTemp

	float sumTemp;

	for(int j=0; j<n; j++)
	{
		for(int i=0; i<m; i++)
		{
			sumTemp = 0.0;
			for(int k=0; k<n; k++)
				sumTemp += (double)var->data[i][k]*(double)eigenVectorTemp->data[k][j];
				//sumTemp += (double)P->data[i][k]*(double)var->data[k][j];
			eigenVector->data[i][j] = sumTemp;
		}
		
		///normalize eigenvectors
		sumTemp = 0.0;
		for(int i=0; i<m; i++)
			sumTemp += (double)eigenVector->data[i][j]*(double)eigenVector->data[i][j];
		sumTemp = sqrt(sumTemp);
		for(int i=0; i<m; i++)
			eigenVector->data[i][j] /= sumTemp ;

		eigenValue[j] = eigenValueTemp[j] ;

      }
	//printf("PCA done!\n");

	FreeMatrix(eigenVectorTemp);
	free(eigenValueTemp);
}

void canonicalCorrelationAnalysis(MatrixShen *varX, MatrixShen *varY, MatrixShen *corrUV)
{
	/*this programe is to calculate the canonical correlation between two sets of variables. 
	In our application, the two sets of variables are 
	the principle components of the two ROIs.*/

	char filename[200];

	int mx=varX->height; 
	int nx=varX->width;  

	int my=varY->height; 
	int ny=varY->width;  

	if(mx!=my)
	{
		printf("error: the two set of variables don't have the same number of observations!\n");
	}

	/////////////////////////////////////////calculate cov MatrixShen of varX, varY//////////////////////////////////////////////////
	////the cov MatrixShen includes: CovXX, CovYY, CovXY, CovYX
	//step1:calculate the column-average of varX and varY.
	//

	float *averX;
	averX = Falloc1d(nx);
	for(int j=0; j<nx; j++)
	{
		float tmp=0.0;
		for(int i=0; i<mx; i++)
			tmp+= varX->data[i][j];
		averX[j] = tmp/mx;
	}

	float *averY;
	averY = Falloc1d(ny);
	for(int j=0; j<ny; j++)
	{
		float tmp=0.0;
		for(int i=0; i<my; i++)
			tmp+= varY->data[i][j];
		averY[j] = tmp/my;
	}

	//step2:substract varX and varY from its column-average.
	for(int i=0;i<mx;i++)
		for(int j=0;j<nx;j++)
			varX->data[i][j]=varX->data[i][j]-averX[j];

	for(int i=0;i<my;i++)
		for(int j=0;j<ny;j++)
			varY->data[i][j]=varY->data[i][j]-averY[j];

	free(averX);
	free(averY);

	//////////step3: calculate Cov////////////
	////1)calculate Covxx
	MatrixShen *CovXX;
	CreateMatrix(&CovXX, nx, nx);
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<nx;j++)
		{
			float tmp=0.0;
			for(int k=0; k<mx; k++)
				tmp+=varX->data[k][i]*varX->data[k][j];
			CovXX->data[i][j]=tmp/mx;
		}
	}
	//printMatrixShen(CovXX);

	////2)calculate CovYY
	MatrixShen *CovYY;
	CreateMatrix(&CovYY,ny,ny);
	for(int i=0;i<ny;i++)
	{
		for(int j=0;j<ny;j++)
		{
			float tmp=0.0;
			for(int k=0; k<my; k++)
				tmp+=varY->data[k][i]*varY->data[k][j];
			CovYY->data[i][j]=tmp/my;
		}
	}
	//printMatrixShen(CovYY);

	////3)calculate CovXY
	MatrixShen *CovXY;
	CreateMatrix(&CovXY, nx,ny);
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			float tmp=0.0;
			for(int k=0;k<mx;k++)
				tmp+=varX->data[k][i]*varY->data[k][j];
			CovXY->data[i][j]=tmp/mx;
		}
	}
	//printMatrixShen(CovXY);

	////4)calculate CovYX
	MatrixShen *CovYX;
	CreateMatrix(&CovYX,ny,nx);
	for(int i=0;i<ny;i++)
	{
		for(int j=0;j<nx;j++)
		{
			float tmp=0.0;
			for(int k=0;k<my;k++)
				tmp+=varY->data[k][i]*varX->data[k][j];
			CovYX->data[i][j]=tmp/my;
		}
	}

	//printMatrixShen(CovYX);

	////step4: diagonize those Cov MatrixShen for the calculation of Cov^(-0.5), Cov^(-1)////////////////
	//1) calculate CovXX^(-0.5), CovXX^(-1)
	float *w,**u,**v;
	w = vectorSHEN(1,nx);
	u = matrixSHEN(1,nx,1,nx);
	v = matrixSHEN(1,nx,1,nx);

	for(int i=1;i<=nx;i++)
		w[i]=0;

	for(int i=1;i<=nx;i++)
		for(int j=1;j<=nx;j++)
			v[i][j]=0;

	for(int k=1; k<=nx; k++)
		for(int l=1; l<=nx; l++) 
			u[k][l] = CovXX->data[k-1][l-1];

	svdcmp(u,nx,nx,w,v); //svd decomposition

	///CovXX=PXpre*DiagX*PXpost
	MatrixShen *PXpreTemp;
	CreateMatrix(&PXpreTemp, nx, nx);
	for(int k=1; k<=nx; k++)
	{
		for(int l=1;l<=nx;l++) 
			PXpreTemp->data[l-1][k-1] = u[l][k];
	}

	//printMatrixShen(PXpreTemp);

	MatrixShen *PXpostTemp;
	CreateMatrix(&PXpostTemp, nx, nx);
	for(int k=0; k<nx; k++)
	{
		for(int l=0;l<nx;l++)
			PXpostTemp->data[l][k] = PXpreTemp->data[k][l];
	}
	//printMatrixShen(PXpostTemp);
	
	float eValueSum=0.0;
	MatrixShen *DiagXTemp;
	CreateMatrix(&DiagXTemp, nx, nx);
	for(int k=0; k<nx; k++)
		for(int l=0;l<nx;l++)
		{
			if(k==l)
			{
				DiagXTemp->data[k][l] = w[k+1];
				eValueSum += w[k+1];
			}
			else
				DiagXTemp->data[k][l] = 0;
		}
	//printMatrixShen(DiagXTemp);


	float eValueTemp = 0.0;
	int nEValue = 0;
	float corthreshold2 = 0.9;
	for(int k=0; k<nx; k++)
	{
		eValueTemp += DiagXTemp->data[k][k]; 
		if(eValueTemp>corthreshold2*eValueSum)
		{
			nEValue = k+1; 
			break;
		}
	}

	MatrixShen *PXpre;
	CreateMatrix(&PXpre, nx, nEValue);
	for(int k=0; k<PXpre->height; k++)
	{
		for(int l=0;l<PXpre->width;l++) 
			PXpre->data[k][l] = PXpreTemp->data[k][l] ;
	}

	MatrixShen *PXpost;
	CreateMatrix(&PXpost, nEValue, nx);
	for(int k=0; k<PXpost->height; k++)
	{
		for(int l=0;l<PXpost->width;l++) 
			PXpost->data[k][l] = PXpostTemp->data[k][l] ;
	}
	
	MatrixShen *DiagX;
	CreateMatrix(&DiagX, nEValue, nEValue);
	for(int k=0; k<DiagX->height; k++)
	{
		for(int l=0;l<DiagX->width;l++) 
			DiagX->data[k][l] = DiagXTemp->data[k][l] ;
	}

	//calculate DiagX^(-1)
	MatrixShen *DiagXminusone;
	CreateMatrix(&DiagXminusone, nEValue, nEValue);
	for(int k=0; k<nEValue; k++)
		for(int l=0;l<nEValue;l++)
		{
			   if(k==l)
				   DiagXminusone->data[k][l] = 1/DiagX->data[k][l];
			   else
				   DiagXminusone->data[k][l] = 0;
		}
	//printMatrixShen(DiagXminusone);

	//calculate DiagX^(-0.5)	
	MatrixShen *DiagXminushalf;
	CreateMatrix(&DiagXminushalf, nEValue, nEValue);
	for(int k=0; k<nEValue; k++)
		for(int l=0;l<nEValue;l++)
		{
			if(k==l)
				DiagXminushalf->data[k][l] = 1/sqrt(DiagX->data[k][l]);
			else
				DiagXminushalf->data[k][l] = 0;
		}

	//printMatrixShen(DiagXminushalf);

	MatrixShen *CovXXminushalf, *MatrixShentemp;
	CreateMatrix(&CovXXminushalf, nx, nx);
	CreateMatrix(&MatrixShentemp,nx, nEValue);

	Mat_A_equal_BxC(MatrixShentemp, PXpre, DiagXminushalf);
	Mat_A_equal_BxC(CovXXminushalf, MatrixShentemp, PXpost); 
	FreeMatrix(MatrixShentemp);
	//printMatrixShen(CovXXminushalf);

	//calculate CovXX^(-1)
	MatrixShen *CovXXminusone;
	CreateMatrix(&CovXXminusone, nx, nx);
	CreateMatrix(&MatrixShentemp, nx, nEValue);
	Mat_A_equal_BxC(MatrixShentemp, PXpre, DiagXminusone);
	Mat_A_equal_BxC(CovXXminusone, MatrixShentemp, PXpost); 
	//printMatrixShen(CovXXminusone);

	FreeMatrix(MatrixShentemp);
	free_vectorSHEN(w,1,nx); 
	free_matrixSHEN(u,1,nx,1,nx);
	free_matrixSHEN(v,1,nx,1,nx);
	//////////////////////////////////////end of 1)/////////////////////////////////////////

	
	//2) calculate CovYY^(-0.5)
	w = vectorSHEN(1,ny);
	u = matrixSHEN(1,ny,1,ny);
	v = matrixSHEN(1,ny,1,ny);

	for(int i=1;i<=ny;i++)
		w[i]=0;

	for(int i=1;i<=ny;i++)
		for(int j=1;j<=ny;j++)
			v[i][j]=0;

	for(int k=1; k<=ny; k++)
		for(int l=1; l<=ny; l++) 
			u[k][l] = CovYY->data[k-1][l-1];

	svdcmp(u,ny,ny,w,v); //svd decomposition


	///CovYY=PYpre*DiagY*PYpost
	MatrixShen *PYpreTemp;
	CreateMatrix(&PYpreTemp, ny, ny);
	for(int k=1; k<=ny; k++)
	{
		for(int l=1;l<=ny;l++) 
			PYpreTemp->data[l-1][k-1] = u[l][k];
	}
	//printMatrixShen(PYpreTemp);

	MatrixShen *PYpostTemp;
	CreateMatrix(&PYpostTemp, ny, ny);
	for(int k=0; k<ny; k++)
	{
		for(int l=0;l<ny;l++)
			PYpostTemp->data[l][k] = PYpreTemp->data[k][l];
	}
	//printMatrixShen(PYpostTemp);

	float eValueSumY=0.0;
	MatrixShen *DiagYTemp;
	CreateMatrix(&DiagYTemp, ny, ny);
	for(int k=0; k<ny; k++)
		for(int l=0;l<ny;l++)
		{
			if(k==l)
			{
				DiagYTemp->data[k][l] = w[k+1];
				eValueSumY += w[k+1];
			}
			else
				DiagYTemp->data[k][l] = 0;
		}
	//printMatrixShen(DiagYTemp);

	float eValueTempY = 0.0;
	int nEValueY = 0;
	for(int k=0; k<ny; k++)
	{
		eValueTempY += DiagYTemp->data[k][k]; 
		if(eValueTempY>corthreshold2*eValueSumY)
		{
			nEValueY = k+1; 
			break;
		}
	}

	MatrixShen *PYpre;
	CreateMatrix(&PYpre, ny, nEValueY);
	for(int k=0; k<PYpre->height; k++)
	{
		for(int l=0;l<PYpre->width;l++) 
			PYpre->data[k][l] = PYpreTemp->data[k][l] ;
	}

	MatrixShen *PYpost;
	CreateMatrix(&PYpost, nEValueY, ny);
	for(int k=0; k<PYpost->height; k++)
	{
		for(int l=0;l<PYpost->width;l++) 
			PYpost->data[k][l] = PYpostTemp->data[k][l] ;
	}
	
	MatrixShen *DiagY;
	CreateMatrix(&DiagY, nEValueY, nEValueY);
	for(int k=0; k<DiagY->height; k++)
	{
		for(int l=0;l<DiagY->width;l++) 
			DiagY->data[k][l] = DiagYTemp->data[k][l] ;
	}

	//calculate DiagY^(-1)
	MatrixShen *DiagYminusone;
	CreateMatrix(&DiagYminusone, nEValueY, nEValueY);
	for(int k=0; k<nEValueY; k++)
		for(int l=0;l<nEValueY;l++)
		{
			   if(k==l)
				   DiagYminusone->data[k][l] = 1/DiagY->data[k][l];
			   else
				   DiagYminusone->data[k][l] = 0;
		}
	//printMatrixShen(DiagYminusone);

	//calculate DiagY^(-0.5)	
	MatrixShen *DiagYminushalf;
	CreateMatrix(&DiagYminushalf, nEValueY, nEValueY);
	for(int k=0; k<nEValueY; k++)
		for(int l=0;l<nEValueY;l++)
		{
			if(k==l)
				DiagYminushalf->data[k][l] = 1/sqrt(DiagY->data[k][l]);
			else
				DiagYminushalf->data[k][l] = 0;
		}

	//printMatrixShen(DiagYminushalf);

	//calculate CovYY^(-0.5)
	MatrixShen *CovYYminushalf;
	CreateMatrix(&CovYYminushalf, ny, ny);
	CreateMatrix(&MatrixShentemp,ny, nEValueY);

	Mat_A_equal_BxC(MatrixShentemp, PYpre, DiagYminushalf);
	Mat_A_equal_BxC(CovYYminushalf, MatrixShentemp, PYpost); 
	FreeMatrix(MatrixShentemp);
	//printMatrixShen(CovYYminushalf);

	//calculate CovYY^(-1)
	MatrixShen *CovYYminusone;
	CreateMatrix(&CovYYminusone, ny, ny);
	CreateMatrix(&MatrixShentemp, ny, nEValueY);
	Mat_A_equal_BxC(MatrixShentemp, PYpre, DiagYminusone);
	Mat_A_equal_BxC(CovYYminusone, MatrixShentemp, PYpost); 
	//printMatrixShen(CovYYminusone);

	FreeMatrix(MatrixShentemp);
	free_vectorSHEN(w,1,ny); 
	free_matrixSHEN(u,1,ny,1,ny);
	free_matrixSHEN(v,1,ny,1,ny);
	//////////////////////////////////////end of 2)/////////////////////////////////////////
	

	//////////////////////Step5: calculate U (latent variable of X)////////////////////////
	////1) calculate  CovXX^(-0.5)*CovXY*CovYY^(-1)*CovYX*CovXX^(-0.5) ///////////////

	MatrixShen *UMatrixShenforCanonicalVars;
	MatrixShen *MatrixShentemp1, *MatrixShentemp2, *MatrixShentemp3; 

	CreateMatrix(&MatrixShentemp1,nx, ny);
	Mat_A_equal_BxC(MatrixShentemp1, CovXXminushalf, CovXY); 

	CreateMatrix(&MatrixShentemp2,nx, ny);
	Mat_A_equal_BxC(MatrixShentemp2, MatrixShentemp1, CovYYminusone); 

	CreateMatrix(&MatrixShentemp3,nx, nx);
	Mat_A_equal_BxC(MatrixShentemp3, MatrixShentemp2, CovYX);

	CreateMatrix(&UMatrixShenforCanonicalVars, nx, nx);
	Mat_A_equal_BxC(UMatrixShenforCanonicalVars, MatrixShentemp3, CovXXminushalf); 

	FreeMatrix(MatrixShentemp1);
	FreeMatrix(MatrixShentemp2);
	FreeMatrix(MatrixShentemp3);

	////2): solve the enginevalue problem of MatrixShenforCanonicalVars*enginevectors=enginevalues*enginvectors////////////
	MatrixShen *XEngineVectors;
	float  *XEngineValues;

	CreateMatrix(&XEngineVectors,nx,nx);
	XEngineValues = vectorSHEN(0, nx-1);
	Mat_Calculate_EigenVectors_EigenValues(UMatrixShenforCanonicalVars, XEngineValues, XEngineVectors, 0);

	////3): calculate U from eigenvectors: U=X*(XEnginvectors*CovXX(-0.5))
	MatrixShen *UTemp;
	CreateMatrix(&UTemp, mx, nx);
	CreateMatrix(&MatrixShentemp1, nx, nx);
	Mat_A_equal_BxC(MatrixShentemp1, CovXXminushalf, XEngineVectors);
	Mat_A_equal_BxC(UTemp, varX, MatrixShentemp1);

	//sprintf(filename, "canoniclaVariablesU.txt");
	//fprintfMatrixShen(UTemp, filename);

	FreeMatrix(MatrixShentemp1);
	FreeMatrix(XEngineVectors);
	free(XEngineValues);

	//////////////////////Step6: calculate V (latent variable of Y)////////////////////////
	////1) calculate  CovYY^(-0.5)*CovYX*CovXX^(-1)*CovXY*CovYY^(-0.5) ///////////////

	MatrixShen *VMatrixShenforCanonicalVars;

	CreateMatrix(&MatrixShentemp1,ny, nx);
	Mat_A_equal_BxC(MatrixShentemp1, CovYYminushalf, CovYX); 

	CreateMatrix(&MatrixShentemp2,ny, nx);
	Mat_A_equal_BxC(MatrixShentemp2, MatrixShentemp1, CovXXminusone); 

	CreateMatrix(&MatrixShentemp3,ny, ny);
	Mat_A_equal_BxC(MatrixShentemp3, MatrixShentemp2, CovXY); 

	CreateMatrix(&VMatrixShenforCanonicalVars, ny, ny);
	Mat_A_equal_BxC(VMatrixShenforCanonicalVars, MatrixShentemp3, CovYYminushalf); 

	FreeMatrix(MatrixShentemp1);
	FreeMatrix(MatrixShentemp2);
	FreeMatrix(MatrixShentemp3);

	////2): solve the enginevalue problem of VMatrixShenforCanonicalVars*enginevectors=enginevalues*enginvectors////////////
	MatrixShen *YEngineVectors;
	float  *YEngineValues;

	CreateMatrix(&YEngineVectors,ny,ny);
	YEngineValues = vectorSHEN(0, ny-1);
	Mat_Calculate_EigenVectors_EigenValues(VMatrixShenforCanonicalVars, YEngineValues, YEngineVectors, 0);

	////3): calculate V from eigenvectors: V=YEnginvectors*CovYY(-0.5)*Y'
	MatrixShen *VTemp;
	CreateMatrix(&VTemp, my, ny);
	CreateMatrix(&MatrixShentemp1, ny, ny);
	Mat_A_equal_BxC(MatrixShentemp1, CovYYminushalf, YEngineVectors);
	Mat_A_equal_BxC(VTemp, varY, MatrixShentemp1);

	//sprintf(filename, "canoniclaVariablesV.txt");
	//fprintfMatrixShen(VTemp, filename);

	FreeMatrix(MatrixShentemp1);
	FreeMatrix(YEngineVectors);
	free(YEngineValues);

	FreeMatrix(PXpreTemp);
	FreeMatrix(PXpostTemp);
	FreeMatrix(DiagXTemp);
	FreeMatrix(PXpre);
	FreeMatrix(PXpost);
	FreeMatrix(DiagX);
	FreeMatrix(DiagXminusone);
	FreeMatrix(DiagXminushalf);
	FreeMatrix(CovXXminushalf);
	FreeMatrix(CovXXminusone);
	FreeMatrix(PYpreTemp);
	FreeMatrix(PYpostTemp);
	FreeMatrix(DiagYTemp);
	FreeMatrix(PYpre);
	FreeMatrix(PYpost);
	FreeMatrix(DiagY);
	FreeMatrix(DiagYminusone);
	FreeMatrix(DiagYminushalf);
	FreeMatrix(CovYYminushalf);
	FreeMatrix(CovYYminusone);
	FreeMatrix(UMatrixShenforCanonicalVars);
	FreeMatrix(VMatrixShenforCanonicalVars);


	///////////////////////////setp7: calculate the pearson correlation between the column of U and V////////////////
	//corrUV=corre(U(:,i), V(:,j)
	///
	MatrixShen *U, *V;

	int mp=0, np=0;
	mp = UTemp->height;

	if(UTemp->width >= VTemp->width)
	{
		np = VTemp->width;
		CreateMatrix(&U, mp, np);
		CreateMatrix(&V, mp, np);
		for(int i=0;i<U->height;i++)
			for(int j=0;j<U->width;j++)
			{
				U->data[i][j] = UTemp->data[i][j];
				V->data[i][j] = VTemp->data[i][j];
			}
	}
	else
	{
		np = UTemp->width;
		CreateMatrix(&U, mp, np);
		CreateMatrix(&V, mp, np);
		for(int i=0;i<U->height;i++)
			for(int j=0;j<U->width;j++)
			{
				U->data[i][j] = UTemp->data[i][j];
				V->data[i][j] = VTemp->data[i][j];
			}
	}

	for(int i=0;i<np;i++)
		for(int j=0;j<np;j++)
		{
				corrUV->data[i][j]=pearsonCorrelation(i, U, j, V); 
				if(corrUV->data[i][j] > 1)
					corrUV->data[i][j] = 1; //remove numerical errors
		}

	//printf("CCA done!\n");

	FreeMatrix(U);
	FreeMatrix(V);
	FreeMatrix(UTemp);
	FreeMatrix(VTemp);

	/////////////////////////////////////////the end of CCA ////////////////////////////////////////////////////////////////////
}



}

			


