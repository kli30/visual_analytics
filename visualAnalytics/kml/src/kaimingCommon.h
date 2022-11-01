/*
 * kaimingCommon.h
 *
 *Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.
 */

#ifndef KAIMINGFILES_H_
#define KAIMINGFILES_H_

#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include "newmat.h"
#include <stdio.h>
#include <malloc.h>
#include <cmath>
#include "Vector3D.h"
#include <vector>
#include <sstream>

using namespace NEWMAT;
using namespace std;
namespace KML {

#ifndef ByteSwapAllType
#define ByteSwapAllType(x) ByteSwap((unsigned char *) &x,sizeof(x))
#endif

#ifndef kmPut
#define kmPut(x) cout<< #x <<" "<<x<<endl
#endif

#ifndef foreach
#define foreach         BOOST_FOREACH
#endif

#ifndef PAI
#define PAI 3.1415926
#endif

void OpenReadStrmAscii(fstream& strm, const char* fileName);
void OpenWriteStrmAscii(fstream& strm, const char* fileName);
void OpenReadStrmBinary(fstream& strm, const char* fileName);
void OpenWriteStrmBinary(fstream& strm, const char* fileName);
void OpenReadStrmAscii(fstream& strm, const string& fileName);
void OpenWriteStrmAscii(fstream& strm,  const string& fileName);
void OpenReadStrmBinary(fstream& strm,  const string& fileName);
void OpenWriteStrmBinary(fstream& strm,  const string& fileName);
bool IsBigEndian();
void ByteSwap(unsigned char * b, int n);
float FindMaxValueMatrix(Matrix & matrix, int& line, int& col);
float FindMinValueMatrix(Matrix & matrix, int& line, int& col);
string operator + ( string  base, int index);
string& operator += ( string&  base, int index);
void VisualizePointsUsingSphereByVTK(vector<Vector3D<float> >& roiCenters, string fileName, float radius=4, float resolution=30);
void VisualizePointsUsingVerticeByVTKWithColors(vector< KML::Vector3D< float > >& allInOnePoints, vector< KML::Vector3D< float > >& allPointColors, string fileName);
void VisualizePointsUsingSphereByVTKWithColors(vector< KML::Vector3D< float > >& allcenters, vector< KML::Vector3D< float > >& allPointColors, string fileName,float radius=4, float resolution=30);
void VisualizeLinesByVTKWithColors(vector<Vector3D<float> >& p1, vector<Vector3D<float> >& p2, vector<Vector3D<float> >& colors, int resolution=10);
void VisualizePointsUsingVerticeByVTK(vector<Vector3D<float> >& roiCenters, string fileName);
void GetOffsetFromMHD( const char* fileName, Vector3D<float>& offset);
void GetDimsFromMHD( const char* fileName, Vector3D<float>& dims);
void GetSizesFromMHD(const char* fileName, Vector3D<size_t>& sizes);
void GetSizesFromMHD(const char* fileName, Vector3D<float>& sizes);
void GetMetaInfoFromMHD(string fileName, Vector3D<float>& offset,Vector3D<float>& dims,Vector3D<float>& sizes);
void GetImgGridFromPhysicalCoord( const Vector3D<float>& theCoord, Vector3D<size_t>&  imgGrid, const Vector3D<float>& offsets, const Vector3D<float>& dims);

void ReadIntMatrix(string fileName, vector<vector<int> >& mat);
void SaveIntMatrix(string fileName, vector<vector<int> >& mat);
void ReadFloatMatrix(string fileName, vector<vector<float> >& mat);
void SaveFloatMatrix(string fileName, vector<vector<float> >& mat);
void ReadNameList(string fileName, vector<string > & allNames);
void SaveNameList(string fileName,vector<string>& allNames);
void TransposeIntMatrix(string infile,string outFile);
void TransposeFloatMatrix(string infile,string outFile);

float StatsCovariance(std::vector< float >& vector_X, std::vector< float >& vector_Y);
float StatsMean(std::vector< float >& vector_Input);
float StatsSum(std::vector< float >& vector_Input);
void StatsNormalizeHist(std::vector< float >& vector_Input); 
float StatsVariance(std::vector< float >& vector_Input);
float StatsStdev(std::vector< float >& vector_Input);
void StatsNormalizeFisher(vector<float>& vector_Input);
void StatsTTestEqualVariance(vector<float>& group1, vector<float>& group2,double& pvalueG1Bigger, double& pvalueG2Bigger);


template<typename T>
void SaveStdVectorAsLine(string fileName,vector<T>& vec) {
    fstream outStrm;
    KML::OpenWriteStrmAscii(outStrm,fileName);
    for (int idx = 0; idx < vec.size(); ++idx )
    {
        outStrm<< vec[idx]<<" ";
    } // end for::idx;
    outStrm.close();
}


template<typename T>
void SaveStdVectorAsColumn(string fileName,vector<T>& vec) {
    fstream outStrm;
    KML::OpenWriteStrmAscii(outStrm,fileName);
    for (int idx = 0; idx < vec.size(); ++idx )
    {
        outStrm<< vec[idx]<<endl;
    } // end for::idx;
    outStrm.close();
}

template<typename T>
void ReadColumnVector(string fileName, vector<T>& dst)
{
    fstream srcStrm;
    srcStrm.open(fileName.c_str(), ios::in);
    if (NULL==srcStrm)
    {
        cout<<"error reading vector: "<<fileName<<endl;
        exit(EXIT_FAILURE);
    }
    dst.clear();
    T tmpValue=0;
    string tmpLine;
    while (  ! srcStrm.eof() && getline(srcStrm,tmpLine) ) {
        if (tmpLine.size())
        {
            T tmpValue;
            stringstream mySStrm;
            mySStrm<<tmpLine;
            mySStrm>>tmpValue;
            dst.push_back(tmpValue);
        }
    }
}



template<typename T>
void AssignValue2Matrix(vector<vector<T> >& mat, T value)
{
    for (int idx1 = 0; idx1 < mat.size(); ++idx1 )
    {
        for (int idx2 = 0; idx2 < mat[idx1].size(); ++idx2 )
        {
            mat[idx1][idx2] = value;
        } // end for::idx2;
    } // end for::idx1;
}

template<typename T>
void AssignValue2Vector(vector< T  >& vec, T value)
{
    for (int idx1 = 0; idx1 < vec.size(); ++idx1 )
    {
        vec[idx1] = value;
    } // end for::idx1;
}


template<typename T>
void CheckValueRange(T value2Check, T minValue, T maxValue)
{
    if (value2Check< minValue || value2Check > maxValue)
    {
        cout<<"Invalid value, should be within range [minValue, maxValue],I will quit."<<endl;
        exit(1);
    }
}





template<typename T>
inline void findMaxValueAndIndex(const vector<T>& vec, T& value, int& index)
{
    value = vec[0];
    index=0;

    for (int i=1; i<vec.size(); ++i)
    {
        if (vec[i]>value)
        {
            value=vec[i];
            index=i;
        }
    }
}

template<typename T>
inline void findMinValueAndIndex(const vector<T>& vec, T& value, int& index)
{
    value = vec[0];
    index=0;

    for (int i=1; i<vec.size(); ++i)
    {
        if (vec[i]<value)
        {
            value=vec[i];
            index=i;
        }
    }
}

template<typename T>
inline float CoorelOfTwoCenteredSeries(const vector<T>& centeredSerie1, const vector<T>& centeredSerie2, int start, int end)
{
    float corr=0;

    double ss1=0;
    double ss2=0;

    for (int index=start; index<end; ++index)
    {
        ss1+=centeredSerie1[index]*centeredSerie1[index];
        ss2+=centeredSerie2[index]*centeredSerie2[index];
    }

    double denominator=sqrt(ss1)*sqrt(ss2)+0.000001;

    float numerator=0;
    for (int index=start; index<end; ++index)
    {
        numerator+=centeredSerie1[index]*centeredSerie2[index];
    }
    corr=numerator/denominator;
    return (corr);
}



template<typename T>
inline float AbsCoorelOfTwoCenteredSeries(const vector<T>& centeredSerie1, const vector<T>& centeredSerie2, int start, int end)
{
    return abs(CoorelOfTwoCenteredSeries(centeredSerie1,centeredSerie2,start,end));
}

template<typename T>
inline float AbsCoorelOfTwoCenteredSeries(const vector<T>& centeredSerie1, const vector<T>& centeredSerie2)
{
    return abs(CoorelOfTwoCenteredSeries(centeredSerie1,centeredSerie2,0,centeredSerie2.size()));
}

template<typename T>
ostream&  operator << (ostream& out, Vector3D<T>& src)
{
    out<<src.x<<" "<<src.y<<" "<<src.z;
    return out;
}
template<typename T>
ostream&  operator << (ostream& out, const Vector3D<T>& src)
{
    out<<src.x<<" "<<src.y<<" "<<src.z;
    return out;
}
template<typename T>
ostream&  operator << (ostream& out, vector<T>& src)
{

    for (int i=0; i< src.size(); ++i)
        out<<src[i]<<" ";
    return out;
}

template<typename T>
vector<T>& operator*(vector<T>& vec, float ratio) {
    for (int i = 0 ; i < vec.size(); ++i) {
        vec[i]=(T) ( vec[i]*ratio);
    } //end of for loop:: i
    return vec;
}

template<typename  T>
vector<T>& operator+=(vector<T>& vec1, vector<T>& vec2) {
    for (int i = 0 ; i < vec1.size(); ++i) {
        vec1[i]=vec1[i]+vec2[i];
    } //end of for loop:: i
    return vec1;
}

template<typename T>
vector<T>& operator*=(vector<T>& vec, float den) {
    for (int i = 0 ; i < vec.size(); ++i) {
        vec[i]*=den;
    } //end of for loop:: i
    return vec;
}

template<typename T>
vector<T>& operator/=(vector<T>& vec, float den) {
    for (int i = 0 ; i < vec.size(); ++i) {
        vec[i]/=den;
    } //end of for loop:: i
    return vec;
}
template<typename T>
float FirstOrderDisofTwoSeries(vector<T>& s1, vector<T>& s2) {
    if (s1.size()!=s2.size()) {
        cout<<"Error! unbalanced vector size."<<endl;
        cout<<"FirstOrderDisofTwoSeries"<<endl;
        exit(EXIT_FAILURE);
    }
    float dis=0;
    for (int index = 0; index < s2.size(); ++index) {
        dis+=abs(s1[index]-s2[index]);
    } //end of for loop:: index
    return dis;
}


template<class T>
inline float CoorelOfTwoSeries( vector<T>& centeredSerie1,  vector<T>& centeredSerie2, int start, int end)
{

    //de-mean;

    float mean1=0;
    float mean2=0;

    for (int index=start; index<end; ++index)
    {
        mean1+=centeredSerie1[index];
        mean2+=centeredSerie2[index];
    }

    mean1/=(end-start);
    mean2/=(end-start);

    for (int index=start; index<end; ++index)
    {
        centeredSerie1[index]-=mean1;
        centeredSerie2[index]-=mean2;
    }

    float corr=0;

    double ss1=0;

    double ss2=0;

    for (int index=start; index<end; ++index)
    {
        ss1+=centeredSerie1[index]*centeredSerie1[index];
        ss2+=centeredSerie2[index]*centeredSerie2[index];
    }

    double denominator=sqrt(ss1)*sqrt(ss2)+0.000001;

    float numerator=0;
    for (int index=start; index<end; ++index)
    {
        numerator+=centeredSerie1[index]*centeredSerie2[index];
    }
    corr=numerator/denominator;
    return (corr);
}

template<class T>
void ReadVectorFromOneLineFile(vector<T>& dst, const string& fileName) {
    fstream srcStrm;
    srcStrm.open(fileName.c_str(), ios::in);
    if (NULL==srcStrm)
    {
        cout<<"error reading vector: "<<fileName<<endl;
        exit(EXIT_FAILURE);
    }
    dst.clear();
    T tmpValue=0;
    while (  ! srcStrm.eof() && srcStrm>>tmpValue ) {
        dst.push_back(tmpValue);
    }
}


template<class T>
void ReadVectorFromSStream(vector<T>& dst, stringstream& fileStrm) {
    dst.clear();
    T tmpValue=0;
    while ( ! fileStrm.eof() && fileStrm>>tmpValue ) {
        dst.push_back(tmpValue);
    }
}


template<class T>
void ReadVectorFromString(vector<T>& dst, string& values) {
    stringstream fileStrm;
    fileStrm<<values;
    dst.clear();
    T tmpValue=0;
    while (  ! fileStrm.eof() && fileStrm>>tmpValue ) {
        dst.push_back(tmpValue);
    }
}

template<class T>
float L1Distance(const vector<T>& v1,const vector<T>& v2) {
    if (v1.size()!=v2.size())
    {
        cerr<<"length is not equal!"<<endl;
        exit(1);
    }

    float dis=0;
    for ( int idx=0; idx < v1.size() ; ++idx ) {
        dis+= abs(v1[idx]-v2[idx]);
    }//end of for::idx

    return dis;
}


template <typename T>
string NumberToString ( T Number )
{
    stringstream ss;
    ss << Number;
    return ss.str();
}

template<typename T>
void TransposeMatrix(vector<vector<T> >& mat, vector<vector<T> >& matT)
{
    matT.clear();
    matT.resize(mat[0].size(), vector<T>(mat.size(),T(0)));
    for (int line = 0; line < mat.size(); ++line )
    {
        for (int col = 0; col < mat[line].size(); ++col )
        {
            matT[col][line] = mat[line][col];
        } // end for::col;
    } // end for::line;
}
template<typename T>
ostream&  operator << (ostream& out, const vector<vector<T> >& src)
{
    for (int idx = 0; idx < src.size(); ++idx )
    {

        for (int idx2 = 0; idx2 < src[idx].size(); ++idx2 )
        {
            out<<src[idx][idx2]<<" ";
        } // end for::idx2;

        out<<endl;
    } // end for::idx;
    return out;
}

template<typename T>
void ReadLineVector(string fileName, vector<T>& vec)
{
    KML::ReadVectorFromOneLineFile<T>(vec,fileName);
}







///////////////////////////////////////////////////
//ding gang libs
typedef struct
{
    float x;
    float y;
    float z;
} Fvector3d;
Fvector3d ***Fvector3dalloc3d(int i_size,int j_size,int k_size);
void Fvector3dfree3d(Fvector3d ***array,int k_size,int i_size);

/* Feb 2002, for warping DTI */
typedef struct
{
    float v1 ;
    float v2 ;
    float v3 ;
    float v4 ;
    float v5 ;
    float v6 ;
} DTIattribute ;



/* Nov 2001, for Head brain image */
typedef struct
{
    unsigned char Edge;
    unsigned char Tiss;
    unsigned char Geom;
    unsigned char BGvlm;
    unsigned char CSFvlm;
    unsigned char VNvlm;
} HeadImgAttribute ;



/* June 2001, for skull-stripped brain image */
typedef struct
{
    unsigned char Edge;
    unsigned char Tiss;
    unsigned char Geom;
    unsigned char VNvlm;
    unsigned char CSFBG;
} ImgAttribute ;



typedef struct MatrixStruct
{
    double **data;
    int height, width;
} MatrixShen;

typedef struct uc_MatrixStruct
{
    unsigned char **data;
    int height, width;
} uc_Matrix;

typedef struct i_MatrixStruct
{
    int **data;
    int height, width;
} i_Matrix;

#ifndef TRUE
#define TRUE     1
#endif

#define FREE_ARG char*
#define NR_END 1

#ifndef PI
#define PI 3.141592653589793115997963468544185161590576171875
#endif

static float tempr;
#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?  (iminarg1) : (iminarg2))

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


float  **matrixSHEN(int nrl,int nrh,int ncl,int nch);

void   nrerrorSHEN(char *error_text);
double *dvectorSHEN(int nl, int nh);
void   free_dvectorSHEN(double *v, int nl, int nh);
float  *vectorSHEN(int nl, int nh);
void   free_vectorSHEN(float *v, int nl, int nh);
float  **matrixSHEN(int nrl,int nrh,int ncl,int nch);
void   free_matrixSHEN(float **m,int nrl,int nrh,int ncl,int nch);
double **dmatrixSHEN(int nrl, int nrh, int ncl, int nch);
void   free_dmatrixSHEN(double **m, int nrl, int nrh, int ncl, int nch);

//double log2(double a);
void   sort(double *Y, int *I, double *A, int length);
void   minimun(double *Y, int *I, double *A, int length);
void   Mat_Abs(MatrixShen *A);
void   Mat_Mean(double *mean, MatrixShen *A);
void   Mat_Variance(double *variance, MatrixShen *A) ;
void   Mat_Vector(MatrixShen *A, float *a);
void   Mat_Shift(MatrixShen *A, MatrixShen *B, int side);
void   Mat_Zeros(MatrixShen *A);
void   Mat_Zeros_uc(uc_Matrix *A);
void   Mat_Zeros_i(i_Matrix *A);
void   CreateMatrix(MatrixShen **M, int hei, int wid);
void   FreeMatrix(MatrixShen *M);
void   Create_i_Matrix(i_Matrix **M, int hei, int wid);
void   Free_i_Matrix(i_Matrix *M);
void   Create_uc_Matrix(uc_Matrix **M, int hei, int wid);
void   Free_uc_Matrix(uc_Matrix *M);
void   Mat_FFT2(MatrixShen *Output_real, MatrixShen *Output_imag, MatrixShen *Input_real, MatrixShen *Input_imag);
void   Mat_IFFT2(MatrixShen *Output_real, MatrixShen *Output_imag, MatrixShen *Input_real, MatrixShen *Input_imag);
void   four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign);
void   four1(double *data, int nn, int isign);
void   Mat_Copy(MatrixShen *A, MatrixShen *B, int h_target, int w_target, int h_begin, int w_begin, int h_end,int w_end);
void   Mat_uc_Copy(uc_Matrix *A, uc_Matrix *B, int h_target, int w_target, int h_begin, int w_begin,int h_end, int w_end);
void   Mat_i_Copy(i_Matrix *A, i_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end);
void   Mat_Product(MatrixShen *A, MatrixShen *B, MatrixShen *C);
void   Mat_Sum(MatrixShen *A, MatrixShen *B, MatrixShen *C);
void   Mat_Substract(MatrixShen *A, MatrixShen *B, MatrixShen *C);
void   Mat_Fliplr(MatrixShen *A);
void   Mat_Flipud(MatrixShen *A);
void   Mat_uc_Fliplr(uc_Matrix *A);
void   Mat_uc_Flipud(uc_Matrix *A);

/*by SHEN in JHU*/
int    *ivectorSHEN(long nl, long nh) ;
void   free_ivectorSHEN(int *v, long nl, long nh) ;
void   Mat_Inverse(MatrixShen *A, MatrixShen *B) ; /*A in, B out*/
void   Mat_times_Vector(float *Vout, MatrixShen *A, float *Vin) ;
void   Mat_A_equal_BxC(MatrixShen *A, MatrixShen *B, MatrixShen *C) ;
void   Mat_EqualCopy(MatrixShen *A, MatrixShen *B) ;
void   Mat_Print(MatrixShen *A) ;
void   Mat_Calculate_EigenVectors_EigenValues(MatrixShen *C, float *EigenValue, MatrixShen *EigenVector, int PRNorNOT) ;
void   svdcmp(float **a, int m, int n, float w[], float **v) ;
void   vector_Print(float *v, int size) ;
float  gasdev(long *idum) ;

/* June 2001 */
ImgAttribute ****ImgAttributealloc4d(int i_size,int j_size,int k_size, int t_size);
ImgAttribute ***ImgAttributealloc3d(int i_size,int j_size,int k_size) ;
ImgAttribute *ImgAttributealloc1d(int k_size) ;
void ImgAttributefree4d(ImgAttribute ****array,int t_size,int k_size,int i_size);
void ImgAttributefree3d(ImgAttribute ***array,int k_size,int i_size) ;

/* June 2001 */
HeadImgAttribute ***HeadImgAttributealloc3d(int i_size,int j_size,int k_size) ;
HeadImgAttribute *HeadImgAttributealloc1d(int k_size) ;
void HeadImgAttributefree3d(HeadImgAttribute ***array,int k_size,int i_size) ;

/* Feb 2002 */
DTIattribute ***DTIattributealloc3d(int i_size,int j_size,int k_size) ;
DTIattribute *DTIattributealloc1d(int k_size) ;
void DTIattributefree3d(DTIattribute ***array,int k_size,int i_size) ;


//by kaiming Mar 7, 2010;
void ConvertMatrix2Shen(Matrix& newmat, MatrixShen& shenmat);
void ConvertMatrix2Shen(Matrix& newmat, MatrixShen* shenmat);
void OutputMatShen ( MatrixShen& shenmat);
void OutputMatShen  ( MatrixShen* shenmat);

float *Falloc1d(int);
int *Ialloc1d(int i_size);
float **Falloc2d(int i_size,int j_size);

void Ffree2d(float **array,int i_size);

unsigned char ***UCalloc3d(int i_size,int j_size,int k_size);

void UCfree3d(unsigned char ***array,int k_size,int i_size);

int ***Ialloc3d(int i_size,int j_size,int k_size);

void Ifree3d(int ***array,int k_size,int i_size);

int ****Ialloc4d(int i_size,int j_size,int k_size, int t_size);

void Ifree4d(int ****array,int t_size,int k_size,int i_size);

float ***Falloc3d(int i_size, int y_size, int z_size);

void Ffree3d(float ***array,int k_size,int i_size);

float ****Falloc4d(int i_size,int j_size,int k_size, int t_size);

void Ffree4d(float ****array, int i_size, int j_size, int k_size, int t_size);

unsigned short ***Salloc3d(int i_size,int j_size,int k_size);

void Sfree3d(unsigned short ***array,int k_size,int i_size);
void sort(float *Y, int *I, float *A, int length);

} //end of KML
#endif /* KAIMINGFILES_H_ */
