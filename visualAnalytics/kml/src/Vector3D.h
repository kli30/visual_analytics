/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef __Vector3D__H
#define __Vector3D__H
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "newmat.h"
#include <iomanip>
#include <sstream>

using namespace NEWMAT;
using namespace std;

namespace KML
{
	template<typename T>
	class Vector3D
	{
	public:

		T x;
		T y;
		T z;
		Vector3D():x(0),y(0),z(0){};
		Vector3D(T a,T b,T c):x(a),y(b),z(c){};

		template<typename U>
		Vector3D(U dx[3]):x(T(dx[0])), y(T(dx[1])),z(T(dx[2])){};
		template<typename U>
		Vector3D(Vector3D<U>& oth):x((T)oth.x), y((T)oth.y),z((T)oth.z){};

		//const vector3D can only access to const methods

		inline float Norm1stOrder() {
			float big= abs(x) > abs(y) ? abs(x) : abs(y);
			return abs(z) > big ? abs(z) : big ;
		}
		inline void Normalize() { float norm=this->Norm(); x/=norm; y/=norm; z/=norm; };
		inline float Norm() const { return sqrt(x*x+y*y+z*z);};
		inline float NormSquare() const { return x*x+y*y+z*z;};
		template<typename U>
		inline void operator = (const Vector3D<U>& oth)  { this->x=(T)oth.x; this->y=(T)oth.y; this->z=(T)oth.z;}
		template<typename U>
		inline void operator = (U value)  { this->x=(T) value; this->y=this->x; this->z=this->x;}
		inline bool operator==(const Vector3D<T>& oth) const { return x==oth.x && y==oth.y && z==oth.z ;};
		inline bool operator!=(const Vector3D<T>& oth) const {return !( this->operator ==(oth) );} ;
		Vector3D<T>&  operator-=(const Vector3D<T>& oth){
			x-=oth.x;
			y-=oth.y;
			z-=oth.z;
			return *this;
		}
		inline bool operator <(const Vector3D<T>& oth) const
		{
			if (x<oth.x)
				return true;
			else if(x==oth.x && y<oth.y)
				return true;
			else if(x==oth.x && y==oth.y &&z<oth.z)
				return true;
			else return false;
		}
		inline bool operator <=(const Vector3D<T>& oth) const
		{
			return this->operator ==(oth) || this->operator <(oth);
		}

		inline bool operator >(const Vector3D<T>& oth) const
		{
			return ! (this->operator <=(oth));
		}
		inline Vector3D<T>& operator += (const Vector3D<T>& oth){
			x+=oth.x;
			y+=oth.y;
			z+=oth.z;
			return *this;
		}
		inline Vector3D<T>& operator /= (float den){
			x/=den;
			y/=den;
			z/=den;
			return *this;
		}

 		inline bool operator >=(const Vector3D<T>& oth) const
		{
			return ! (this->operator <(oth	));
		}
		inline Vector3D<T>& operator *= (float scalar)
		{
			x*=scalar;
			y*=scalar;
			z*=scalar;
			return *this;
		}

	};

	template<typename T, typename S>
	inline Vector3D<float> operator-(const Vector3D<T>& v1, const Vector3D<S>& v2)
	{
		Vector3D<float> diff;
		diff.x=v1.x+0.f-v2.x;
		diff.y=v1.y+0.f-v2.y;
		diff.z=v1.z+0.f-v2.z;
		return diff;
	}
	template<typename T>
	inline Vector3D<T> operator+ (const  Vector3D<T>& v1, const Vector3D<T>& v2)
	{
		Vector3D<T> sum;
		sum.x=v1.x+v2.x;
		sum.y=v1.y+v2.y;
		sum.z=v1.z+v2.z;
		return sum;
	}

	template<typename T>
	inline Vector3D<T> operator* (float factor, const Vector3D<T>& v2)
	{
		Vector3D<T> mul;
		mul.x=factor*v2.x;
		mul.y=factor*v2.y;
		mul.z=factor*v2.z;
		return mul;
	}
	template<typename T>
	inline Vector3D<T> operator/ (const Vector3D<T>& v2, float den)
	{
		float factor=1/den;
		Vector3D<T> mul;
		mul.x=factor*v2.x;
		mul.y=factor*v2.y;
		mul.z=factor*v2.z;
		return mul;
	}

	template<typename T>
	inline  Vector3D<T> CrossProduct( const Vector3D<T>& v1, const Vector3D<T>& v2)
	{
		Vector3D<T> crossProduct;
		crossProduct.x=v1.y*v2.z-v2.y*v1.z;
		crossProduct.y=v1.z*v2.x-v1.x*v2.z;
		crossProduct.z=v1.x*v2.y-v1.y*v2.x;
		return crossProduct;
	}
	template<typename T>
	inline float DotProduct(const Vector3D<T>& v1, const Vector3D<T>& v2)
	{
		return (float)(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
	}

	template<typename T>
	inline float DistanceVector3D(const Vector3D<T>& v1, const Vector3D<T>& v2)
	{
		Vector3D<T> diff=v2-v1;
		return float(diff.Norm());
	}

	template<typename T>
	inline float AbsDistanceVector3D(const Vector3D<T>& v1, const Vector3D<T>& v2)
	{
		Vector3D<T> diff=v2-v1;
		return float(abs(diff.x)+abs(diff.y)+abs(diff.z));
	}



	template<typename T>
	inline float CosAngleBetweenVectors(const Vector3D<T>& vector1,  const Vector3D<T>& vector2)
	{
		Vector3D<T> diff=vector1-vector2;
		return (vector1.NormSquare()+vector2.NormSquare()-diff.NormSquare())/(2*vector1.Norm()*vector2.Norm());
	}


} // end of namespace kml; 


#endif
