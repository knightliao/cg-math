/////////////////////////////////////////////////////////////////////////////////
//       ___          ___                     ___          ___
//      /__/|        /__/\       ___         /  /\        /__/\         ___
//     |  |:|        \  \:\     /  /\       /  /:/_       \  \:\       /  /\
//     |  |:|         \  \:\   /  /:/      /  /:/ /\       \__\:\     /  /:/
//   __|  |:|     _____\__\:\ /__/::\     /  /:/_/::\  ___ /  /::\   /  /:/
//  /__/\_|:|____/__/::::::::\\__\/\:\__ /__/:/__\/\:\/__/\  /:/\:\ /  /::\
//  \  \:\/:::::/\  \:\~~\~~\/   \  \:\/\\  \:\ /~~/:/\  \:\/:/__\//__/:/\:\
//   \  \::/~~~~  \  \:\  ~~~     \__\::/ \  \:\  /:/  \  \::/     \__\/  \:\
//    \  \:\       \  \:\         /__/:/   \  \:\/:/    \  \:\          \  \:\
//     \  \:\       \  \:\        \__\/     \  \::/      \  \:\          \__\/
//      \__\/        \__\/                   \__\/        \__\/
//
//
// author		:	KNIGHT
// Description	:	
//
//
// Time			:					
// Email		:	KnightLiao@gmail.com	liaoqiqi@cad.zju.edu.cn
// Blog			:   http://knightliao.blogspot.com
////////////////////////////////////////////////////////////////////////////////////


#ifndef KNIGHT_VECTOR_2R_H
#define KNIGHT_VECTOR_2R_H

#include <assert.h>
#include <iostream>
#include "KnightMath.h"

using namespace std;

//****************************************************************************
// ��άʸ����
//****************************************************************************
class VectorR2 {

public:
	double x, y;		// The x & y  coordinates.

public:
	// constructor
	// Ĭ������� x=0 y=0
	VectorR2( ) : x(0.0), y(0.0) {}		// Ĭ�Ϲ��캯��
	VectorR2( double xVal, double yVal ): x(xVal), y(yVal) {}

	// set value
	VectorR2& SetZero() { x=0.0; y=0.0; return *this;}	// zero
	VectorR2& Set( double xx, double yy ) { x=xx; y=yy; return *this;}	// set value defined by user
	VectorR2& SetUnitX() { x=1.0; y=0.0; return *this;}	// unit one in x axis
	VectorR2& SetUnitY() { x=0.0; y=1.0; return *this;}	// unit one in y axis
	VectorR2& SetNegUnitX() { x=-1.0; y=0.0; return *this;}	// unit one in negative x axis
	VectorR2& SetNegUnitY() { x=0.0; y=-1.0; return *this;}	// unit one in negative y axis


	// load value 
	VectorR2& Load( const double* v );
	VectorR2& Load( const float* v );

	// pull out the value
	void Dump( double* v ) const;	// no allow to modify the internal value
	void Dump( float* v ) const;

	// ���������
	VectorR2& operator+= ( const VectorR2& v ) { x+=v.x; y+=v.y; return(*this); } 
	VectorR2& operator-= ( const VectorR2& v ) { x-=v.x; y-=v.y; return(*this); }
	VectorR2& operator*= ( double m ) { x*=m; y*=m; return(*this); }
	VectorR2& operator/= ( double m ) { 
		assert(m!=0.0);	// ���ɳ�0
		register double mInv = 1.0/m; 
		x*=mInv; y*=mInv;
		return(*this); 
	}
	// ȡ��
	VectorR2 operator- () const { return ( VectorR2(-x, -y) ); }


	// modify
	// ��������������ĳ˻�
	VectorR2& ArrayProd(const VectorR2&);		// Component-wise product
	// T1+s*u
	VectorR2& AddScaled( const VectorR2& u, double s );

	// get value
	// ģ
	double Norm() const { return ( sqrt( x*x + y*y ) ); }
	// ���X��Y����ֵ�����ֵ(����fabs�ķ���)
	double L1Norm() const { return (Max(fabs(x),fabs(y))); }
	// ��ָ�������ľ���
	double Dist( const VectorR2& u ) const;	
	// ��ָ�������ľ����ƽ��
	double DistSq( const VectorR2& u ) const;	
	// ���ģ��ƽ��
	double NormSq() const { return ( x*x + y*y ); }
	// �������ֵ�����ֵ(������fabs�ķ���)
	double MaxAbs() const;
	// ��һ������(û�д�����)
	VectorR2& Normalize () { *this /= Norm(); return *this;}
	// ��һ������(�д�����)
	VectorR2& MakeUnit();

	// judge
	// �ж��Ƿ��ǵ�λ����(ֻҪ����Ϊ1��OK)
	bool IsUnit( double tolerance = 1.0e-15 ) const{ 
		register double norm = Norm();
		return ( 1.0+tolerance>=norm && norm>=1.0-tolerance ); 
	}
	// һ��ҪΪ0�ſ��Է���TRUE
	bool IsZero() const { return ( x==0.0 && y==0.0 ); }
	// ����Ϊ0�Ϳ��Է���TRUE  tolerance should be non-negative
	bool NearZero(double tolerance) const { return( MaxAbs()<=tolerance );}					

	// rotate
	// ���뻡��
	VectorR2& Rotate( double theta );	// rotate through angle theta
	// �����Ѿ�����õĻ��ȵ�cos��sinֵ
	VectorR2& Rotate( double costheta, double sintheta );
};


// �ⲿ�������÷��� �������  -----------------------------------------------------------------------------------------------------------------------------

inline VectorR2 operator+( const VectorR2& u, const VectorR2& v );
inline VectorR2 operator-( const VectorR2& u, const VectorR2& v ); 
inline VectorR2 operator*( const VectorR2& u, double m); 
inline VectorR2 operator*( double m, const VectorR2& u); 
inline VectorR2 operator/( const VectorR2& u, double m); 
inline bool operator==( const VectorR2& u, const VectorR2& v ); 

// ���
inline double operator^ (const VectorR2& u, const VectorR2& v ); // Dot Product
inline double InnerProduct(const VectorR2& u, const VectorR2& v ) { return (u^v); }
inline VectorR2 ArrayProd ( const VectorR2& u, const VectorR2& v );

// ���
inline double CrossR2( const VectorR2& u, const VectorR2& v );	// A scalar valued cross product on two R2 vectors

// �������ģ
inline double Mag(const VectorR2& u) { return u.Norm(); }
// ������������ľ���
inline double Dist(const VectorR2& u, const VectorR2& v) { return u.Dist(v); }
// �������������ƽ��
inline double DistSq(const VectorR2& u, const VectorR2& v) { return u.DistSq(v); }
// �����
ostream& operator<< ( ostream& os, const VectorR2& u );
















inline VectorR2& VectorR2::Load( const double* v ) 
{
	x = *v; 
	y = *(v+1);
	return *this;
}

inline VectorR2& VectorR2::Load( const float* v ) 
{
	x = *v; 
	y = *(v+1);
	return *this;
}

inline 	void VectorR2::Dump( double* v ) const
{
	*v = x; 
	*(v+1) = y;
}

inline 	void VectorR2::Dump( float* v ) const
{
	*v = (float)x; 
	*(v+1) = (float)y;
}

inline VectorR2& VectorR2::ArrayProd (const VectorR2& v)		// Component-wise Product
{
	x *= v.x;
	y *= v.y;
	return ( *this );
}

inline VectorR2& VectorR2::MakeUnit ()	 // Convert to unit vector (or leave zero).
{
	double nSq = NormSq();
	if (nSq != 0.0) {
		*this /= sqrt(nSq);
	}
	return *this;
}


// Rotate through angle theta
inline VectorR2& VectorR2::Rotate( double theta )
{
	double costheta = cos(theta);
	double sintheta = sin(theta);
	double tempx = x*costheta - y*sintheta;
	y = y*costheta + x*sintheta;
	x = tempx;
	return *this;
}

inline VectorR2& VectorR2::Rotate( double costheta, double sintheta )
{
	double tempx = x*costheta + y*sintheta;
	y = y*costheta - x*sintheta;
	x = tempx;
	return *this;
}

inline double VectorR2::MaxAbs() const
{
	register double m;
	m = (x>=0.0) ? x : -x;
	if ( y>m ) m=y;
	else if ( -y >m ) m = -y;
	return m;
}

inline VectorR2 ArrayProd ( const VectorR2& u, const VectorR2& v )
{
	return ( VectorR2( u.x*v.x, u.y*v.y ) );
}

#endif	
