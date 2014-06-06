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
// 二维矢量类
//****************************************************************************
class VectorR2 {

public:
	double x, y;		// The x & y  coordinates.

public:
	// constructor
	// 默认情况下 x=0 y=0
	VectorR2( ) : x(0.0), y(0.0) {}		// 默认构造函数
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

	// 重载运算符
	VectorR2& operator+= ( const VectorR2& v ) { x+=v.x; y+=v.y; return(*this); } 
	VectorR2& operator-= ( const VectorR2& v ) { x-=v.x; y-=v.y; return(*this); }
	VectorR2& operator*= ( double m ) { x*=m; y*=m; return(*this); }
	VectorR2& operator/= ( double m ) { 
		assert(m!=0.0);	// 不可除0
		register double mInv = 1.0/m; 
		x*=mInv; y*=mInv;
		return(*this); 
	}
	// 取负
	VectorR2 operator- () const { return ( VectorR2(-x, -y) ); }


	// modify
	// 按分量逐个做出的乘积
	VectorR2& ArrayProd(const VectorR2&);		// Component-wise product
	// T1+s*u
	VectorR2& AddScaled( const VectorR2& u, double s );

	// get value
	// 模
	double Norm() const { return ( sqrt( x*x + y*y ) ); }
	// 求出X和Y绝对值的最大值(调用fabs的方法)
	double L1Norm() const { return (Max(fabs(x),fabs(y))); }
	// 与指定向量的距离
	double Dist( const VectorR2& u ) const;	
	// 与指定向量的距离的平方
	double DistSq( const VectorR2& u ) const;	
	// 求出模的平方
	double NormSq() const { return ( x*x + y*y ); }
	// 求出绝对值的最大值(不调用fabs的方法)
	double MaxAbs() const;
	// 归一化向量(没有错误检查)
	VectorR2& Normalize () { *this /= Norm(); return *this;}
	// 归一化向量(有错误检查)
	VectorR2& MakeUnit();

	// judge
	// 判断是否是单位向量(只要近似为1就OK)
	bool IsUnit( double tolerance = 1.0e-15 ) const{ 
		register double norm = Norm();
		return ( 1.0+tolerance>=norm && norm>=1.0-tolerance ); 
	}
	// 一定要为0才可以返回TRUE
	bool IsZero() const { return ( x==0.0 && y==0.0 ); }
	// 近似为0就可以返回TRUE  tolerance should be non-negative
	bool NearZero(double tolerance) const { return( MaxAbs()<=tolerance );}					

	// rotate
	// 输入弧度
	VectorR2& Rotate( double theta );	// rotate through angle theta
	// 输入已经计算好的弧度的cos和sin值
	VectorR2& Rotate( double costheta, double sintheta );
};


// 外部函数调用方法 非类调用  -----------------------------------------------------------------------------------------------------------------------------

inline VectorR2 operator+( const VectorR2& u, const VectorR2& v );
inline VectorR2 operator-( const VectorR2& u, const VectorR2& v ); 
inline VectorR2 operator*( const VectorR2& u, double m); 
inline VectorR2 operator*( double m, const VectorR2& u); 
inline VectorR2 operator/( const VectorR2& u, double m); 
inline bool operator==( const VectorR2& u, const VectorR2& v ); 

// 点乘
inline double operator^ (const VectorR2& u, const VectorR2& v ); // Dot Product
inline double InnerProduct(const VectorR2& u, const VectorR2& v ) { return (u^v); }
inline VectorR2 ArrayProd ( const VectorR2& u, const VectorR2& v );

// 叉乘
inline double CrossR2( const VectorR2& u, const VectorR2& v );	// A scalar valued cross product on two R2 vectors

// 求出向量模
inline double Mag(const VectorR2& u) { return u.Norm(); }
// 求出两个向量的距离
inline double Dist(const VectorR2& u, const VectorR2& v) { return u.Dist(v); }
// 两个向量距离的平方
inline double DistSq(const VectorR2& u, const VectorR2& v) { return u.DistSq(v); }
// 输出流
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
