///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *	OPCODE - Optimized Collision Detection
 *	Copyright (C) 2001 Pierre Terdiman
 *	Homepage: http://www.codercorner.com/Opcode.htm
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for a tree collider.
 *	\file		OPC_TreeCollider.cpp
 *	\author		Pierre Terdiman
 *	\date		March, 20, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains an AABB tree collider.
 *
 *	\class		AABBTreeCollider
 *	\author		Pierre Terdiman
 *	\version	1.0
 *	\date		March, 20, 2001
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Precompiled Header
#include "Stdafx.h"

using namespace Opcode;

//! Quickly rotates & translates a vector
__forceinline void TransformPoint(Point& dest, const Point* source, const Matrix3x3& rot, const Point& trans)
{
	dest.x = trans.x + source->x * rot.m[0][0] + source->y * rot.m[1][0] + source->z * rot.m[2][0];
	dest.y = trans.y + source->x * rot.m[0][1] + source->y * rot.m[1][1] + source->z * rot.m[2][1];
	dest.z = trans.z + source->x * rot.m[0][2] + source->y * rot.m[1][2] + source->z * rot.m[2][2];
}

//! Use CPU comparisons (comment that line to use standard FPU compares)
#define CPU_COMPARE

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	OBB-OBB overlap test using the separating axis theorem.
 *	- original code by Gomez / Gamasutra (similar to Gottschalk's one in RAPID)
 *	- optimized for AABB trees by computing the rotation matrix once (SOLID-fashion)
 *	- the fabs matrix is precomputed as well and epsilon-tweaked (RAPID-style, we found this almost mandatory)
 *	- Class III axes can be disabled... (SOLID & Intel fashion)
 *	- ...or enabled to perform some profiling
 *	- CPU comparisons used when appropriate
 *	- lazy evaluation sometimes saves some work in case of early exits (unlike SOLID)
 *
 *	\param		a			[in] extent from box A
 *	\param		Pa			[in] center from box A
 *	\param		b			[in] extent from box B
 *	\param		Pb			[in] center from box B
 *	\return		true if boxes overlap
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__forceinline bool AABBTreeCollider::BoxBoxOverlap(const Point& a, const Point& Pa, const Point& b, const Point& Pb)
{
	// Stats
	mNbBVBVTests++;

	float t,t2;

	// Class I : A's basis vectors
#ifdef CPU_COMPARE
	float Tx = (mR1to0.m[0][0]*Pb.x + mR1to0.m[1][0]*Pb.y + mR1to0.m[2][0]*Pb.z) + mT1to0.x - Pa.x;
	t = a.x + b.x*mAR.m[0][0] + b.y*mAR.m[1][0] + b.z*mAR.m[2][0];
	if(AIR(Tx) > IR(t))	return false;

	float Ty = (mR1to0.m[0][1]*Pb.x + mR1to0.m[1][1]*Pb.y + mR1to0.m[2][1]*Pb.z) + mT1to0.y - Pa.y;
	t = a.y + b.x*mAR.m[0][1] + b.y*mAR.m[1][1] + b.z*mAR.m[2][1];
	if(AIR(Ty) > IR(t))	return false;

	float Tz = (mR1to0.m[0][2]*Pb.x + mR1to0.m[1][2]*Pb.y + mR1to0.m[2][2]*Pb.z) + mT1to0.z - Pa.z;
	t = a.z + b.x*mAR.m[0][2] + b.y*mAR.m[1][2] + b.z*mAR.m[2][2];
	if(AIR(Tz) > IR(t))	return false;
#else
	float Tx = (mR1to0.m[0][0]*Pb.x + mR1to0.m[1][0]*Pb.y + mR1to0.m[2][0]*Pb.z) + mT1to0.x - Pa.x;
	t = a.x + b.x*mAR.m[0][0] + b.y*mAR.m[1][0] + b.z*mAR.m[2][0];
	if(fabsf(Tx) > t)	return false;

	float Ty = (mR1to0.m[0][1]*Pb.x + mR1to0.m[1][1]*Pb.y + mR1to0.m[2][1]*Pb.z) + mT1to0.y - Pa.y;
	t = a.y + b.x*mAR.m[0][1] + b.y*mAR.m[1][1] + b.z*mAR.m[2][1];
	if(fabsf(Ty) > t)	return false;

	float Tz = (mR1to0.m[0][2]*Pb.x + mR1to0.m[1][2]*Pb.y + mR1to0.m[2][2]*Pb.z) + mT1to0.z - Pa.z;
	t = a.z + b.x*mAR.m[0][2] + b.y*mAR.m[1][2] + b.z*mAR.m[2][2];
	if(fabsf(Tz) > t)	return false;
#endif

	// Class II : B's basis vectors
#ifdef CPU_COMPARE
	t = Tx*mR1to0.m[0][0] + Ty*mR1to0.m[0][1] + Tz*mR1to0.m[0][2];	t2 = a.x*mAR.m[0][0] + a.y*mAR.m[0][1] + a.z*mAR.m[0][2] + b.x;
	if(AIR(t)>IR(t2))	return false;

	t = Tx*mR1to0.m[1][0] + Ty*mR1to0.m[1][1] + Tz*mR1to0.m[1][2];	t2 = a.x*mAR.m[1][0] + a.y*mAR.m[1][1] + a.z*mAR.m[1][2] + b.y;
	if(AIR(t)>IR(t2))	return false;

	t = Tx*mR1to0.m[2][0] + Ty*mR1to0.m[2][1] + Tz*mR1to0.m[2][2];	t2 = a.x*mAR.m[2][0] + a.y*mAR.m[2][1] + a.z*mAR.m[2][2] + b.z;
	if(AIR(t)>IR(t2))	return false;
#else
	t = Tx*mR1to0.m[0][0] + Ty*mR1to0.m[0][1] + Tz*mR1to0.m[0][2];	t2 = a.x*mAR.m[0][0] + a.y*mAR.m[0][1] + a.z*mAR.m[0][2] + b.x;
	if(fabsf(t) > t2)	return false;

	t = Tx*mR1to0.m[1][0] + Ty*mR1to0.m[1][1] + Tz*mR1to0.m[1][2];	t2 = a.x*mAR.m[1][0] + a.y*mAR.m[1][1] + a.z*mAR.m[1][2] + b.y;
	if(fabsf(t) > t2)	return false;

	t = Tx*mR1to0.m[2][0] + Ty*mR1to0.m[2][1] + Tz*mR1to0.m[2][2];	t2 = a.x*mAR.m[2][0] + a.y*mAR.m[2][1] + a.z*mAR.m[2][2] + b.z;
	if(fabsf(t) > t2)	return false;
#endif

	// Class III : 9 cross products
	// Cool trick: always perform the full test for first level, regardless of settings.
	// That way pathological cases (such as the pencils scene) are quickly rejected anyway !
	if(mFullBoxBoxTest || mNbBVBVTests==1)
	{
#ifdef CPU_COMPARE
		t = Tz*mR1to0.m[0][1] - Ty*mR1to0.m[0][2];	t2 = a.y*mAR.m[0][2] + a.z*mAR.m[0][1] + b.y*mAR.m[2][0] + b.z*mAR.m[1][0];	if(AIR(t) > IR(t2))	return false;	// L = A0 x B0
		t = Tz*mR1to0.m[1][1] - Ty*mR1to0.m[1][2];	t2 = a.y*mAR.m[1][2] + a.z*mAR.m[1][1] + b.x*mAR.m[2][0] + b.z*mAR.m[0][0];	if(AIR(t) > IR(t2))	return false;	// L = A0 x B1
		t = Tz*mR1to0.m[2][1] - Ty*mR1to0.m[2][2];	t2 = a.y*mAR.m[2][2] + a.z*mAR.m[2][1] + b.x*mAR.m[1][0] + b.y*mAR.m[0][0];	if(AIR(t) > IR(t2))	return false;	// L = A0 x B2
		t = Tx*mR1to0.m[0][2] - Tz*mR1to0.m[0][0];	t2 = a.x*mAR.m[0][2] + a.z*mAR.m[0][0] + b.y*mAR.m[2][1] + b.z*mAR.m[1][1];	if(AIR(t) > IR(t2))	return false;	// L = A1 x B0
		t = Tx*mR1to0.m[1][2] - Tz*mR1to0.m[1][0];	t2 = a.x*mAR.m[1][2] + a.z*mAR.m[1][0] + b.x*mAR.m[2][1] + b.z*mAR.m[0][1];	if(AIR(t) > IR(t2))	return false;	// L = A1 x B1
		t = Tx*mR1to0.m[2][2] - Tz*mR1to0.m[2][0];	t2 = a.x*mAR.m[2][2] + a.z*mAR.m[2][0] + b.x*mAR.m[1][1] + b.y*mAR.m[0][1];	if(AIR(t) > IR(t2))	return false;	// L = A1 x B2
		t = Ty*mR1to0.m[0][0] - Tx*mR1to0.m[0][1];	t2 = a.x*mAR.m[0][1] + a.y*mAR.m[0][0] + b.y*mAR.m[2][2] + b.z*mAR.m[1][2];	if(AIR(t) > IR(t2))	return false;	// L = A2 x B0
		t = Ty*mR1to0.m[1][0] - Tx*mR1to0.m[1][1];	t2 = a.x*mAR.m[1][1] + a.y*mAR.m[1][0] + b.x*mAR.m[2][2] + b.z*mAR.m[0][2];	if(AIR(t) > IR(t2))	return false;	// L = A2 x B1
		t = Ty*mR1to0.m[2][0] - Tx*mR1to0.m[2][1];	t2 = a.x*mAR.m[2][1] + a.y*mAR.m[2][0] + b.x*mAR.m[1][2] + b.y*mAR.m[0][2];	if(AIR(t) > IR(t2))	return false;	// L = A2 x B2
#else
		t = Tz*mR1to0.m[0][1] - Ty*mR1to0.m[0][2];	t2 = a.y*mAR.m[0][2] + a.z*mAR.m[0][1] + b.y*mAR.m[2][0] + b.z*mAR.m[1][0];	if(fabsf(t) > t2)	return false;
		t = Tz*mR1to0.m[1][1] - Ty*mR1to0.m[1][2];	t2 = a.y*mAR.m[1][2] + a.z*mAR.m[1][1] + b.x*mAR.m[2][0] + b.z*mAR.m[0][0];	if(fabsf(t) > t2)	return false;
		t = Tz*mR1to0.m[2][1] - Ty*mR1to0.m[2][2];	t2 = a.y*mAR.m[2][2] + a.z*mAR.m[2][1] + b.x*mAR.m[1][0] + b.y*mAR.m[0][0];	if(fabsf(t) > t2)	return false;
		t = Tx*mR1to0.m[0][2] - Tz*mR1to0.m[0][0];	t2 = a.x*mAR.m[0][2] + a.z*mAR.m[0][0] + b.y*mAR.m[2][1] + b.z*mAR.m[1][1];	if(fabsf(t) > t2)	return false;
		t = Tx*mR1to0.m[1][2] - Tz*mR1to0.m[1][0];	t2 = a.x*mAR.m[1][2] + a.z*mAR.m[1][0] + b.x*mAR.m[2][1] + b.z*mAR.m[0][1];	if(fabsf(t) > t2)	return false;
		t = Tx*mR1to0.m[2][2] - Tz*mR1to0.m[2][0];	t2 = a.x*mAR.m[2][2] + a.z*mAR.m[2][0] + b.x*mAR.m[1][1] + b.y*mAR.m[0][1];	if(fabsf(t) > t2)	return false;
		t = Ty*mR1to0.m[0][0] - Tx*mR1to0.m[0][1];	t2 = a.x*mAR.m[0][1] + a.y*mAR.m[0][0] + b.y*mAR.m[2][2] + b.z*mAR.m[1][2];	if(fabsf(t) > t2)	return false;
		t = Ty*mR1to0.m[1][0] - Tx*mR1to0.m[1][1];	t2 = a.x*mAR.m[1][1] + a.y*mAR.m[1][0] + b.x*mAR.m[2][2] + b.z*mAR.m[0][2];	if(fabsf(t) > t2)	return false;
		t = Ty*mR1to0.m[2][0] - Tx*mR1to0.m[2][1];	t2 = a.x*mAR.m[2][1] + a.y*mAR.m[2][0] + b.x*mAR.m[1][2] + b.y*mAR.m[0][2];	if(fabsf(t) > t2)	return false;
#endif
	}
	return true;
}

//! Use FCOMI / FCMOV on Pentium-Pro based processors (comment that line to use plain C++)
#define USE_FCOMI

//! This macro quickly finds the min & max values among 3 variables
#define FINDMINMAX(x0, x1, x2, min, max)	\
	min = max = x0;							\
	if(x1<min) min=x1;						\
	if(x1>max) max=x1;						\
	if(x2<min) min=x2;						\
	if(x2>max) max=x2;

//! TO BE DOCUMENTED
__forceinline bool planeBoxOverlap(const Point& normal, const float d, const Point& maxbox)
{
	Point vmin, vmax;
	for(udword q=0;q<=2;q++)
	{
		if(normal[q]>0.0f)	{ vmin[q]=-maxbox[q]; vmax[q]=maxbox[q]; }
		else				{ vmin[q]=maxbox[q]; vmax[q]=-maxbox[q]; }
	}
	if((normal|vmin)+d>0.0f) return false;
	if((normal|vmax)+d>0.0f) return true;

	return false;
}

//! TO BE DOCUMENTED
#define AXISTEST_X01(a, b, fa, fb)							\
	min = a*v0.y - b*v0.z;									\
	max = a*v2.y - b*v2.z;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.y + fb * extents.z;					\
	if(min>rad || max<-rad) return false;

//! TO BE DOCUMENTED
#define AXISTEST_X2(a, b, fa, fb)							\
	min = a*v0.y - b*v0.z;									\
	max = a*v1.y - b*v1.z;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.y + fb * extents.z;					\
	if(min>rad || max<-rad) return false;

//! TO BE DOCUMENTED
#define AXISTEST_Y02(a, b, fa, fb)							\
	min = b*v0.z - a*v0.x;									\
	max = b*v2.z - a*v2.x;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.x + fb * extents.z;					\
	if(min>rad || max<-rad) return false;

//! TO BE DOCUMENTED
#define AXISTEST_Y1(a, b, fa, fb)							\
	min = b*v0.z - a*v0.x;									\
	max = b*v1.z - a*v1.x;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.x + fb * extents.z;					\
	if(min>rad || max<-rad) return false;

//! TO BE DOCUMENTED
#define AXISTEST_Z12(a, b, fa, fb)							\
	min = a*v1.x - b*v1.y;									\
	max = a*v2.x - b*v2.y;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.x + fb * extents.y;					\
	if(min>rad || max<-rad) return false;

//! TO BE DOCUMENTED
#define AXISTEST_Z0(a, b, fa, fb)							\
	min = a*v0.x - b*v0.y;									\
	max = a*v1.x - b*v1.y;									\
	if(min>max) {const float tmp=max; max=min; min=tmp;	}	\
	rad = fa * extents.x + fb * extents.y;					\
	if(min>rad || max<-rad) return false;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Triangle-Box overlap test using the separating axis theorem.
 *	This is the code from Tomas Möller, a bit optimized:
 *	- with some more lazy evaluation (faster path on PC)
 *	- with a tiny bit of assembly
 *	- with "SAT-lite" applied if needed
 *	- and perhaps with some more minor modifs...
 *
 *	\param		center		[in] box center
 *	\param		extents		[in] box extents
 *	\return		true if triangle & box overlap
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__forceinline bool AABBTreeCollider::TriBoxOverlap(const Point& center, const Point& extents)
{
	// Stats
	mNbBVPrimTests++;

	// use separating axis theorem to test overlap between triangle and box 
	// need to test for overlap in these directions: 
	// 1) the {x,y,z}-directions (actually, since we use the AABB of the triangle 
	//    we do not even need to test these) 
	// 2) normal of the triangle 
	// 3) crossproduct(edge from tri, {x,y,z}-directin) 
	//    this gives 3x3=9 more tests 

	// move everything so that the boxcenter is in (0,0,0) 
	Point v0, v1, v2;
	v0.x = mLeafVerts[0].x - center.x;
	v1.x = mLeafVerts[1].x - center.x;
	v2.x = mLeafVerts[2].x - center.x;

	// First, test overlap in the {x,y,z}-directions
#ifdef USE_FCOMI
	// find min, max of the triangle in x-direction, and test for overlap in X
	if(FCMin3(v0.x, v1.x, v2.x)>extents.x)	return false;
	if(FCMax3(v0.x, v1.x, v2.x)<-extents.x)	return false;

	// same for Y
	v0.y = mLeafVerts[0].y - center.y;
	v1.y = mLeafVerts[1].y - center.y;
	v2.y = mLeafVerts[2].y - center.y;

	if(FCMin3(v0.y, v1.y, v2.y)>extents.y)	return false;
	if(FCMax3(v0.y, v1.y, v2.y)<-extents.y)	return false;

	// same for Z
	v0.z = mLeafVerts[0].z - center.z;
	v1.z = mLeafVerts[1].z - center.z;
	v2.z = mLeafVerts[2].z - center.z;

	if(FCMin3(v0.z, v1.z, v2.z)>extents.z)	return false;
	if(FCMax3(v0.z, v1.z, v2.z)<-extents.z)	return false;
#else
	float min,max;
	// Find min, max of the triangle in x-direction, and test for overlap in X
	FINDMINMAX(v0.x, v1.x, v2.x, min, max);
	if(min>extents.x || max<-extents.x) return false;

	// Same for Y
	v0.y = mLeafVerts[0].y - center.y;
	v1.y = mLeafVerts[1].y - center.y;
	v2.y = mLeafVerts[2].y - center.y;

	FINDMINMAX(v0.y, v1.y, v2.y, min, max);
	if(min>extents.y || max<-extents.y) return false;

	// Same for Z
	v0.z = mLeafVerts[0].z - center.z;
	v1.z = mLeafVerts[1].z - center.z;
	v2.z = mLeafVerts[2].z - center.z;

	FINDMINMAX(v0.z, v1.z, v2.z, min, max);
	if(min>extents.z || max<-extents.z) return false;
#endif
	// 2) Test if the box intersects the plane of the triangle
	// compute plane equation of triangle: normal*x+d=0
	// ### could be precomputed since we use the same leaf triangle several times
	const Point e0 = v1 - v0;
	const Point e1 = v2 - v1;
	const Point normal = e0 ^ e1;
	const float d = -normal|v0;
	if(!planeBoxOverlap(normal, d, extents)) return false;

	// 3) "Class III" tests
	if(mFullPrimBoxTest)
	{
		float rad;
		float min, max;
		// compute triangle edges
		// - edges lazy evaluated to take advantage of early exits
		// - fabs precomputed (half less work, possible since extents are always >0)
		// - customized macros to take advantage of the null component
		// - axis vector discarded, possibly saves useless movs

		const float fey0 = fabsf(e0.y);
		const float fez0 = fabsf(e0.z);
		AXISTEST_X01(e0.z, e0.y, fez0, fey0);
		const float fex0 = fabsf(e0.x);
		AXISTEST_Y02(e0.z, e0.x, fez0, fex0);
		AXISTEST_Z12(e0.y, e0.x, fey0, fex0);

		const float fey1 = fabsf(e1.y);
		const float fez1 = fabsf(e1.z);
		AXISTEST_X01(e1.z, e1.y, fez1, fey1);
		const float fex1 = fabsf(e1.x);
		AXISTEST_Y02(e1.z, e1.x, fez1, fex1);
		AXISTEST_Z0(e1.y, e1.x, fey1, fex1);

		const Point e2 = mLeafVerts[0] - mLeafVerts[2];
		const float fey2 = fabsf(e2.y);
		const float fez2 = fabsf(e2.z);
		AXISTEST_X2(e2.z, e2.y, fez2, fey2);
		const float fex2 = fabsf(e2.x);
		AXISTEST_Y1(e2.z, e2.x, fez2, fex2);
		AXISTEST_Z12(e2.y, e2.x, fey2, fex2);
	}
	return true;
}

//! if USE_EPSILON_TEST is true then we do a check (if |dv|<EPSILON then dv=0.0;) else no check is done (which is less robust, but faster)
#define USE_EPSILON_TEST
#define LOCAL_EPSILON 0.000001f

//! sort so that a<=b
#define SORT(a,b)			\
	if(a>b)					\
	{						\
		const float c=a;	\
		a=b;				\
		b=c;				\
	}

//! Edge to edge test based on Franlin Antonio's gem: "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202
#define EDGE_EDGE_TEST(V0, U0, U1)						\
	Bx = U0[i0] - U1[i0];								\
	By = U0[i1] - U1[i1];								\
	Cx = V0[i0] - U0[i0];								\
	Cy = V0[i1] - U0[i1];								\
	f  = Ay*Bx - Ax*By;									\
	d  = By*Cx - Bx*Cy;									\
	if((f>0.0f && d>=0.0f && d<=f) || (f<0.0f && d<=0.0f && d>=f))	\
	{													\
		const float e=Ax*Cy - Ay*Cx;					\
		if(f>0.0f)										\
		{												\
			if(e>=0.0f && e<=f) return 1;				\
		}												\
		else											\
		{												\
			if(e<=0.0f && e>=f) return 1;				\
		}												\
	}

//! TO BE DOCUMENTED
#define EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2)		\
{														\
	float Bx,By,Cx,Cy,d,f;								\
	const float Ax = V1[i0] - V0[i0];					\
	const float Ay = V1[i1] - V0[i1];					\
	/* test edge U0,U1 against V0,V1 */					\
	EDGE_EDGE_TEST(V0, U0, U1);							\
	/* test edge U1,U2 against V0,V1 */					\
	EDGE_EDGE_TEST(V0, U1, U2);							\
	/* test edge U2,U1 against V0,V1 */					\
	EDGE_EDGE_TEST(V0, U2, U0);							\
}

//! TO BE DOCUMENTED
#define POINT_IN_TRI(V0, U0, U1, U2)					\
{														\
	/* is T1 completly inside T2? */					\
	/* check if V0 is inside tri(U0,U1,U2) */			\
	float a  = U1[i1] - U0[i1];							\
	float b  = -(U1[i0] - U0[i0]);						\
	float c  = -a*U0[i0] - b*U0[i1];					\
	float d0 = a*V0[i0] + b*V0[i1] + c;					\
														\
	a  = U2[i1] - U1[i1];								\
	b  = -(U2[i0] - U1[i0]);							\
	c  = -a*U1[i0] - b*U1[i1];							\
	const float d1 = a*V0[i0] + b*V0[i1] + c;			\
														\
	a  = U0[i1] - U2[i1];								\
	b  = -(U0[i0] - U2[i0]);							\
	c  = -a*U2[i0] - b*U2[i1];							\
	const float d2 = a*V0[i0] + b*V0[i1] + c;			\
	if(d0*d1>0.0f)										\
	{													\
		if(d0*d2>0.0f) return 1;						\
	}													\
}

//! TO BE DOCUMENTED
bool CoplanarTriTri(const Point& n, const Point& v0, const Point& v1, const Point& v2, const Point& u0, const Point& u1, const Point& u2)
{
	float A[3];
	short i0,i1;
	/* first project onto an axis-aligned plane, that maximizes the area */
	/* of the triangles, compute indices: i0,i1. */
	A[0] = fabsf(n[0]);
	A[1] = fabsf(n[1]);
	A[2] = fabsf(n[2]);
	if(A[0]>A[1])
	{
		if(A[0]>A[2])
		{
			i0=1;      /* A[0] is greatest */
			i1=2;
		}
		else
		{
			i0=0;      /* A[2] is greatest */
			i1=1;
		}
	}
	else   /* A[0]<=A[1] */
	{
		if(A[2]>A[1])
		{
			i0=0;      /* A[2] is greatest */
			i1=1;
		}
		else
		{
			i0=0;      /* A[1] is greatest */
			i1=2;
		}
	}

	/* test all edges of triangle 1 against the edges of triangle 2 */
	EDGE_AGAINST_TRI_EDGES(v0, v1, u0, u1, u2);
	EDGE_AGAINST_TRI_EDGES(v1, v2, u0, u1, u2);
	EDGE_AGAINST_TRI_EDGES(v2, v0, u0, u1, u2);

	/* finally, test if tri1 is totally contained in tri2 or vice versa */
	POINT_IN_TRI(v0, u0, u1, u2);
	POINT_IN_TRI(u0, v0, v1, v2);

	return 0;
}

//! TO BE DOCUMENTED
#define NEWCOMPUTE_INTERVALS(VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2, A, B, C, X0, X1)	\
{																						\
	if(D0D1>0.0f)																		\
	{																					\
		/* here we know that D0D2<=0.0 */												\
		/* that is D0, D1 are on the same side, D2 on the other or on the plane */		\
		A=VV2; B=(VV0 - VV2)*D2; C=(VV1 - VV2)*D2; X0=D2 - D0; X1=D2 - D1;				\
	}																					\
	else if(D0D2>0.0f)																	\
	{																					\
		/* here we know that d0d1<=0.0 */												\
		A=VV1; B=(VV0 - VV1)*D1; C=(VV2 - VV1)*D1; X0=D1 - D0; X1=D1 - D2;				\
	}																					\
	else if(D1*D2>0.0f || D0!=0.0f)														\
	{																					\
		/* here we know that d0d1<=0.0 or that D0!=0.0 */								\
		A=VV0; B=(VV1 - VV0)*D0; C=(VV2 - VV0)*D0; X0=D0 - D1; X1=D0 - D2;				\
	}																					\
	else if(D1!=0.0f)																	\
	{																					\
		A=VV1; B=(VV0 - VV1)*D1; C=(VV2 - VV1)*D1; X0=D1 - D0; X1=D1 - D2;				\
	}																					\
	else if(D2!=0.0f)																	\
	{																					\
		A=VV2; B=(VV0 - VV2)*D2; C=(VV1 - VV2)*D2; X0=D2 - D0; X1=D2 - D1;				\
	}																					\
	else																				\
	{																					\
		/* triangles are coplanar */													\
		return CoplanarTriTri(N1, V0, V1, V2, U0, U1, U2);								\
	}																					\
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Triangle/triangle intersection test routine,
 *	by Tomas Moller, 1997.
 *	See article "A Fast Triangle-Triangle Intersection Test",
 *	Journal of Graphics Tools, 2(2), 1997
 *
 *	Updated June 1999: removed the divisions -- a little faster now!
 *	Updated October 1999: added {} to CROSS and SUB macros 
 *
 *	int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
 *                      float U0[3],float U1[3],float U2[3])
 *
 *	\param		V0		[in] triangle 0, vertex 0
 *	\param		V1		[in] triangle 0, vertex 1
 *	\param		V2		[in] triangle 0, vertex 2
 *	\param		U0		[in] triangle 1, vertex 0
 *	\param		U1		[in] triangle 1, vertex 1
 *	\param		U2		[in] triangle 1, vertex 2
 *	\return		true if triangles overlap
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__forceinline bool AABBTreeCollider::TriTriOverlap(const Point& V0, const Point& V1, const Point& V2, const Point& U0, const Point& U1, const Point& U2)
{
	// Stats
	mNbPrimPrimTests++;

	// Compute plane equation of triangle(V0,V1,V2)
	Point E1 = V1 - V0;
	Point E2 = V2 - V0;
	const Point N1 = E1 ^ E2;
	const float d1 =-N1 | V0;
	// Plane equation 1: N1.X+d1=0

	// Put U0,U1,U2 into plane equation 1 to compute signed distances to the plane
	float du0 = (N1|U0) + d1;
	float du1 = (N1|U1) + d1;
	float du2 = (N1|U2) + d1;

	// Coplanarity robustness check
#ifdef USE_EPSILON_TEST
	if(fabsf(du0)<LOCAL_EPSILON) du0 = 0.0f;
	if(fabsf(du1)<LOCAL_EPSILON) du1 = 0.0f;
	if(fabsf(du2)<LOCAL_EPSILON) du2 = 0.0f;
#endif
	const float du0du1 = du0 * du1;
	const float du0du2 = du0 * du2;

	if(du0du1>0.0f && du0du2>0.0f)	// same sign on all of them + not equal 0 ?
		return 0;					// no intersection occurs

	// Compute plane of triangle (U0,U1,U2)
	E1 = U1 - U0;
	E2 = U2 - U0;
	const Point N2 = E1 ^ E2;
	const float d2=-N2 | U0;
	// plane equation 2: N2.X+d2=0

	// put V0,V1,V2 into plane equation 2
	float dv0 = (N2|V0) + d2;
	float dv1 = (N2|V1) + d2;
	float dv2 = (N2|V2) + d2;

#ifdef USE_EPSILON_TEST
	if(fabsf(dv0)<LOCAL_EPSILON) dv0 = 0.0f;
	if(fabsf(dv1)<LOCAL_EPSILON) dv1 = 0.0f;
	if(fabsf(dv2)<LOCAL_EPSILON) dv2 = 0.0f;
#endif

	const float dv0dv1 = dv0 * dv1;
	const float dv0dv2 = dv0 * dv2;

	if(dv0dv1>0.0f && dv0dv2>0.0f)	// same sign on all of them + not equal 0 ?
		return 0;					// no intersection occurs

	// Compute direction of intersection line
	const Point D = N1^N2;

	// Compute and index to the largest component of D
	float max=fabsf(D[0]);
	short index=0;
	float bb=fabsf(D[1]);
	float cc=fabsf(D[2]);
	if(bb>max) max=bb,index=1;
	if(cc>max) max=cc,index=2;

	// This is the simplified projection onto L
	const float vp0 = V0[index];
	const float vp1 = V1[index];
	const float vp2 = V2[index];

	const float up0 = U0[index];
	const float up1 = U1[index];
	const float up2 = U2[index];

	// Compute interval for triangle 1
	float a,b,c,x0,x1;
	NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1);

	// Compute interval for triangle 2
	float d,e,f,y0,y1;
	NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);

	const float xx=x0*x1;
	const float yy=y0*y1;
	const float xxyy=xx*yy;

	float isect1[2], isect2[2];

	float tmp=a*xxyy;
	isect1[0]=tmp+b*x1*yy;
	isect1[1]=tmp+c*x0*yy;

	tmp=d*xxyy;
	isect2[0]=tmp+e*xx*y1;
	isect2[1]=tmp+f*xx*y0;

	SORT(isect1[0],isect1[1]);
	SORT(isect2[0],isect2[1]);

	if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
	return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Constructor.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AABBTreeCollider::AABBTreeCollider()
 :	mNbBVBVTests(0),
	mNbPrimPrimTests(0),
	mNbBVPrimTests(0),
	mFullBoxBoxTest(true),
	mFullPrimBoxTest(true),
	mFirstContact(false),
	mTemporalCoherence(false),
	mUserData0(0),
	mUserData1(0),
	mObj0Callback(null),
	mObj1Callback(null)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Destructor.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AABBTreeCollider::~AABBTreeCollider()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Generic collision query for generic OPCODE models. After the call, access the results with:
 *	- GetContactStatus()
 *	- GetNbPairs()
 *	- GetPairs()
 *
 *	\param		cache			[in] collision cache for model pointers and a colliding pair of primitives
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\return		true if success
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::Collide(BVTCache& cache, const Matrix4x4& world0, const Matrix4x4& world1)
{
	// Checkings
	if(!cache.Model0 || !cache.Model1)								return false;
	if(cache.Model0->HasLeafNodes()!=cache.Model1->HasLeafNodes())	return false;
	if(cache.Model0->IsQuantized()!=cache.Model1->IsQuantized())	return false;

	// Simple double-dispatch
	if(!cache.Model0->HasLeafNodes())
	{
		if(cache.Model0->IsQuantized())
		{
			const AABBQuantizedNoLeafTree* T0 = (const AABBQuantizedNoLeafTree*)cache.Model0->GetTree();
			const AABBQuantizedNoLeafTree* T1 = (const AABBQuantizedNoLeafTree*)cache.Model1->GetTree();
			return Collide(T0, T1, world0, world1, &cache);
		}
		else
		{
			const AABBNoLeafTree* T0 = (const AABBNoLeafTree*)cache.Model0->GetTree();
			const AABBNoLeafTree* T1 = (const AABBNoLeafTree*)cache.Model1->GetTree();
			return Collide(T0, T1, world0, world1, &cache);
		}
	}
	else
	{
		if(cache.Model0->IsQuantized())
		{
			const AABBQuantizedTree* T0 = (const AABBQuantizedTree*)cache.Model0->GetTree();
			const AABBQuantizedTree* T1 = (const AABBQuantizedTree*)cache.Model1->GetTree();
			return Collide(T0, T1, world0, world1, &cache);
		}
		else
		{
			const AABBCollisionTree* T0 = (const AABBCollisionTree*)cache.Model0->GetTree();
			const AABBCollisionTree* T1 = (const AABBCollisionTree*)cache.Model1->GetTree();
			return Collide(T0, T1, world0, world1, &cache);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	A method to initialize a collision query :
 *	- reset stats & contact status
 *	- setup matrices
 *
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::InitQuery(const Matrix4x4& world0, const Matrix4x4& world1)
{
	// Reset stats & contact status
	mContact			= false;
	mNbBVBVTests		= 0;
	mNbPrimPrimTests	= 0;
	mNbBVPrimTests		= 0;
	mPairs.Reset();

	// Setup matrices
	Matrix4x4 InvWorld0, InvWorld1;
	InvertPRMatrix(InvWorld0, world0);
	InvertPRMatrix(InvWorld1, world1);

	Matrix4x4 World0to1 = world0 * InvWorld1;
	Matrix4x4 World1to0 = world1 * InvWorld0;

	mR0to1 = World0to1;		World0to1.GetTrans(mT0to1);
	mR1to0 = World1to0;		World1to0.GetTrans(mT1to0);

	// Precompute absolute 1-to-0 rotation matrix
	for(udword i=0;i<3;i++)
	{
		for(udword j=0;j<3;j++)
		{
			// Epsilon value prevents floating-point inaccuracies (strategy borrowed from RAPID)
			mAR.m[i][j] = 1e-6f + fabsf(mR1to0.m[i][j]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	A method to take advantage of temporal coherence.
 *	\param		cache	[in] cache for a pair of previously colliding primitives
 *	\warning	only works for "First Contact" mode
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::CheckTemporalCoherence(Pair* cache)
{
	// Checkings
	if(!cache)	return false;

	// Test previously colliding primitives first
	if(mTemporalCoherence && mFirstContact)
	{
		PrimTest(cache->id0, cache->id1);
		if(mContact)	return true;
	}
	return false;
}

#define UPDATE_CACHE						\
	if(cache && mContact)					\
	{										\
		cache->id0 = mPairs.GetEntry(0);	\
		cache->id1 = mPairs.GetEntry(1);	\
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Collision query for normal AABB trees.
 *	\param		tree0			[in] AABB tree from first object
 *	\param		tree1			[in] AABB tree from second object
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\param		cache			[in/out] cache for a pair of previously colliding primitives
 *	\return		true if success
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::Collide(const AABBCollisionTree* tree0, const AABBCollisionTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache)
{
	// Checkings
	if(!tree0 || !tree1)					return false;
	if(!mObj0Callback || !mObj1Callback)	return false;

	// Init collision query
	InitQuery(world0, world1);

	// Check previous state
	if(CheckTemporalCoherence(cache))		return true;

	// Perform collision query
	_Collide(tree0->GetNodes(), tree1->GetNodes());

	UPDATE_CACHE

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Collision query for no-leaf AABB trees.
 *	\param		tree0			[in] AABB tree from first object
 *	\param		tree1			[in] AABB tree from second object
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\param		cache			[in/out] cache for a pair of previously colliding primitives
 *	\return		true if success
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::Collide(const AABBNoLeafTree* tree0, const AABBNoLeafTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache)
{
	// Checkings
	if(!tree0 || !tree1)					return false;
	if(!mObj0Callback || !mObj1Callback)	return false;

	// Init collision query
	InitQuery(world0, world1);

	// Check previous state
	if(CheckTemporalCoherence(cache))		return true;

	// Perform collision query
	_Collide(tree0->GetNodes(), tree1->GetNodes());

	UPDATE_CACHE

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Collision query for quantized AABB trees.
 *	\param		tree0			[in] AABB tree from first object
 *	\param		tree1			[in] AABB tree from second object
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\param		cache			[in/out] cache for a pair of previously colliding primitives
 *	\return		true if success
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::Collide(const AABBQuantizedTree* tree0, const AABBQuantizedTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache)
{
	// Checkings
	if(!tree0 || !tree1)					return false;
	if(!mObj0Callback || !mObj1Callback)	return false;

	// Init collision query
	InitQuery(world0, world1);

	// Check previous state
	if(CheckTemporalCoherence(cache))		return true;

	// Setup dequantization coeffs
	mCenterCoeff0	= tree0->mCenterCoeff;
	mExtentsCoeff0	= tree0->mExtentsCoeff;
	mCenterCoeff1	= tree1->mCenterCoeff;
	mExtentsCoeff1	= tree1->mExtentsCoeff;

	// Dequantize box A
	const AABBQuantizedNode* N0 = tree0->GetNodes();
	const Point a(float(N0->mAABB.mExtents[0]) * mExtentsCoeff0.x, float(N0->mAABB.mExtents[1]) * mExtentsCoeff0.y, float(N0->mAABB.mExtents[2]) * mExtentsCoeff0.z);
	const Point Pa(float(N0->mAABB.mCenter[0]) * mCenterCoeff0.x, float(N0->mAABB.mCenter[1]) * mCenterCoeff0.y, float(N0->mAABB.mCenter[2]) * mCenterCoeff0.z);
	// Dequantize box B
	const AABBQuantizedNode* N1 = tree1->GetNodes();
	const Point b(float(N1->mAABB.mExtents[0]) * mExtentsCoeff1.x, float(N1->mAABB.mExtents[1]) * mExtentsCoeff1.y, float(N1->mAABB.mExtents[2]) * mExtentsCoeff1.z);
	const Point Pb(float(N1->mAABB.mCenter[0]) * mCenterCoeff1.x, float(N1->mAABB.mCenter[1]) * mCenterCoeff1.y, float(N1->mAABB.mCenter[2]) * mCenterCoeff1.z);

	// Perform collision query
	_Collide(N0, N1, a, Pa, b, Pb);

	UPDATE_CACHE

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Collision query for quantized no-leaf AABB trees.
 *	\param		tree0			[in] AABB tree from first object
 *	\param		tree1			[in] AABB tree from second object
 *	\param		world0			[in] world matrix for first object
 *	\param		world1			[in] world matrix for second object
 *	\param		cache			[in/out] cache for a pair of previously colliding primitives
 *	\return		true if success
 *	\warning	SCALE NOT SUPPORTED. The matrices must contain rotation & translation parts only.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool AABBTreeCollider::Collide(const AABBQuantizedNoLeafTree* tree0, const AABBQuantizedNoLeafTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache)
{
	// Checkings
	if(!tree0 || !tree1)					return false;
	if(!mObj0Callback || !mObj1Callback)	return false;

	// Init collision query
	InitQuery(world0, world1);

	// Check previous state
	if(CheckTemporalCoherence(cache))		return true;

	// Setup dequantization coeffs
	mCenterCoeff0	= tree0->mCenterCoeff;
	mExtentsCoeff0	= tree0->mExtentsCoeff;
	mCenterCoeff1	= tree1->mCenterCoeff;
	mExtentsCoeff1	= tree1->mExtentsCoeff;

	// Perform collision query
	_Collide(tree0->GetNodes(), tree1->GetNodes());

	UPDATE_CACHE

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard trees
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The normal AABB tree can use 2 different descent rules (with different performances)
//#define ORIGINAL_CODE			//!< UNC-like descent rules
#define ALTERNATIVE_CODE		//!< Alternative descent rules

#ifdef ORIGINAL_CODE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision query for normal AABB trees.
 *	\param		b0		[in] collision node from first tree
 *	\param		b1		[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_Collide(const AABBCollisionNode* b0, const AABBCollisionNode* b1)
{
	// Perform BV-BV overlap test
	if(!BoxBoxOverlap(b0->mAABB.mExtents, b0->mAABB.mCenter, b1->mAABB.mExtents, b1->mAABB.mCenter))	return;

	if(b0->IsLeaf() && b1->IsLeaf()) { PrimTest(b0->GetPrimitive(), b1->GetPrimitive()); return; }

	if(b1->IsLeaf() || (!b0->IsLeaf() && (b0->GetSize() > b1->GetSize())))
	{
		_Collide(b0->GetNeg(), b1);
		if(mFirstContact && mContact) return;
		_Collide(b0->GetPos(), b1);
	}
	else
	{
		_Collide(b0, b1->GetNeg());
		if(mFirstContact && mContact) return;
		_Collide(b0, b1->GetPos());
	}
}
#endif

#ifdef ALTERNATIVE_CODE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision query for normal AABB trees.
 *	\param		b0		[in] collision node from first tree
 *	\param		b1		[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_Collide(const AABBCollisionNode* b0, const AABBCollisionNode* b1)
{
	// Perform BV-BV overlap test
	if(!BoxBoxOverlap(b0->mAABB.mExtents, b0->mAABB.mCenter, b1->mAABB.mExtents, b1->mAABB.mCenter))	return;

	if(b0->IsLeaf())
	{
		if(b1->IsLeaf())
		{
			PrimTest(b0->GetPrimitive(), b1->GetPrimitive());
		}
		else
		{
			_Collide(b0, b1->GetNeg());
			if(mFirstContact && mContact) return;
			_Collide(b0, b1->GetPos());
		}
	}
	else if(b1->IsLeaf())
	{
		_Collide(b0->GetNeg(), b1);
		if(mFirstContact && mContact) return;
		_Collide(b0->GetPos(), b1);
	}
	else
	{
		_Collide(b0->GetNeg(), b1->GetNeg());
		if(mFirstContact && mContact) return;
		_Collide(b0->GetNeg(), b1->GetPos());
		if(mFirstContact && mContact) return;
		_Collide(b0->GetPos(), b1->GetNeg());
		if(mFirstContact && mContact) return;
		_Collide(b0->GetPos(), b1->GetPos());
	}
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// No-leaf trees
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Leaf-leaf test for two primitive indices.
 *	\param		id0		[in] index from first leaf-triangle
 *	\param		id1		[in] index from second leaf-triangle
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::PrimTest(udword id0, udword id1)
{
	// Request vertices from the app
	VertexPointers VP0;	(mObj0Callback)(id0, VP0, mUserData0);
	VertexPointers VP1;	(mObj1Callback)(id1, VP1, mUserData1);

	// Transform from space 1 to space 0
	Point u0,u1,u2;
	TransformPoint(u0, VP1.Vertex[0], mR1to0, mT1to0);
	TransformPoint(u1, VP1.Vertex[1], mR1to0, mT1to0);
	TransformPoint(u2, VP1.Vertex[2], mR1to0, mT1to0);

	// Perform triangle-triangle overlap test
	if(TriTriOverlap(*VP0.Vertex[0], *VP0.Vertex[1], *VP0.Vertex[2], u0, u1, u2))
	{
		// Keep track of colliding pairs
		mPairs.Add(id0).Add(id1);
		// Set contact status
		mContact = true;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Leaf-leaf test for a previously fetched triangle from tree A (in B's space) and a new leaf from B.
 *	\param		id1		[in] leaf-triangle index from tree B
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__forceinline void AABBTreeCollider::PrimTestTriIndex(udword id1)
{
	// Request vertices from the app
	VertexPointers VP;	(mObj1Callback)(id1, VP, mUserData1);

	// Perform triangle-triangle overlap test
	if(TriTriOverlap(mLeafVerts[0], mLeafVerts[1], mLeafVerts[2], *VP.Vertex[0], *VP.Vertex[1], *VP.Vertex[2]))
	{
		// Keep track of colliding pairs
		mPairs.Add(mLeafIndex).Add(id1);
		// Set contact status
		mContact = true;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Leaf-leaf test for a previously fetched triangle from tree B (in A's space) and a new leaf from A.
 *	\param		id0		[in] leaf-triangle index from tree A
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__forceinline void AABBTreeCollider::PrimTestIndexTri(udword id0)
{
	// Request vertices from the app
	VertexPointers VP;	(mObj0Callback)(id0, VP, mUserData0);

	// Perform triangle-triangle overlap test
	if(TriTriOverlap(mLeafVerts[0], mLeafVerts[1], mLeafVerts[2], *VP.Vertex[0], *VP.Vertex[1], *VP.Vertex[2]))
	{
		// Keep track of colliding pairs
		mPairs.Add(id0).Add(mLeafIndex);
		// Set contact status
		mContact = true;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision of a leaf node from A and a branch from B.
 *	\param		b		[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_CollideTriBox(const AABBNoLeafNode* b)
{
	// Perform triangle-box overlap test
	if(!TriBoxOverlap(b->mAABB.mCenter, b->mAABB.mExtents))	return;

	// Keep same triangle, deal with first child
	if(b->HasLeaf())	PrimTestTriIndex(b->GetPrimitive());
	else				_CollideTriBox(b->GetPos());

	if(mFirstContact && mContact) return;

	// Keep same triangle, deal with second child
	if(b->HasLeaf2())	PrimTestTriIndex(b->GetPrimitive2());
	else				_CollideTriBox(b->GetNeg());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision of a leaf node from B and a branch from A.
 *	\param		b		[in] collision node from first tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_CollideBoxTri(const AABBNoLeafNode* b)
{
	// Perform triangle-box overlap test
	if(!TriBoxOverlap(b->mAABB.mCenter, b->mAABB.mExtents))	return;

	// Keep same triangle, deal with first child
	if(b->HasLeaf())	PrimTestIndexTri(b->GetPrimitive());
	else				_CollideBoxTri(b->GetPos());

	if(mFirstContact && mContact) return;

	// Keep same triangle, deal with second child
	if(b->HasLeaf2())	PrimTestIndexTri(b->GetPrimitive2());
	else				_CollideBoxTri(b->GetNeg());
}

//! Request triangle vertices from the app and transform them
#define FETCH_LEAF(primindex, callback, userdata, rot, trans)	\
	mLeafIndex = primindex;										\
	/* Request vertices from the app */							\
	VertexPointers VP;	(callback)(primindex, VP, userdata);	\
	/* Transform them in a common space */						\
	TransformPoint(mLeafVerts[0], VP.Vertex[0], rot, trans);	\
	TransformPoint(mLeafVerts[1], VP.Vertex[1], rot, trans);	\
	TransformPoint(mLeafVerts[2], VP.Vertex[2], rot, trans);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision query for no-leaf AABB trees.
 *	\param		a	[in] collision node from first tree
 *	\param		b	[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_Collide(const AABBNoLeafNode* a, const AABBNoLeafNode* b)
{
	// Perform BV-BV overlap test
	if(!BoxBoxOverlap(a->mAABB.mExtents, a->mAABB.mCenter, b->mAABB.mExtents, b->mAABB.mCenter))	return;

	// Catch leaf status
	BOOL BHasLeaf	= b->HasLeaf();
	BOOL BHasLeaf2	= b->HasLeaf2();

	if(a->HasLeaf())
	{
		FETCH_LEAF(a->GetPrimitive(), mObj0Callback, mUserData0, mR0to1, mT0to1)

		if(BHasLeaf)	PrimTestTriIndex(b->GetPrimitive());
		else			_CollideTriBox(b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)	PrimTestTriIndex(b->GetPrimitive2());
		else			_CollideTriBox(b->GetNeg());
	}
	else
	{
		if(BHasLeaf)
		{
			FETCH_LEAF(b->GetPrimitive(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetPos());
		}
		else _Collide(a->GetPos(), b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)
		{
			FETCH_LEAF(b->GetPrimitive2(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetPos());
		}
		else _Collide(a->GetPos(), b->GetNeg());
	}

	if(mFirstContact && mContact) return;

	if(a->HasLeaf2())
	{
		FETCH_LEAF(a->GetPrimitive2(), mObj0Callback, mUserData0, mR0to1, mT0to1)

		if(BHasLeaf)	PrimTestTriIndex(b->GetPrimitive());
		else			_CollideTriBox(b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)	PrimTestTriIndex(b->GetPrimitive2());
		else			_CollideTriBox(b->GetNeg());
	}
	else
	{
		if(BHasLeaf)
		{
			// ### That leaf has possibly already been fetched
			FETCH_LEAF(b->GetPrimitive(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetNeg());
		}
		else _Collide(a->GetNeg(), b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)
		{
			// ### That leaf has possibly already been fetched
			FETCH_LEAF(b->GetPrimitive2(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetNeg());
		}
		else _Collide(a->GetNeg(), b->GetNeg());
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Quantized trees
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision query for quantized AABB trees.
 *	\param		b0		[in] collision node from first tree
 *	\param		b1		[in] collision node from second tree
 *	\param		a		[in] extent from box A
 *	\param		Pa		[in] center from box A
 *	\param		b		[in] extent from box B
 *	\param		Pb		[in] center from box B
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_Collide(const AABBQuantizedNode* b0, const AABBQuantizedNode* b1, const Point& a, const Point& Pa, const Point& b, const Point& Pb)
{
	// Perform BV-BV overlap test
	if(!BoxBoxOverlap(a, Pa, b, Pb))	return;

	if(b0->IsLeaf() && b1->IsLeaf()) { PrimTest(b0->GetPrimitive(), b1->GetPrimitive()); return; }

	if(b1->IsLeaf() || (!b0->IsLeaf() && (b0->GetSize() > b1->GetSize())))
	{
		// Dequantize box
		const QuantizedAABB* Box = &b0->GetNeg()->mAABB;
		const Point negPa(float(Box->mCenter[0]) * mCenterCoeff0.x, float(Box->mCenter[1]) * mCenterCoeff0.y, float(Box->mCenter[2]) * mCenterCoeff0.z);
		const Point nega(float(Box->mExtents[0]) * mExtentsCoeff0.x, float(Box->mExtents[1]) * mExtentsCoeff0.y, float(Box->mExtents[2]) * mExtentsCoeff0.z);
		_Collide(b0->GetNeg(), b1, nega, negPa, b, Pb);

		if(mFirstContact && mContact) return;

		// Dequantize box
		Box = &b0->GetPos()->mAABB;
		const Point posPa(float(Box->mCenter[0]) * mCenterCoeff0.x, float(Box->mCenter[1]) * mCenterCoeff0.y, float(Box->mCenter[2]) * mCenterCoeff0.z);
		const Point posa(float(Box->mExtents[0]) * mExtentsCoeff0.x, float(Box->mExtents[1]) * mExtentsCoeff0.y, float(Box->mExtents[2]) * mExtentsCoeff0.z);
		_Collide(b0->GetPos(), b1, posa, posPa, b, Pb);
	}
	else
	{
		// Dequantize box
		const QuantizedAABB* Box = &b1->GetNeg()->mAABB;
		const Point negPb(float(Box->mCenter[0]) * mCenterCoeff1.x, float(Box->mCenter[1]) * mCenterCoeff1.y, float(Box->mCenter[2]) * mCenterCoeff1.z);
		const Point negb(float(Box->mExtents[0]) * mExtentsCoeff1.x, float(Box->mExtents[1]) * mExtentsCoeff1.y, float(Box->mExtents[2]) * mExtentsCoeff1.z);
		_Collide(b0, b1->GetNeg(), a, Pa, negb, negPb);

		if(mFirstContact && mContact) return;

		// Dequantize box
		Box = &b1->GetPos()->mAABB;
		const Point posPb(float(Box->mCenter[0]) * mCenterCoeff1.x, float(Box->mCenter[1]) * mCenterCoeff1.y, float(Box->mCenter[2]) * mCenterCoeff1.z);
		const Point posb(float(Box->mExtents[0]) * mExtentsCoeff1.x, float(Box->mExtents[1]) * mExtentsCoeff1.y, float(Box->mExtents[2]) * mExtentsCoeff1.z);
		_Collide(b0, b1->GetPos(), a, Pa, posb, posPb);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Quantized no-leaf trees
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision of a leaf node from A and a quantized branch from B.
 *	\param		leaf	[in] leaf triangle from first tree
 *	\param		b		[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_CollideTriBox(const AABBQuantizedNoLeafNode* b)
{
	// Dequantize box
	const QuantizedAABB* bb = &b->mAABB;
	const Point Pb(float(bb->mCenter[0]) * mCenterCoeff1.x, float(bb->mCenter[1]) * mCenterCoeff1.y, float(bb->mCenter[2]) * mCenterCoeff1.z);
	const Point eb(float(bb->mExtents[0]) * mExtentsCoeff1.x, float(bb->mExtents[1]) * mExtentsCoeff1.y, float(bb->mExtents[2]) * mExtentsCoeff1.z);

	// Perform triangle-box overlap test
	if(!TriBoxOverlap(Pb, eb))	return;

	if(b->HasLeaf())	PrimTestTriIndex(b->GetPrimitive());
	else				_CollideTriBox(b->GetPos());

	if(mFirstContact && mContact) return;

	if(b->HasLeaf2())	PrimTestTriIndex(b->GetPrimitive2());
	else				_CollideTriBox(b->GetNeg());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision of a leaf node from B and a quantized branch from A.
 *	\param		b		[in] collision node from first tree
 *	\param		leaf	[in] leaf triangle from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_CollideBoxTri(const AABBQuantizedNoLeafNode* b)
{
	// Dequantize box
	const QuantizedAABB* bb = &b->mAABB;
	const Point Pa(float(bb->mCenter[0]) * mCenterCoeff0.x, float(bb->mCenter[1]) * mCenterCoeff0.y, float(bb->mCenter[2]) * mCenterCoeff0.z);
	const Point ea(float(bb->mExtents[0]) * mExtentsCoeff0.x, float(bb->mExtents[1]) * mExtentsCoeff0.y, float(bb->mExtents[2]) * mExtentsCoeff0.z);

	// Perform triangle-box overlap test
	if(!TriBoxOverlap(Pa, ea))	return;

	if(b->HasLeaf())	PrimTestIndexTri(b->GetPrimitive());
	else				_CollideBoxTri(b->GetPos());

	if(mFirstContact && mContact) return;

	if(b->HasLeaf2())	PrimTestIndexTri(b->GetPrimitive2());
	else				_CollideBoxTri(b->GetNeg());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Recursive collision query for quantized no-leaf AABB trees.
 *	\param		a	[in] collision node from first tree
 *	\param		b	[in] collision node from second tree
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AABBTreeCollider::_Collide(const AABBQuantizedNoLeafNode* a, const AABBQuantizedNoLeafNode* b)
{
	// Dequantize box A
	const QuantizedAABB* ab = &a->mAABB;
	const Point Pa(float(ab->mCenter[0]) * mCenterCoeff0.x, float(ab->mCenter[1]) * mCenterCoeff0.y, float(ab->mCenter[2]) * mCenterCoeff0.z);
	const Point ea(float(ab->mExtents[0]) * mExtentsCoeff0.x, float(ab->mExtents[1]) * mExtentsCoeff0.y, float(ab->mExtents[2]) * mExtentsCoeff0.z);
	// Dequantize box B
	const QuantizedAABB* bb = &b->mAABB;
	const Point Pb(float(bb->mCenter[0]) * mCenterCoeff1.x, float(bb->mCenter[1]) * mCenterCoeff1.y, float(bb->mCenter[2]) * mCenterCoeff1.z);
	const Point eb(float(bb->mExtents[0]) * mExtentsCoeff1.x, float(bb->mExtents[1]) * mExtentsCoeff1.y, float(bb->mExtents[2]) * mExtentsCoeff1.z);

	// Perform BV-BV overlap test
	if(!BoxBoxOverlap(ea, Pa, eb, Pb))	return;

	// Catch leaf status
	BOOL BHasLeaf	= b->HasLeaf();
	BOOL BHasLeaf2	= b->HasLeaf2();

	if(a->HasLeaf())
	{
		FETCH_LEAF(a->GetPrimitive(), mObj0Callback, mUserData0, mR0to1, mT0to1)

		if(BHasLeaf)	PrimTestTriIndex( b->GetPrimitive());
		else			_CollideTriBox(b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)	PrimTestTriIndex(b->GetPrimitive2());
		else			_CollideTriBox(b->GetNeg());
	}
	else
	{
		if(BHasLeaf)
		{
			FETCH_LEAF(b->GetPrimitive(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetPos());
		}
		else _Collide(a->GetPos(), b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)
		{
			FETCH_LEAF(b->GetPrimitive2(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetPos());
		}
		else _Collide(a->GetPos(), b->GetNeg());
	}

	if(mFirstContact && mContact) return;

	if(a->HasLeaf2())
	{
		FETCH_LEAF(a->GetPrimitive2(), mObj0Callback, mUserData0, mR0to1, mT0to1)

		if(BHasLeaf)	PrimTestTriIndex(b->GetPrimitive());
		else			_CollideTriBox(b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)	PrimTestTriIndex(b->GetPrimitive2());
		else			_CollideTriBox(b->GetNeg());
	}
	else
	{
		if(BHasLeaf)
		{
			// ### That leaf has possibly already been fetched
			FETCH_LEAF(b->GetPrimitive(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetNeg());
		}
		else _Collide(a->GetNeg(), b->GetPos());

		if(mFirstContact && mContact) return;

		if(BHasLeaf2)
		{
			// ### That leaf has possibly already been fetched
			FETCH_LEAF(b->GetPrimitive2(), mObj1Callback, mUserData1, mR1to0, mT1to0)
			_CollideBoxTri(a->GetNeg());
		}
		else _Collide(a->GetNeg(), b->GetNeg());
	}
}
