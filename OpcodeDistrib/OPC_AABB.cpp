///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *	OPCODE - Optimized Collision Detection
 *	Copyright (C) 2001 Pierre Terdiman
 *	Homepage: http://www.codercorner.com/Opcode.htm
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains AABB-related code.
 *	\file		OPC_AABB.cpp
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	AABB class.
 *	\class		AABB
 *	\author		Pierre Terdiman
 *	\version	1.0
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Precompiled Header
#include "Stdafx.h"

using namespace Opcode;

#ifndef __MESHMERIZER_H__

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	A method to compute the sum of two AABBs.
 *	\param		aabb		[in] the other AABB
 *	\return		Self-Reference
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
AABB& AABB::Add(const AABB& aabb)
{
	// Compute new min & max values
	Point Min;	GetMin(Min);
	Point Tmp;	aabb.GetMin(Tmp);
	Min.Min(Tmp);

	Point Max;	GetMax(Max);
	aabb.GetMax(Tmp);
	Max.Max(Tmp);

	// Update this
	SetMinMax(Min, Max);
	return *this;
}

#endif // __MESHMERIZER_H__
