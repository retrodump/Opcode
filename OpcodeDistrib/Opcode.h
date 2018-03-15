///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *	OPCODE - Optimized Collision Detection
 *	Copyright (C) 2001 Pierre Terdiman
 *	Homepage: http://www.codercorner.com/Opcode.htm
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Main file for Opcode.dll.
 *	\file		Opcode.h
 *	\author		Pierre Terdiman
 *	\date		March, 20, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef __OPCODE_H__
#define __OPCODE_H__

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compilation messages
#if defined(OPCODEDISTRIB_EXPORTS)
	#pragma message("Compiling OPCODE")
#elif !defined(OPCODEDISTRIB_EXPORTS)
	#pragma message("Using OPCODE")
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Preprocessor
#ifdef OPCODEDISTRIB_EXPORTS
	#define OPCODE_API __declspec(dllexport)
#else
	#define OPCODE_API __declspec(dllimport)
#endif

#ifndef __ICECORE_H__
	#ifdef WIN32
	#include <windows.h>
	#include <windowsx.h>
	#endif // WIN32

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <ctype.h>
	#include <assert.h>
	#include <shlobj.h>
	#include <mmsystem.h>
	#include <stdarg.h>
	#include <time.h>
	#include <float.h>

	#ifndef ASSERT
		#define	ASSERT	assert
	#endif

	#define	Log
	#define	SetIceError
	#define	gEC_OutOfMemory	0
	#define	Alignment

	#include "OPC_Types.h"
	#include "OPC_FPU.h"
	#include "OPC_MemoryMacros.h"
#endif

#ifndef __ICEMATHS_H__
	#include <Math.h>
#endif

	namespace Opcode
	{
#ifndef __ICECORE_H__
		#include "OPC_Container.h"
#endif

#ifndef __ICEMATHS_H__
		#include "OPC_Point.h"
		#include "OPC_Matrix3x3.h"
		#include "OPC_Matrix4x4.h"
#endif

#ifndef __MESHMERIZER_H__
		#include "OPC_Triangle.h"
		#include "OPC_AABB.h"
#endif

		// Bulk-of-the-work
		#include "OPC_Common.h"
		#include "OPC_TreeBuilders.h"
		#include "OPC_AABBTree.h"
		#include "OPC_OptimizedTree.h"
		#include "OPC_Model.h"
		#include "OPC_TreeCollider.h"
	}

#endif // __OPCODE_H__
