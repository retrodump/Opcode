// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__F5D791A3_3FDF_11D5_8B0F_0050BAC83302__INCLUDED_)
#define AFX_STDAFX_H__F5D791A3_3FDF_11D5_8B0F_0050BAC83302__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// Insert your headers here
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

//#define OPCODE_USING_ICE

#ifdef OPCODE_USING_ICE

	#include "IceCore.h"
	using namespace IceCore;

	#include "IceMaths.h"
	using namespace IceMaths;

	#include "Meshmerizer.h"
	using namespace Meshmerizer;

	#include "ZCollide.h"
	using namespace ZCollide;

#endif // OPCODE_USING_ICE

#include "Opcode.h"

// TODO: reference additional headers your program requires here

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__F5D791A3_3FDF_11D5_8B0F_0050BAC83302__INCLUDED_)
