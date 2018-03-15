# Microsoft Developer Studio Project File - Name="OpcodeDistrib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=OpcodeDistrib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OpcodeDistrib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OpcodeDistrib.mak" CFG="OpcodeDistrib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OpcodeDistrib - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OpcodeDistrib - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OpcodeDistrib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OPCODEDISTRIB_EXPORTS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /G6 /MD /W3 /GX- /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OPCODEDISTRIB_EXPORTS" /Yu"stdafx.h" /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 /nologo /dll /machine:I386 /out:"Release/Opcode.dll"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=xcopy release\Opcode.lib y:\lib /Q /Y	xcopy *.h Y:\inc /Q /Y	del Y:\inc\stdafx.h
# End Special Build Tool

!ELSEIF  "$(CFG)" == "OpcodeDistrib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OPCODEDISTRIB_EXPORTS" /Yu"stdafx.h" /FD /GZ  /c
# ADD CPP /nologo /G6 /MDd /W3 /Gm /GX- /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OPCODEDISTRIB_EXPORTS" /FR /Yu"stdafx.h" /FD /GZ  /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 /nologo /dll /incremental:no /debug /machine:I386 /out:"Debug/Opcode_D.dll" /pdbtype:sept
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=xcopy debug\Opcode_D.lib y:\lib /Q /Y	xcopy *.h Y:\inc /Q /Y	del Y:\inc\stdafx.h
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "OpcodeDistrib - Win32 Release"
# Name "OpcodeDistrib - Win32 Debug"
# Begin Group "API"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\OPC_Model.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Model.h
# End Source File
# Begin Source File

SOURCE=.\Opcode.cpp
# End Source File
# Begin Source File

SOURCE=.\Opcode.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.cpp
# ADD CPP /Yc"stdafx.h"
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# End Group
# Begin Group "Trees"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\OPC_AABBTree.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_AABBTree.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Common.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Common.h
# End Source File
# Begin Source File

SOURCE=.\OPC_OptimizedTree.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_OptimizedTree.h
# End Source File
# Begin Source File

SOURCE=.\OPC_TreeBuilders.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_TreeBuilders.h
# End Source File
# Begin Source File

SOURCE=.\OPC_TreeCollider.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_TreeCollider.h
# End Source File
# End Group
# Begin Group "ICE utils"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\OPC_AABB.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_AABB.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Container.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Container.h
# End Source File
# Begin Source File

SOURCE=.\OPC_FPU.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Matrix3x3.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Matrix3x3.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Matrix4x4.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Matrix4x4.h
# End Source File
# Begin Source File

SOURCE=.\OPC_MemoryMacros.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Point.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Point.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Triangle.cpp
# End Source File
# Begin Source File

SOURCE=.\OPC_Triangle.h
# End Source File
# Begin Source File

SOURCE=.\OPC_Types.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project