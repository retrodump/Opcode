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
 *	\file		OPC_TreeCollider.h
 *	\author		Pierre Terdiman
 *	\date		March, 20, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef __OPC_TREECOLLIDER_H__
#define __OPC_TREECOLLIDER_H__

	//! This structure holds cached information used by the algorithm.
	//! Two model pointers and two colliding primitives are cached. Model pointers are assigned
	//! to their respective meshes, and the pair of colliding primitives is used for temporal
	//! coherence. That is, in case temporal coherence is enabled, those two primitives are
	//! tested for overlap before everything else. If they still collide, we're done before
	//! even entering the recursive collision code.
	struct OPCODE_API BVTCache : Pair
	{
		//! Constructor
		__forceinline	BVTCache()
						{
							ResetCache();
						}

						void ResetCache()
						{
							Model0			= null;
							Model1			= null;
							id0				= 0;
							id1				= 1;
						}

		OPCODE_Model*	Model0;	//!< Model for first object
		OPCODE_Model*	Model1;	//!< Model for second object
	};

	class OPCODE_API AABBTreeCollider
	{
		public:
		// Constructor / Destructor
											AABBTreeCollider();
											~AABBTreeCollider();
		// Generic collision query

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
						bool				Collide(BVTCache& cache, const Matrix4x4& world0, const Matrix4x4& world1);

		// Collision queries
						bool				Collide(const AABBCollisionTree* tree0, const AABBCollisionTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache=null);
						bool				Collide(const AABBNoLeafTree* tree0, const AABBNoLeafTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache=null);
						bool				Collide(const AABBQuantizedTree* tree0, const AABBQuantizedTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache=null);
						bool				Collide(const AABBQuantizedNoLeafTree* tree0, const AABBQuantizedNoLeafTree* tree1, const Matrix4x4& world0, const Matrix4x4& world1, Pair* cache=null);
		// Settings

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Settings: select between full box-box tests or "SAT-lite" tests (where Class III axes are discarded)
		 *	\param		flag		[in] true for full tests, false for coarse tests
		 *	\see		SetFullPrimBoxTest(bool flag)
		 *	\see		SetFirstContact(bool flag)
		 *	\see		SetTemporalCoherence(bool flag)
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	void				SetFullBoxBoxTest(bool flag)			{ mFullBoxBoxTest		= flag;					}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Settings: select between full triangle-box tests or "SAT-lite" tests (where Class III axes are discarded)
		 *	\param		flag		[in] true for full tests, false for coarse tests
		 *	\see		SetFullBoxBoxTest(bool flag)
		 *	\see		SetFirstContact(bool flag)
		 *	\see		SetTemporalCoherence(bool flag)
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	void				SetFullPrimBoxTest(bool flag)			{ mFullPrimBoxTest		= flag;					}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Settings: reports all contacts (false) or first contact only (true)
		 *	\param		flag		[in] true for first contact, false for all contacts
		 *	\see		SetFullBoxBoxTest(bool flag)
		 *	\see		SetFullPrimBoxTest(bool flag)
		 *	\see		SetTemporalCoherence(bool flag)
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	void				SetFirstContact(bool flag)				{ mFirstContact			= flag;					}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Settings: test pairs of colliding triangles from previous frame before anything else
		 *	\param		flag		[in] true to enable temporal coherence, false to discard it
		 *	\warning	Only works in "First contact" mode, and currently wouldn't work in an N-body system (last pair is cached in AABBTreeCollider)
		 *	\see		SetFullBoxBoxTest(bool flag)
		 *	\see		SetFullPrimBoxTest(bool flag)
		 *	\see		SetFirstContact(bool flag)
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	void				SetTemporalCoherence(bool flag)			{ mTemporalCoherence	= flag;					}

		// Stats

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Stats: a method to get the number of BV-BV overlap tests after a collision query.
		 *	\see		GetNbPrimPrimTests()
		 *	\see		GetNbBVPrimTests()
		 *	\return		the number of BV-BV tests performed during last query
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	udword				GetNbBVBVTests()				const	{ return mNbBVBVTests;							}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Stats: a method to get the number of Triangle-Triangle overlap tests after a collision query.
		 *	\see		GetNbBVBVTests()
		 *	\see		GetNbBVPrimTests()
		 *	\return		the number of Triangle-Triangle tests performed during last query
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	udword				GetNbPrimPrimTests()			const	{ return mNbPrimPrimTests;						}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Stats: a method to get the number of BV-Triangle overlap tests after a collision query.
		 *	\see		GetNbBVBVTests()
		 *	\see		GetNbPrimPrimTests()
		 *	\return		the number of BV-Triangle tests performed during last query
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	udword				GetNbBVPrimTests()				const	{ return mNbBVPrimTests;						}

		// Data access

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	A method to get the number of contacts after a collision query.
		 *	\see		GetContactStatus()
		 *	\see		GetPairs()
		 *	\return		the number of contacts / colliding pairs.
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	udword				GetNbPairs()					const	{ return mPairs.GetNbEntries()>>1;				}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	A method to get the pairs of colliding triangles after a collision query.
		 *	\see		GetContactStatus()
		 *	\see		GetNbPairs()
		 *	\return		the list of colliding pairs (triangle indices)
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	Pair*				GetPairs()						const	{ return (Pair*)mPairs.GetEntries();			}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	A method to get the last collision status after a collision query.
		 *	\see		GetPairs()
		 *	\see		GetNbPairs()
		 *	\return		true if the objects overlap
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	bool				GetContactStatus()				const	{ return mContact;								}

		// Callback control

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Callback control: a method to setup user-data assigned to callback 0.
		 *	\param		data		[in] user-defined data
		 *	\see		SetUserData1(udword data)
		 *	\return		Self-Reference
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	AABBTreeCollider&	SetUserData0(udword data)				{ mUserData0	= data;			return *this;	}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Callback control: a method to setup user-data assigned to callback 1.
		 *	\param		data		[in] user-defined data
		 *	\see		SetUserData0(udword data)
		 *	\return		Self-Reference
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	AABBTreeCollider&	SetUserData1(udword data)				{ mUserData1	= data;			return *this;	}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Callback control: a method to setup callback 0. Must provide triangle-vertices for a given triangle index.
		 *	\param		callback	[in] user-defined callback
		 *	\see		SetCallbackObj1(OPC_CALLBACK callback)
		 *	\return		Self-Reference
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	AABBTreeCollider&	SetCallbackObj0(OPC_CALLBACK callback)	{ mObj0Callback	= callback;		return *this;	}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		 *	Callback control: a method to setup callback 1. Must provide triangle-vertices for a given triangle index.
		 *	\param		callback	[in] user-defined callback
		 *	\see		SetCallbackObj0(OPC_CALLBACK callback)
		 *	\return		Self-Reference
		 */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		__forceinline	AABBTreeCollider&	SetCallbackObj1(OPC_CALLBACK callback)	{ mObj1Callback	= callback;		return *this;	}

		private:
		// Colliding pairs
						Container			mPairs;				//!< Pairs of colliding primitives
		// User callback
						udword				mUserData0;			//!< User-defined data sent to callbacks
						udword				mUserData1;			//!< User-defined data sent to callbacks
						OPC_CALLBACK		mObj0Callback;		//!< Callback for object 0
						OPC_CALLBACK		mObj1Callback;		//!< Callback for object 1
		// Stats
						udword				mNbBVBVTests;		//!< Number of BV-BV tests
						udword				mNbPrimPrimTests;	//!< Number of Primitive-Primitive tests
						udword				mNbBVPrimTests;		//!< Number of BV-Primitive tests
		// Precomputed data
						Matrix3x3			mAR;				//!< Absolute rotation matrix
						Matrix3x3			mR0to1;				//!< Rotation from object0 to object1
						Matrix3x3			mR1to0;				//!< Rotation from object1 to object0
						Point				mT0to1;				//!< Translation from object0 to object1
						Point				mT1to0;				//!< Translation from object1 to object0
		// Dequantization coeffs
						Point				mCenterCoeff0;
						Point				mExtentsCoeff0;
						Point				mCenterCoeff1;
						Point				mExtentsCoeff1;
		// Leaf description
						Point				mLeafVerts[3];		//!< Triangle vertices
						udword				mLeafIndex;			//!< Triangle index
		// Settings
						bool				mFullBoxBoxTest;	//!< Perform full BV-BV tests (true) or SAT-lite tests (false)
						bool				mFullPrimBoxTest;	//!< Perform full Primitive-BV tests (true) or SAT-lite tests (false)
						bool				mFirstContact;		//!< Report all contacts (false) or only first one (true)
						bool				mTemporalCoherence;	//!< Use temporal coherence or not
		// Collision result
						bool				mContact;			//!< Final contact status after a collision query
		// Internal methods

			// Standard AABB trees
						void				_Collide(const AABBCollisionNode* b0, const AABBCollisionNode* b1);
			// Quantized AABB trees
						void				_Collide(const AABBQuantizedNode* b0, const AABBQuantizedNode* b1, const Point& a, const Point& Pa, const Point& b, const Point& Pb);
			// No-leaf AABB trees
						void				_CollideTriBox(const AABBNoLeafNode* b);
						void				_CollideBoxTri(const AABBNoLeafNode* b);
						void				_Collide(const AABBNoLeafNode* a, const AABBNoLeafNode* b);
			// Quantized no-leaf AABB trees
						void				_CollideTriBox(const AABBQuantizedNoLeafNode* b);
						void				_CollideBoxTri(const AABBQuantizedNoLeafNode* b);
						void				_Collide(const AABBQuantizedNoLeafNode* a, const AABBQuantizedNoLeafNode* b);
			// Overlap tests
						void				PrimTest(udword id0, udword id1);
						void				PrimTestTriIndex(udword id1);
						void				PrimTestIndexTri(udword id0);

						bool				BoxBoxOverlap(const Point& a, const Point& Pa, const Point& b, const Point& Pb);
						bool				TriBoxOverlap(const Point& center, const Point& extents);
						bool				TriTriOverlap(const Point& V0, const Point& V1, const Point& V2, const Point& U0, const Point& U1, const Point& U2);
			// Init methods
						void				InitQuery(const Matrix4x4& world0, const Matrix4x4& world1);
						bool				CheckTemporalCoherence(Pair* cache);
	};

#endif // __OPC_TREECOLLIDER_H__
