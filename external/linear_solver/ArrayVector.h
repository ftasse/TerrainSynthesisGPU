/*
ArrayVector.h
Written by Matthew Fisher

The Vector class contains dynamic arrays of a template type.  This is equivalent to the Standard Template
Library's vector class, but for various reasons I like to have my own.
*/

#ifndef __VECTOR

#define __VECTOR

#pragma once
#include "init.h"

template <class type> class Vector
{
public:
	//Initalization
	Vector();
	Vector(UINT _Size);				//creates a vector that can store the appropriate number of elements
	Vector(const Vector<type> &V);	//creates a vector as a copy of V

	//Destruction
	~Vector();
	void FreeMemory();				//releases all dynamic memory and resizes the Vector to zero

	//Assignment
	Vector<type>& operator = (const Vector<type> &V);

	//Memory
	void Allocate(UINT _Size);		//allocates space for Size elements
	void ReSize(UINT _Size);		//allocates space for Size elements and keeps all previous data that
									//still fits inside the space, starting from data element 0 and going up.

	//Accessors
	inline type& operator [] (UINT k)
	{
#ifdef VECTOR_DEBUG
		MyAssert(k < Size);
#endif
		return Data[k];
	}

	inline type& operator [] (int k)
	{
#ifdef VECTOR_DEBUG
		MyAssert(k >= 0 && k < int(Size));
#endif
		return Data[k];
	}

	inline const type& operator [] (UINT k) const
	{
#ifdef VECTOR_DEBUG
		MyAssert(k < Size);
#endif
		return Data[k];
	}

	inline const type& operator [] (int k) const
	{
#ifdef VECTOR_DEBUG
		MyAssert(k >= 0 && k < Size);
#endif
		return Data[k];
	}

	inline UINT Length() const
	{
		return Size;						//returns the number of elements that can be stored
	}

	void PushEnd(const type &t);	//resizes the Vector to accomidate an extra element, and puts t at the end
									//of the new Vector
	void PopEnd();					//resizes the Vector to accomidate one less element, deleting the last element

private:
	type *Data;		//the local data in a standard C array
	UINT Size;		//the number of elements referred to by Data
};

#include "ArrayVector.cpp"

#endif
