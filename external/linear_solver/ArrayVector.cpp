/*
ArrayVector.cpp
Written by Matthew Fisher

See ArrayVector.h for a full definition of the Vector class.
*/

#ifndef __VECTOR_CPP
#define __VECTOR_CPP

#include "ArrayVector.h"

template <class type> Vector<type>::Vector()
{
	Size = 0;
	Data = NULL;
}

template <class type> Vector<type>::Vector(UINT _Size)
{
	Size = _Size;
	Data = new type[Size];
}

template <class type> Vector<type>::~Vector()
{
	FreeMemory();
}

template <class type> void Vector<type>::FreeMemory()
{
	Size = 0;
	if(Data) delete[] Data;
	Data = NULL;
}

/*template <class type> type& Vector<type>::operator [] (int k)
{
	if(k < 0 || k >= Size)	//make sure k is in bounds
	{
		//normally I set a breakpoint here so I can see the call stack
		cout << "Out of bounds Vector access.  k = " << k << ", Size = " << Size << endl;
		exit(1);
	}
	return Data[k];
}

template <class type> const type& Vector<type>::operator [] (int k) const
{
	if(k < 0 || k >= Size)	//make sure k is in bounds
	{
		//normally I set a breakpoint here so I can see the call stack
		cout << "Out of bounds Vector access.  k = " << k << ", Size = " << Size << endl;
		exit(1);
	}
	return Data[k];
}*/

template <class type> void Vector<type>::Allocate(UINT _Size)
{
	FreeMemory();
	Size = _Size;
	Data = new type[Size];
	if(!Data)
	{
		cout << "Bad malloc" << endl;
		exit(1);	//I've never actually seen this occur, but it's always possible to run out of memory.
	}
}

template <class type> void Vector<type>::ReSize(UINT _Size)
{
	if(Size == 0) Allocate(_Size);	//ReSize is equivalent to Allocate if Size is 0.
	else {
		UINT ElementsToCopy = Size;
		type *Storage = new type[Size];				//create intermediate storage (not necessary if I rewrote this function...)
		for(UINT i=0; i<Size; i++)
		{
			Storage[i] = Data[i];	//load the current data into storage
		}
		Allocate(_Size);							//allocate new space
		if(Size < ElementsToCopy)
		{
			ElementsToCopy = Size;	//don't write more elements than we can store
		}
		for(UINT i=0;i<ElementsToCopy;i++)
		{
			Data[i] = Storage[i];	//load our new data elements from storage
		}
		delete[] Storage;									//free the storage
	}
}

template <class type> Vector<type>::Vector(const Vector<type> &V)
{
	Size = 0;
	Data = NULL;
	Allocate(V.Size);
	for(UINT i=0;i<V.Size;i++)
		Data[i] = V.Data[i];
	//memcpy(Data, V.Data, V.Size * sizeof(type));  //using memcpy is bad if the Vector stores dynamic data
}

template <class type> Vector<type>& Vector<type>::operator = (const Vector<type> &V)
{
	Allocate(V.Size);
	for(UINT i=0;i<V.Size;i++)
		Data[i] = V.Data[i];
	//memcpy(Data, V.Data, V.Size * sizeof(type));	//using memcpy is bad if the Vector stores dynamic data
	return *this;
}

template <class type> void Vector<type>:: PushEnd(const type &t)
{
	ReSize(Size+1);
	Data[Size-1] = t;
}

template <class type> void Vector<type>:: PopEnd()
{
	ReSize(Size-1);
}

#endif