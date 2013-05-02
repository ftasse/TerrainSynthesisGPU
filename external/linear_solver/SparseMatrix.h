#pragma once
#include "init.h"
#include "ArrayVector.h"

struct SparseElement {
	unsigned int Col;
	double Entry;
};

struct SparseRow {
	bool FindElement(UINT Col, SparseElement* &Output);
	Vector<SparseElement> Data;
};

class SparseMatrix {
public:
	void Init(UINT _Size);
	void PushElement(UINT Row, UINT Col, double Entry);
	void Transpose(SparseMatrix &O);

	Vector<SparseRow> Rows;
	UINT Size;
};

void Multiply(Vector<double> &Output, SparseMatrix &M, Vector<double> &V);
void Multiply(Vector<double> &Output, double D, Vector<double> &V);
void Multiply(Vector<double> &Output, Vector<double> &V, double D);
double Multiply(Vector<double> &L, Vector<double> &R);
void Add(Vector<double> &Output, Vector<double> &L, Vector<double> &R);
void Subtract(Vector<double> &Output, Vector<double> &L, Vector<double> &R);
