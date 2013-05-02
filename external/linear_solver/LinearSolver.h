#pragma once
#include "init.h"
#include "ArrayVector.h"
#include "SparseMatrix.h"

class LinearSolver {
public:
	void Init(UINT n);
	void LoadDiagTerms();
	void CGradSolve();
	void BiCGradSolve(double _BMag = 0.0);
	inline void PushElement(UINT Row, UINT Col, double Entry)
	{
		M.PushElement(Row, Col, Entry);
	}

	Vector<double> p, p2, q, q2, r, r2, z, z2, Temp, x, b, DiagTerms;
	SparseMatrix M;
};
