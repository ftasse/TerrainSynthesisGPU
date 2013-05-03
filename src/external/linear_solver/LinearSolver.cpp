#include "LinearSolver.h"

double ErrorTolerance = 0.000001;
UINT MaxConvergenceIter = 10000;

void LinearSolver::Init(UINT n)
{
	r.Allocate(n);
	q.Allocate(n);
	p.Allocate(n);
	z.Allocate(n);

	r2.Allocate(n);
	q2.Allocate(n);
	p2.Allocate(n);
	z2.Allocate(n);

	x.Allocate(n);
	b.Allocate(n);
	Temp.Allocate(n);
	DiagTerms.Allocate(n);
	M.Init(n);

	for(UINT i=0;i<n;i++)
	{
		x[i] = 0.0;
	}
}

void LinearSolver::LoadDiagTerms()
{
	int i, N = DiagTerms.Length();
	SparseElement *E;
	for(i=0;i<N;i++)
	{
		if(M.Rows[i].FindElement(i, E))
		{
			DiagTerms[i] = E->Entry;
			if(DiagTerms[i] == 0.0) DiagTerms[i] = 1.0;
		}
		else
			DiagTerms[i] = 1.0;
	}
}

void LinearSolver::CGradSolve()
{
	bool Converged = false;
	double Result, Numerator, Denominator, Alpha, Beta, Error;
	UINT Iter=0;

	//p = d
	Multiply(Temp, M, x);
	Subtract(r, b, Temp);
	for(UINT i=0;i<p.Length();i++)
	{
		p[i] = r[i];
	}

	while(!(Converged && Iter>=300))
	{
		Iter++;

		Numerator = Multiply(r, r);
		Multiply(q, M, p);
		Denominator = Multiply(p, q);

		if(Denominator == 0.0)
		{
			cout << "Zero denominator." << endl;
			return;
		}
		if(Numerator == 0.0)
		{
			cout << "Numerator Zero; no residuals." << endl;
			return;
		}

		Alpha = Numerator / Denominator;

		Multiply(Temp, Alpha, p);
		Add(x, x, Temp);

		if(Iter % 50 == 49)
		{
			Multiply(Temp, M, x);
			Subtract(r, b, Temp);
		} else {
			Multiply(Temp, Alpha, q);
			Subtract(r, r, Temp);
		}

		Result = Multiply(r, r);
		Beta = Result / Numerator;

		Multiply(Temp, Beta, p);
		Add(p, r, Temp);

		Numerator = Multiply(r, r);
		Denominator = Multiply(b, b);
		if(Denominator == 0.0)
		{
			cout << "b's zero?" << endl;
			return;
		}
		Error = Numerator / Denominator;
		//cout << "Iteration " << Iter << ", Error -> " << Error << endl;
		if(Iter > MaxConvergenceIter || Error < ErrorTolerance)
		{
			Converged = true;
		}
	}
	if(Iter > MaxConvergenceIter) cout << "Convergenced not reached; algorithm halted." << endl;
	else cout << "Convergence, " << Iter << " iterations." << endl;
}

void LinearSolver::BiCGradSolve(double _BMag)
{
	bool Converged = false;
	bool First = true;
	double Beta, Alpha, Row = 0.0, PrevRow, Numerator, Denominator, Error, BMag;
	UINT Iter=0, N = M.Rows.Length();

	LoadDiagTerms();

	BMag = Multiply(b, b);
	if(_BMag != 0.0) BMag = _BMag;
	if(BMag == 0.0)
	{
		cout << "b's zero?" << endl;
		return;
	}

	SparseMatrix M_T;
	M.Transpose(M_T);

	//p = d
	Multiply(Temp, M, x);
	Subtract(r, b, Temp);
	for(UINT i=0;i<N;i++)
	{
		r2[i] = r[i];
	}

	while(!(Converged && Iter>=500))
	{
		Iter++;

		if(Iter % 50 == 49)
		{
			Multiply(Temp, M, x);
			Subtract(r, b, Temp);
			for(UINT i=0;i<N;i++)
			{
				r2[i] = r[i];
			}
			First = true;
		}

		for(UINT i=0;i<N;i++)
		{
			z[i] = r[i] / DiagTerms[i];
			z2[i] = r2[i] / DiagTerms[i];
		}

		PrevRow = Row;
		Row = Multiply(z, r2);
		if(Row == 0.0)
		{
			cout << "Row Zero; Algorithm fails." << endl;
			return;
		}

		if(First)
		{
			for(UINT i=0;i<N;i++)
			{
				p[i] = z[i];
				p2[i] = z2[i];
			}
			First = false;
		} else {
			Beta = Row / PrevRow;
			for(UINT i=0;i<N;i++)
			{
				p[i] = z[i] + Beta * p[i];
				p2[i] = z2[i] + Beta * p2[i];
			}
		}

		Multiply(q, M, p);
		Multiply(q2, M_T, p2);

		Denominator = Multiply(p2, q);
		if(Denominator == 0.0)
		{
			cout << "Denominator zero." << endl;
			return;
		}
		Alpha = Row / Denominator;

		for(UINT i=0;i<N;i++)
		{
			x[i] += Alpha * p[i];
			r[i] -= Alpha * q[i];
			r2[i] -= Alpha * q2[i];
		}

		Numerator = Multiply(r, r);
		Error = Numerator / BMag;
		//cout << Iter << "\t" << Error << endl;
		if(Iter > MaxConvergenceIter || Error < ErrorTolerance) Converged = true;
	}
	if(Iter > MaxConvergenceIter)
	{
		cout << "Convergenced not reached; algorithm halted." << endl;
	}
	else
	{
		//cout << "Convergence, " << Iter << " iterations." << endl;
		//cout << "Error        " << Error <<endl;
	}
}
