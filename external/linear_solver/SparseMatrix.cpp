#include "SparseMatrix.h"

bool SparseRow::FindElement(UINT Col, SparseElement* &Output)
{
	bool MatchFound = false;
	Output = NULL;
	for(UINT i=0;i<Data.Length() && !MatchFound;i++)
	{
		 if(Data[i].Col == Col)
		 {
			 MatchFound = true;
			 Output = &Data[i];
		 }
	}
	return MatchFound;
}

void SparseMatrix::Init(UINT _Size)
{
	Size = _Size;
	Rows.Allocate(Size);
}

void SparseMatrix::PushElement(UINT Row, UINT Col, double Entry)
{
	SparseElement *Element;
	if(Col >= Size) cout << "Non-square element added." << endl;
	if(Rows[Row].FindElement(Col, Element))
	{
		Element->Entry += Entry;
	} else {
		SparseElement NewElement;
		NewElement.Col = Col;
		NewElement.Entry = Entry;
		Rows[Row].Data.PushEnd(NewElement);
	}
}

void SparseMatrix::Transpose(SparseMatrix &O)
{
	SparseElement *Element;

	O.Rows.FreeMemory();
	O.Init(Size);

	for(UINT i=0;i<Rows.Length();i++)
	{
		for(UINT i2=0;i2<Rows[i].Data.Length();i2++)
		{
			Element = &Rows[i].Data[i2];
			O.PushElement(Element->Col, i, Element->Entry);
		}
	}
}

void Multiply(Vector<double> &Output, SparseMatrix &M, Vector<double> &V)
{
	SparseElement *E;

	if(Output.Length() != V.Length())
	{
		Output.Allocate(V.Length());
	}
	for(UINT i=0;i<Output.Length();i++)
	{
		double Total = 0.0;
		for(UINT i2=0;i2<M.Rows[i].Data.Length();i2++)
		{
			E = &M.Rows[i].Data[i2];
			Total += E->Entry * V[E->Col];
		}
		Output[i] = Total;
	}
}

void Multiply(Vector<double> &Output, double D, Vector<double> &V)
{
	for(UINT i=0;i<V.Length();i++)
	{
		Output[i] = V[i] * D;
	}
}

void Multiply(Vector<double> &Output, Vector<double> &V, double D)
{
	Multiply(Output, D, V);
}
double Multiply(Vector<double> &L, Vector<double> &R)
{
	double Total = 0.0;
	for(UINT i=0;i<L.Length();i++)
	{
		Total += L[i] * R[i];
	}
	return Total;
}

void Add(Vector<double> &Output, Vector<double> &L, Vector<double> &R)
{
	for(UINT i=0;i<L.Length();i++)
	{
		Output[i] = L[i] + R[i];
	}
}

void Subtract(Vector<double> &Output, Vector<double> &L, Vector<double> &R)
{
	for(UINT i=0;i<L.Length();i++)
	{
		Output[i] = L[i] - R[i];
	}
}
