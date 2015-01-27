module LinearAlgebra.Matrix;

class Matrix(size_t r, size_t c, T = real)
{
	this(in T[r*c] values...)
	{
		foreach(size_t i, ref value; values)
		{
			mData[i] = value;
		}
	}


private:
	size_t mRows = r;
	size_t mCols = c;

	T[r*c] mData;
}

