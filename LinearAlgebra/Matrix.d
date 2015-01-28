module LinearAlgebra.Matrix;

import std.stdio;
import std.conv;

class Matrix(size_t r, size_t c, T = real)
{
	alias Matrix!(r, c, T) ThisType;

	this(in T[r*c] values...)
	{
		mData[] = values;
	}

	this()
	{
		mData[] = 0;
	}

	ThisType opBinary(string op, rhsType)(Matrix!(r, c, rhsType) rhs)
	{
		static assert(is (T : rhsType));
		static if(op == "+" || op == "-")
		{
			auto res = new ThisType;
			foreach(size_t i, ref element; res.mData)
				mixin("element = mData[i] "~op~" rhs.mData[i];");
			return res;
		}
		else static assert(0, "Operator not implimented");
	}
	
	Matrix!(ir, c, T) opBinaryRight(string op, size_t ir, size_t ic, lhsType)(Matrix!(ir, ic, lhsType) lhs)
	{
		static assert(is (T : lhsType));
		static assert(ic == r, "Incompatible matricies");
		static if(op == "*")
		{
			auto res = new Matrix!(ir, c, T);
			for(int i = 0; i < ir; i++)
				for(int j = 0; j < c; j++)
					for(int k = 0; k < ic; k++)
						res.mData[i*c + j] += lhs.mData[ic*i + k]*mData[k*c + j];

			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	string ToString()
	{
		string matStr;

		for(int i = 0; i < r; i++)
		{
			matStr ~= "[";
			for(int j = 0; j < c; j++)
				matStr ~= " " ~ to!string(mData[i*c + j]);

			matStr ~= " ]\n";
		}
		return matStr;
	}

	override bool opEquals(Object o)
	{
		if(typeid(this) != typeid(o)) return false;

		ThisType rhs = cast(ThisType)o;
		for(int i = 0; i < r*c; i++)
		{
			if(mData[i] != rhs.mData[i])
				return false;
		}
		return true;
	}

	unittest
	{
		auto m1 = new Matrix!(3, 2)(1, 2,
									3, 4,
									5, 6);
		auto m2 = new Matrix!(2, 2)(7, 8,
									9, 10);

		auto m3 = m1 * m2;

		auto expected = new Matrix!(3, 2)(25, 28,
											57, 64,
											89, 100);
			
		assert(m3 == expected);
	}

	unittest
	{
		auto m1 = new Matrix!(3, 2)(1, 2,
			3, 4,
			5, 6);
		
		auto m2 = new Matrix!(3, 2)(2, 3,
			4, 5,
			6, 7);

		auto m3 = m1 + m2;

		auto expected = new Matrix!(3, 2)(3, 5,
			7, 9,
			11, 13);

		assert(m3 == expected);
	}

private:
	size_t mRows = r;
	size_t mCols = c;

	T[r*c] mData;
}

