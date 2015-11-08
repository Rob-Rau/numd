module LinearAlgebra.Matrix;

import std.conv;
import std.stdio;
import std.complex;
import std.exception;

struct Mcomplex(T)
{
	alias mData this;
	Complex!T mData;
	alias Mcomplex!T ThisType;

	void opAssign(T rhs)
	{
		mData.re = rhs;
	}

	void opAssign(Complex!T rhs)
	{
		mData = rhs;
	}

	ThisType[] opBinary(string op)(ThisType[] rhs)
	{
		static if((op == "*") || (op == "-") || (op == "+") || (op == "/"))
		{
			ThisType[] res;
			foreach(int i, ref ele; rhs)
				mixin("res[i].mData = mData"~op~"ele;");
				//res[i] = mData*ele;

			return res;
		}
	}
}

class Matrix(size_t r, size_t c, T = real)
{
	alias Matrix!(r, c, T) ThisType;

	this(in T[r*c] values...)
	{
		mData[] = values;
	}

	this()
	{
		static if(is(T : Mcomplex!real))
		{
			mData[] = Mcomplex!real(Complex!real(0));
		}
		else
		{
			mData[] = 0;
		}
	}

	this(ThisType mat)
	{
		mData[] = mat.mData[];
	}

	ThisType opBinary(string op, rhsType)(Matrix!(r, c, rhsType) rhs)
	{
		//static assert(is (T : rhsType));
		static if(op == "+" || op == "-")
		{
			auto res = new ThisType;
			foreach(size_t i, ref element; res.mData)
				mixin("element = mData[i] "~op~" rhs.mData[i];");
			return res;
		}
		else static assert(0, "Operator not implimented");
	}

	ThisType opBinary(string op)(real rhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = new ThisType;
			foreach(size_t i, ref element; mData)
				mixin("res.mData[i] = element"~op~"rhs;");

			return res;
		}
		else static assert(0, "Operator not implemented");
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

	ThisType opBinaryRight(string op)(real lhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = new ThisType;
			foreach(size_t i, element; mData)
				mixin("res.mData[i] = lhs"~op~"element;");
			
			return res;
		}
		else static assert(0, "Operator not implemented");
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

	ThisType rref()
	{
		auto rref = new ThisType(mData);

		int currRow = 0, currCol = 0;
		while((currRow < r) && (currCol < c))
		{
			if(rref.mData[currRow*c + currCol] != 0)
			{
				T[] row;
				row = rref.mData[currRow*c..currRow*c+c];

				// Normalize row so first entry is 1
				if(rref.mData[currRow*c + currCol] != 1)
				{
					auto first = rref.mData[currRow*c + currCol];
					row[] /= first;
					debug writeln("Normalizing row " ~ to!string(currRow));
					debug writeln(rref.ToString());
				}
				// Find all non-zero rows in this column and scale the current 
				// row by the value in other row, then subtract current row from
				// from other
				for(int j = 0; j < r; j++)
				{
					if(j != currRow)
					{
						T[] nextRow;
						nextRow = rref.mData[j*c..j*c+c];
						if(rref.mData[j*c + currCol] != 0)
						{
							auto first = rref.mData[j*c + currCol];
							nextRow[] -= (first*row[])[];
							debug writeln("Subtracting row " ~ to!string(currRow) ~ " from " ~ to!string(j));
							debug writeln(rref.ToString());
						}
					}
				}
				currRow++;
				currCol++;
			}
			else if((currRow == r-1) && (currCol == c - 1))
			{
				currCol++;
			}
			// See if there is a row that has a non-zero entry further down
			// if so swap it with current row
			else
			{
				for(int j = (currRow+1); j < r; j++)
				{
					// Swap this row with current row
					if(rref.mData[j*c + currCol] != 0)
					{
						debug writeln("swapping row " ~ to!string(currRow) ~ " with " ~ to!string(j));
						T[c] swapRow;
						swapRow[] = rref.mData[j*c..j*c+c];
						rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
						rref.mData[currRow*c..currRow*c+c] = swapRow[];
						debug writeln(rref.ToString());
						break;
					}
					// All entries in this column are 0, move on
					else if(j == r-1)
					{
						//currRow++;
						currCol++;
					}
				}
				//currCol++;
			}
		}
		return rref;
	}
	// Things that only apply to squar matricies.
	static if(c == r)
	{
		static ThisType Identity()
		{
			//assert(r == c, "Not a square matrix, no identity");
			auto ident = new ThisType;
			//ident.mData[] = 0;
			for(int i = 0; i < r*c; i += r+1)
			{
				ident.mData[i] = 1;
			}
			return ident;
		}

		ThisType Inverse()
		{
			auto inv = Identity();
			auto rref = new ThisType(mData);

			int currRow = 0, currCol = 0;
			while((currRow < r) && (currCol < c))
			{
				if(rref.mData[currRow*c + currCol] != 0)
				{
					T[] row, invRow;
					row = rref.mData[currRow*c..currRow*c+c];
					invRow = inv.mData[currRow*c..currRow*c+c];
					// Normalize row so first entry is 1
					if(rref.mData[currRow*c + currCol] != 1)
					{
						auto first = rref.mData[currRow*c + currCol];
						row[] /= first;
						invRow[] /= first;
						debug writeln("Normalizing row " ~ to!string(currRow));
						debug writeln(rref.ToString());
					}
					// Find all non-zero rows in this column and scale the current 
					// row by the value in other row, then subtract current row from
					// from other
					for(int j = 0; j < r; j++)
					{
						if(j != currRow)
						{
							T[] nextRow, invNexRow;
							nextRow = rref.mData[j*c..j*c+c];
							invNexRow = inv.mData[j*c..j*c+c];

							if(rref.mData[j*c + currCol] != 0)
							{
								auto first = rref.mData[j*c + currCol];
								nextRow[] -= (first*row[])[];
								invNexRow[] -= (first*invRow[])[];
								debug writeln("Subtracting row " ~ to!string(currRow) ~ " from " ~ to!string(j));
								debug writeln(rref.ToString());
							}
						}
					}
					currRow++;
					currCol++;
				}
				else if((currRow == r-1) && (currCol == c - 1))
				{
					currCol++;
				}
				// See if there is a row that has a non-zero entry further down
				// if so swap it with current row
				else
				{
					for(int j = (currRow+1); j < r; j++)
					{
						// Swap this row with current row
						if(rref.mData[j*c + currCol] != 0)
						{
							debug writeln("swapping row " ~ to!string(currRow) ~ " with " ~ to!string(j));
							T[c] swapRow, invSwapRow;
							swapRow[] = rref.mData[j*c..j*c+c];
							invSwapRow[] = inv.mData[j*c..j*c+c];

							rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
							inv.mData[j*c..j*c+c] = inv.mData[currRow*c..currRow*c+c];

							rref.mData[currRow*c..currRow*c+c] = swapRow[];
							inv.mData[currRow*c..currRow*c+c] = invSwapRow[];
							debug writeln(rref.ToString());
							break;
						}
						// All entries in this column are 0, move on
						else if(j == r-1)
						{
							//currRow++;
							currCol++;
						}
					}
					//currCol++;
				}
			}

			enforce(rref == ThisType.Identity(), "Matrix not invertible");

			return inv;
		}
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

	@property size_t rows() { return mRows; };
	@property size_t columns() { return mCols; };

package:
	size_t mRows = r;
	size_t mCols = c;

	T[r*c] mData;
}

// opEquals
unittest
{
	auto m1 = new Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	
	auto m2 = new Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	
	assert(m1 == m2, "Matrix equality test failed");
}

// Identity matrix generation test.
unittest
{
	auto ident = Matrix!(3, 3).Identity();
	auto testIdent = new Matrix!(3, 3)(1, 0, 0, 0, 1, 0, 0, 0, 1);

	assert(ident == testIdent, "Matrix identity failed");
}

// op =
unittest
{
	auto m1 = new Matrix!(2, 2)(1, 2,
		3, 4);

	auto m2 = new Matrix!(2,2)(m1);
	
	assert(m1 == m2, "Matrix assignment test failed");

	// Force failure
	m1.mData[0] = 3;
	assert(m1 != m2, "Matrix assignment test failed");
}

// op *
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
	
	assert(m3 == expected, "Matrix multiplication test failed");
}

// op +
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
	
	assert(m3 == expected, "Matrix addition test failed");
}

// op scalar *
unittest
{
	auto m1 = new Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	real scalar = 2;
	
	auto m3 = scalar*m1;
	
	auto expected = new Matrix!(3, 2)(2, 4,
		6, 8,
		10, 12);
	
	assert(m3 == expected, "Matrix scalar multiplication test failed");
}

// op scalar /
unittest
{
	auto m1 = new Matrix!(3, 2)(2, 4,
		6, 8,
		10, 12);
	
	real scalar = 2;
	
	auto m2 = m1/2;
	
	auto expected = new Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	
	assert(m2 == expected, "Matrix scalar division test failed");
}

unittest
{

	auto m1 = new Matrix!(2, 2)(2, 4, 6, 8);

	auto m2 = m1.rref();

	writeln(m2.ToString());

	auto expected1 = Matrix!(2, 2).Identity();

	assert(m2 == expected1, "Matrix rref test failed");

	auto m3 = new Matrix!(3, 3)(2, 4, 6, 8, 10, 12, 14, 16, 18);

	//writeln(m3.ToString());

	auto m4 = m3.rref();

	auto expected2 = new Matrix!(3, 3)(1, 0, -1, 0, 1, 2, 0, 0, 0);

	writeln(m4.ToString());

	assert(m4 == expected2, "Matrix rref test failed");

	auto m5 = new Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m6 = m5.rref();
	auto expected3 = Matrix!(3, 3).Identity();
	writeln(m6.ToString());
	assert(m6 == expected3, "Matrix rref test failed");

	/*
	 0     4     5
     0     0     7
     0     3     1
	*/
	auto m7 = new Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m8 = m7.rref();
	auto expected4 = new Matrix!(3, 3)(0, 1, 0, 0, 0, 1, 0, 0, 0);
	writeln(m8.ToString());
	assert(m8 == expected4, "Matrix rref test failed");


	auto m9 = new Matrix!(3, 6)(0, 0, 2, 1, 0, 0, 0, 5, 1, 0, 1, 0, 6, 3, 1, 0, 0, 1);
	auto m10 = m9.rref();
	writeln(m10.ToString());
	/*
	-0.0333   -0.1000    0.1667
	-0.1000    0.2000        0
	0.5000         0         0
	*/
	auto expected5 = new Matrix!(3, 6, real)(1, 0, 0, -1.0/30.0, -1.0/10.0, 1.0/6.0, 0, 1, 0, -1.0/10.0, 1.0/5.0, 0, 0, 0, 1.0, 1.0/2.0, 0, 0);
	//writeln(expected5.ToString());
	assert(m10 == expected5, "Matrix rref test failed");
}

unittest
{
	auto m1 = new Matrix!(3, 3, Mcomplex!real)();
}

unittest
{
	auto m1 = new Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m2 = m1.Inverse();
	auto expected = new Matrix!(3, 3)(-1.0/30.0, -1.0/10.0, 1.0/6.0, -1.0/10.0, 1.0/5.0, 0, 1.0/2.0, 0, 0);
	writeln(m2.ToString());
	assert(m2 == expected, "Matrix inverse test failed");

	auto m3 = new Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	Matrix!(3,3) m4;
	try
	{
		m4 = m3.Inverse();
		assert(false, "Matrix inverse test failed");
	}
	catch(Exception e)
	{
		writeln(e.msg);
	}
}