module numd.linearalgebra.matrix;

import std.conv;
import std.stdio;
import std.math;
import std.meta;
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

			return res;
		}
	}
}

alias Vector(size_t l, T = double) = Matrix!(l, 1, T);

struct Matrix(size_t r, size_t c, T = double)
{
	import std.experimental.allocator.mallocator;

	alias ThisType = Matrix!(r, c, T);

@nogc
{
	// if matrix is small enough, stack allocate it.
	static if(r*c*T.sizeof > 512)
	{
		private int* referenceCount = null;
		T[] mData;
		
		private void allocData()
		{
			mData = cast(T[])Mallocator.instance.allocate(r*c*T.sizeof);
			referenceCount = cast(int*)Mallocator.instance.allocate(int.sizeof);
			(*referenceCount) = 1;
		}
		
		this(const T[] values)
		{
			if(referenceCount is null)
			{
				allocData();
				mData[] = values[];
			}
			else
			{
				mData[] = values[];
			}
		}

		this(this)
		{
			if(referenceCount !is null)
			{
				(*referenceCount)++;
			}
			else
			{
				//printf("Error Matrix referenceCount somehow null in postblit\n");
			}
		}

		~this()
		{
			if(referenceCount !is null)
			{
				if((*referenceCount) == 1)
				{
					if(mData !is null)
					{
						Mallocator.instance.deallocate(mData);
						mData = null;
					}
					Mallocator.instance.deallocate(referenceCount[0..int.sizeof]);
				}
				else if((*referenceCount) > 1)
				{
					(*referenceCount)--;
				}
				else if((*referenceCount) == 0)
				{
					Mallocator.instance.deallocate(referenceCount[0..int.sizeof]);
				}
			}
			else
			{
				//printf("Error Matrix referenceCount somehow null in destructor\n");
			}
		}
	}
	else
	{
		T[r*c] mData;
		
		private void allocData()
		{

		}

		this(const T[] values)
		{
			mData[] = values[];
		}
	}
	
	this(T[r*c] values)
	{
		allocData();
		mData[] = values;
	}

	this(immutable T[r*c] values...)
	{
		allocData();
		mData[] = values;
	}

	this(double init)
	{
		allocData();
		static if(is(T : Mcomplex!double))
		{
			mData[] = Mcomplex!double(Complex!double(init));
		}
		else
		{
			mData[] = init;
		}
	}

	ref T opIndex(size_t index)
	{
		return mData[index];
	}

	ref T opIndex(size_t row, size_t col)
	{
		size_t idx = row*c + col;
		return mData[idx];
	}
	
	T[] opIndex()
	{
		return mData[];
	}

	size_t opDollar(size_t pos)()
	{
		return mData.length;
	}

	T[] opSlice(size_t a, size_t b)
	{
		return mData[a..b];
	}

	Matrix!(r, ic, rhsType) opBinary(string op, size_t ic, rhsType)(ref Matrix!(r, ic, rhsType) rhs)
	{
		static if(op == "+" || op == "-")
		{
			static assert(ic == c, "Incompatible matricies");
			auto res = ThisType(0);
			mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
			return res;
		}
		else static if(op == "*")
		{
			auto res = Matrix!(r, ic, rhsType)(0);
			for(int i = 0; i < r; i++)
				for(int j = 0; j < ic; j++)
					for(int k = 0; k < r; k++)
						res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
			return res;
		}
		else static assert(0, "Operator not implimented");
	}

	Matrix!(r, ic, rhsType) opBinary(string op, size_t ic, rhsType)(Matrix!(r, ic, rhsType) rhs)
	{
		static if(op == "+" || op == "-")
		{
			static assert(ic == c, "Incompatible matricies");
			auto res = ThisType(0);
			mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
			return res;
		}
		else static if(op == "*")
		{
			auto res = Matrix!(r, ic, rhsType)(0);
			for(int i = 0; i < r; i++)
				for(int j = 0; j < ic; j++)
					for(int k = 0; k < r; k++)
						res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
			return res;
		}
		else static assert(0, "Operator not implimented");
	}

	ThisType opBinary(string op)(double rhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = ThisType(0);
			mixin("res.mData[] = mData[]"~op~"rhs;");
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	void mult(size_t ir, size_t ic, lhsType)(ref Matrix!(ir, ic, lhsType) rhs, ref Matrix!(ir, c, lhsType) output)
	{
		static assert(is (T : lhsType));
		static assert(ic == r, "Incompatible matricies");

		for(int i = 0; i < ir; i++)
			for(int j = 0; j < c; j++)
				for(int k = 0; k < ic; k++)
					output.mData[i*c + j] += mData[ic*i + k]*rhs.mData[k*c + j];
	}

	Matrix!(ir, c, T) opBinaryRight(string op, size_t ir, size_t ic, lhsType)(ref Matrix!(ir, ic, lhsType) lhs)
	{
		static assert(is (T : lhsType));
		static assert(ic == r, "Incompatible matricies");
		static if(op == "*")
		{
			auto res = Matrix!(ir, c, T)(0);
			for(int i = 0; i < ir; i++)
				for(int j = 0; j < c; j++)
					for(int k = 0; k < ic; k++)
						res.mData[i*c + j] += lhs.mData[ic*i + k]*mData[k*c + j];

			return res;
		}
		else static if((op == "+") || (op == "-"))
		{
			auto res = Matrix!(ir, c, T)(0);
			mixin("res.mData[] = lhs.mData[] "~op~" mData[];");
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	ThisType opBinaryRight(string op)(double lhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = ThisType(0);
			mixin("res.mData[] = lhs"~op~"mData[];");
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	bool opEquals(ref ThisType o)
	{
		for(int i = 0; i < r*c; i++)
		{
			if(mData[i] != o.mData[i])
				return false;
		}
		return true;
	}

	ThisType rref()
	{
		//auto rref = new ThisType(mData);
		auto rref = ThisType(mData);

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
					//debug writeln("Normalizing row " ~ to!string(currRow));
					//debug writeln(rref.ToString());
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
							//debug writeln("Subtracting row " ~ to!string(currRow) ~ " from " ~ to!string(j));
							//debug writeln(rref.ToString());
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
						//debug writeln("swapping row " ~ to!string(currRow) ~ " with " ~ to!string(j));
						T[c] swapRow;
						swapRow[] = rref.mData[j*c..j*c+c];
						rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
						rref.mData[currRow*c..currRow*c+c] = swapRow[];
						//debug writeln(rref.ToString());
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

	Matrix!(c, r, T) transpose()
	{
		auto ret = Matrix!(c, r, T)(0);
		// Loop over new rows
		for(int i = 0; i < c; i++)
		{
			// Loop over new cols
			for(int j = 0; j < r; j++)
			{
				ret.mData[i*r+j] = mData[j*c + i];
			}
		}
		return ret;
	}

	// Things that only apply to squar matricies.
	static if(c == r)
	{
		static ThisType Identity()
		{
			auto ident = ThisType(0);
			for(int i = 0; i < r*c; i += r+1)
			{
				ident.mData[i] = 1;
			}
			return ident;
		}

		ThisType Inverse()
		{
			auto inv = Identity();
			auto rref = ThisType(mData);

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
							//debug writeln("swapping row " ~ to!string(currRow) ~ " with " ~ to!string(j));
							T[c] swapRow, invSwapRow;
							swapRow[] = rref.mData[j*c..j*c+c];
							invSwapRow[] = inv.mData[j*c..j*c+c];

							rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
							inv.mData[j*c..j*c+c] = inv.mData[currRow*c..currRow*c+c];

							rref.mData[currRow*c..currRow*c+c] = swapRow[];
							inv.mData[currRow*c..currRow*c+c] = invSwapRow[];
							//debug writeln(rref.ToString());
							break;
						}
						// All entries in this column are 0, move on
						else if(j == r-1)
						{
							currCol++;
						}
					}
				}
			}
			return inv;
		}
	}

	static if(r*c*T.sizeof > 512)
	{
		ref ThisType opAssign(ThisType rhs)
		{
			if((referenceCount is null) && (rhs.referenceCount !is null))
			{
				referenceCount = rhs.referenceCount;
				(*referenceCount)++;
				mData = rhs.mData;
			}
			else
			{
				mData[] = rhs.mData[];
			}
			return this;
		}
	}

	ref ThisType opAssign(const T[r*c] rhs)
	{
		mData[] = rhs[];
		return this;
	}

	ref ThisType opAssign(const T[] rhs)
	{
		mData[] = rhs[];
		return this;
	}

	ref ThisType opAssign(const T val)
	{
		mData[] = val;
		return this;
	}
	
	void opOpAssign(string op)(T rhs)
	{
		static assert((op == "*") || (op == "/"), "operator not implimented");
		mixin("mData[] "~op~"= rhs;");	
	}

	void opOpAssign(string op)(ThisType rhs)
	{
		static assert((op == "+") || (op == "-"), "operator not implimented");
		mixin("mData[] "~op~"= rhs.mData[];");
	}

	ThisType opUnary(string s)() if (s == "-")
	{
		auto neg = ThisType(0);
		neg.mData[] = -mData[];
		return neg;
	}

	static if(c == 1)
	{
		static if(r == 3)
		{
			Vector!(3) cross(Vector!(3) rhs)
			{
				auto res = Vector!(3)(0);
				
				res.mData[0] = mData[1]*rhs.mData[2] - mData[2]*rhs.mData[1];
				res.mData[1] = mData[2]*rhs.mData[0] - mData[0]*rhs.mData[2];
				res.mData[2] = mData[0]*rhs.mData[1] - mData[1]*rhs.mData[0];
				
				return res;
			}
		}
		
		T dot(S)(S rhs)
		{
			assert(rows == rhs.rows);
			T res = 0;
			for(size_t i = 0; i < rhs.rows; i++)
				res += mData[i]*rhs.mData[i];
			
			return res;
		}
		
		T magnitude()
		{
			T res = 0;
			foreach(ref element; mData)
				res += element^^2;
			
			res = sqrt(res);
			return res;
		}
		
		ThisType normalize()
		{
			//auto res = new Vector!(l, T);
			auto res = Vector!(r, T)(0);
			res.mData[] = mData[]/magnitude();
			return res;
		}
	}

	T[] getData()
	{
		return mData;
	}

	@property size_t rows() { return r; };
	@property size_t columns() { return c; };
}
	string ToString()
	{
		import std.string;
		string matStr;
		
		for(int i = 0; i < r; i++)
		{
			matStr ~= "[";
			for(int j = 0; j < c; j++)
				matStr ~= " " ~ to!string(mData[i*c + j]).rightJustify(3, ' ');
			
			matStr ~= " ]\n";
		}
		return matStr;
	}
}

static void WriteArrayCSV(size_t ri, size_t ci, T)(ref File f, Matrix!(ri, ci, T) x)
{
	for(int i = 0; i < x.mData.length-1; i++)
		f.writef("%40.40f, ", x.mData[i]);

	f.writef("%40.40f\n", x.mData[$-1]);
}

// opEquals
unittest
{
	auto m1 = Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	
	auto m2 = Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	
	assert(m1 == m2, "Matrix equality test failed");
}

// Identity matrix generation test.
unittest
{
	auto ident = Matrix!(3, 3).Identity();
	auto testIdent = Matrix!(3, 3)(1, 0, 0, 0, 1, 0, 0, 0, 1);

	assert(ident == testIdent, "Matrix identity failed");
}

// op =
unittest
{
	auto m1 = Matrix!(2, 2)(1, 2,
		3, 4);

	auto m2 = Matrix!(2, 2)(0);
	m2 = m1;
	//auto m2 = m1; // This doesn't work. copy ctor not called for some reason
	assert(m1 == m2, "Matrix assignment test failed");

	//writeln(m1.ToString());
	//writeln(m2.ToString());
	// Force failure
	m1.mData[0] = 3;

	//writeln(m1.ToString());
	//writeln(m2.ToString());

	assert(m1 != m2, "Matrix assignment test failed");
}

// opBinaryRight * (non square)
unittest
{
	auto m1 = Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	auto m2 = Matrix!(2, 2)(7, 8,
		9, 10);

	auto m3 = m1 * m2;
	
	auto expected = Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	
	assert(m3 == expected, "Matrix opBinaryRight * test failed");
}

// opBinary * (square)
unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	auto m2 = Matrix!(2, 2)(5, 6, 7, 8);
	auto m3 = Matrix!(2, 2)(0);

	m3 = m1*m2;

	auto expected = Matrix!(2, 2)(19, 22, 43, 50);

	assert(m3 == expected, "Matrix opBinary * test failed");
}

// opBinary * (non-square)
unittest
{
	auto m1 = Matrix!(2, 3)(1, 2, 3, 4, 5, 6);
	auto m2 = Matrix!(3, 2)(5, 6, 7, 8,  9, 10);
	auto m3 = Matrix!(2, 2)(0);
	
	m3 = m1*m2;

	auto expected = Matrix!(2, 2)(46, 52, 109, 124);
	
	assert(m3 == expected, "Matrix opBinary * test failed");
}


// op +
unittest
{
	auto m1 = Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	
	auto m2 = Matrix!(3, 2)(2, 3,
		4, 5,
		6, 7);
	
	auto m3 = m1 + m2;
	
	auto expected = Matrix!(3, 2)(3, 5,
		7, 9,
		11, 13);
	
	assert(m3 == expected, "Matrix addition test failed");
}

// op scalar *
unittest
{
	auto m1 = Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	double scalar = 2;
	
	auto m3 = scalar*m1;
	
	auto expected = Matrix!(3, 2)(2, 4,
		6, 8,
		10, 12);
	
	assert(m3 == expected, "Matrix scalar multiplication test failed");
}

// op scalar /
unittest
{
	auto m1 = Matrix!(3, 2)(2, 4,
		6, 8,
		10, 12);
	
	double scalar = 2;
	
	auto m2 = m1/2;
	
	auto expected = Matrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	
	assert(m2 == expected, "Matrix scalar division test failed");
}

// opUnary -
unittest
{
	auto m1 = Matrix!(3, 3)(1, 2, 3, 4, 5, 6, 7, 8, 9);
	auto m2 = -m1;
	auto expected = Matrix!(3, 3)(-1, -2, -3, -4, -5, -6, -7, -8, -9);

	assert(m2 == expected, "Matrix negative test failed");
}

// op *=
unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	m1 *= 2;
	auto expected = Matrix!(2, 2)(2, 4, 6, 8);

	assert(m1 == expected, "Matrix *= test failed");
}

// op +=
unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	m1 += Matrix!(2, 2).Identity();
	auto expected = Matrix!(2, 2)(2, 2, 3, 5);

	assert(m1 == expected, "Matrix += test failed");
}

// transpose
unittest
{
	auto m1 = Matrix!(3, 3)(1, 2, 3, 4, 5, 6, 7, 8, 9);
	auto m2 = m1.transpose();
	/*
	 1     4     7
     2     5     8
     3     6     9
     */
	auto expected1 = Matrix!(3, 3)(1, 4, 7, 2, 5, 8, 3, 6, 9);

	assert(m2 == expected1, "Matrix transpose test failed");

	auto m3 = Matrix!(2, 3)(1, 2, 3, 4, 5, 6);
	auto m4 = m3.transpose();
	/*
     1     4
     2     5
     3     6
     */
	auto expected2 = Matrix!(3, 2)(1, 4, 2, 5, 3, 6);
	assert(m4 == expected2, "Matrix transpose test failed");
}

// Row reduction
unittest
{

	auto m1 = Matrix!(2, 2)(2, 4, 6, 8);

	auto m2 = m1.rref();

	//writeln(m2.ToString());

	auto expected1 = Matrix!(2, 2).Identity();

	assert(m2 == expected1, "Matrix rref test failed");

	auto m3 = Matrix!(3, 3)(2, 4, 6, 8, 10, 12, 14, 16, 18);

	//writeln(m3.ToString());

	auto m4 = m3.rref();

	auto expected2 = Matrix!(3, 3)(1, 0, -1, 0, 1, 2, 0, 0, 0);

	//writeln(m4.ToString());

	assert(m4 == expected2, "Matrix rref test failed");

	auto m5 = Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m6 = m5.rref();
	auto expected3 = Matrix!(3, 3).Identity();
	//writeln(m6.ToString());
	assert(m6 == expected3, "Matrix rref test failed");

	/*
	 0     4     5
     0     0     7
     0     3     1
	*/
	auto m7 = Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m8 = m7.rref();
	auto expected4 = Matrix!(3, 3)(0, 1, 0, 0, 0, 1, 0, 0, 0);
	//writeln(m8.ToString());
	assert(m8 == expected4, "Matrix rref test failed");


	auto m9 = Matrix!(3, 6)(0, 0, 2, 1, 0, 0, 0, 5, 1, 0, 1, 0, 6, 3, 1, 0, 0, 1);
	auto m10 = m9.rref();
	//writeln(m10.ToString());
	/*
	-0.0333   -0.1000    0.1667
	-0.1000    0.2000        0
	0.5000         0         0
	*/
	auto expected5 = Matrix!(3, 6)(1, 0, 0, -1.0/30.0, -1.0/10.0, 1.0/6.0, 0, 1, 0, -1.0/10.0, 1.0/5.0, 0, 0, 0, 1.0, 1.0/2.0, 0, 0);
	//writeln(expected5.ToString());
	assert(m10 == expected5, "Matrix rref test failed");
}

unittest
{
	auto m1 = Matrix!(3, 3, Mcomplex!double)();
}

unittest
{
	auto m1 = Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m2 = m1.Inverse();
	auto expected = Matrix!(3, 3)(-1.0/30.0, -1.0/10.0, 1.0/6.0, -1.0/10.0, 1.0/5.0, 0, 1.0/2.0, 0, 0);
	//writeln(m2.ToString());
	assert(m2 == expected, "Matrix inverse test failed");

	auto m3 = Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
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

unittest
{
	/+
	auto m1 = Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m2 = Matrix!(3, 3)(0, 1, 2, 0, 5, 1, 4, 3, 1);

	m2.swap(m1);

	assert(m1.mData[1] == 1, "Matrix swap failed");
	assert(m2.mData[1] == 0, "Matrix swap failed");
	assert(m1.mData[6] == 4, "Matrix swap failed");
	assert(m2.mData[6] == 6, "Matrix swap failed");
	+/
}

unittest
{
	auto vec = Vector!(3)(1.0, 1.0, 0.0);
	auto mag = vec.magnitude();
	auto expected = sqrt(cast(double)2);
	
	assert(mag == expected, "Vector magnitude test failed");
}

unittest
{
	auto vec1 = Vector!(3)(1, 1, 0);
	auto vec2 = vec1.normalize();
	auto expected = Vector!(3)(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0);
	
	assert(vec2 == expected, "Vector normalize test failed");
}

unittest
{
	auto vec1 = Vector!(3)(1, 2, 3);
	//writeln(vec1.ToString());
	auto vec2 = Vector!(3)(4, 5, 6);
	//writeln(vec1.ToString());
	auto res = vec1.dot(vec2);
	
	double expected = 32;
	
	assert(res == expected, "Vector dot product test failed");
}

unittest
{
	auto vec1 = Vector!(3)(1, 2, 3);
	auto vec2 = Vector!(3)(4, 5, 6);
	
	auto vec3 = vec1.cross(vec2);
	
	auto expected = Vector!(3)(-3, 6, -3);
	
	assert(vec3 == expected, "Vector cross product failed");
}