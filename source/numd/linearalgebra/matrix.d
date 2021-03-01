module numd.linearalgebra.matrix;

import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.range;
import std.stdio;
import std.traits;
import std.typecons;

/+
alias RMatrix(size_t r, size_t c, T = double) = RefCounted!(Matrix!(r, c, T));
alias RVector(size_t l, T = double) = RefCounted!(Vector!(l, T));

// opEquals
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	auto m2 = RMatrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	assert(m1 == m2, "Matrix equality test failed");
}
// identity matrix generation test.
@nogc @system unittest
{
	auto ident = RMatrix!(3, 3).identity();
	auto testIdent = RMatrix!(3, 3)(1, 0, 0, 0, 1, 0, 0, 0, 1);
	assert(ident == testIdent, "RMatrix identity failed");
}

// opBinaryRight * (non square)
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	auto m2 = RMatrix!(2, 2)(7, 8,
		9, 10);
	auto m3 = m1 * m2;
	auto expected = RMatrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	assert(m3 == expected, "RMatrix opBinaryRight * test failed");
}

// opBinary * (square)
@nogc @system unittest
{
	auto m1 = RMatrix!(2, 2)(1, 2, 3, 4);
	auto m2 = RMatrix!(2, 2)(5, 6, 7, 8);
	auto m3 = RMatrix!(2, 2)(0);
	m3 = m1*m2;
	auto expected = RMatrix!(2, 2)(19, 22, 43, 50);
	assert(m3 == expected, "RMatrix opBinary * test failed");
}

// opBinary * (non-square)
@nogc @system unittest
{
	auto m1 = RMatrix!(2, 3)(1, 2, 3, 4, 5, 6);
	auto m2 = RMatrix!(3, 2)(5, 6, 7, 8,  9, 10);
	auto m3 = RMatrix!(2, 2)(0);
	m3 = m1*m2;
	auto expected = RMatrix!(2, 2)(46, 52, 109, 124);
	assert(m3 == expected, "RMatrix opBinary * test failed");
}


// op +
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	auto m2 = RMatrix!(3, 2)(2, 3,
		4, 5,
		6, 7);
	auto m3 = m1 + m2;
	auto expected = RMatrix!(3, 2)(3, 5,
		7, 9,
		11, 13);
	assert(m3 == expected, "RMatrix addition test failed");
}

// op scalar *
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 2)(1, 2,
		3, 4,
		5, 6);
	double scalar = 2;
	auto m3 = scalar*m1;
	auto expected = RMatrix!(3, 2)(2, 4,
		6, 8,
		10, 12);
	assert(m3 == expected, "RMatrix scalar multiplication test failed");
}

// op scalar /
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 2)(2, 4,
							6, 8,
							10, 12);
	double scalar = 2;
	auto m2 = m1/scalar;
	auto expected = RMatrix!(3, 2)(1, 2,
								  3, 4,
								  5, 6);
	assert(m2 == expected, "RMatrix scalar division test failed");
}

// opUnary -
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 3)(1, 2, 3, 4, 5, 6, 7, 8, 9);
	auto m2 = -m1;
	auto expected = RMatrix!(3, 3)(-1, -2, -3, -4, -5, -6, -7, -8, -9);
	assert(m2 == expected, "RMatrix negative test failed");
}

// op *=
@nogc @system unittest
{
	auto m1 = RMatrix!(2, 2)(1, 2, 3, 4);
	m1 *= 2;
	auto expected = RMatrix!(2, 2)(2, 4, 6, 8);
	assert(m1 == expected, "RMatrix *= test failed");
}

// op +=
@nogc @system unittest
{
	auto m1 = RMatrix!(2, 2)(1, 2, 3, 4);
	m1 += RMatrix!(2, 2).identity();
	auto expected = RMatrix!(2, 2)(2, 2, 3, 5);
	assert(m1 == expected, "RMatrix += test failed");
}

// transpose
@nogc @system unittest
{
	auto m1 = RMatrix!(3, 3)(1, 2, 3, 4, 5, 6, 7, 8, 9);
	auto m2 = m1.transpose();
	/*
	 1     4     7
     2     5     8
     3     6     9
     */
	auto expected1 = RMatrix!(3, 3)(1, 4, 7, 2, 5, 8, 3, 6, 9);

	assert(m2 == expected1, "RMatrix transpose test failed");

	auto m3 = RMatrix!(2, 3)(1, 2, 3, 4, 5, 6);
	auto m4 = m3.transpose();
	/*
     1     4
     2     5
     3     6
     */
	auto expected2 = RMatrix!(3, 2)(1, 4, 2, 5, 3, 6);
	assert(m4 == expected2, "RMatrix transpose test failed");
}

// Row reduction
@nogc @system unittest
{
	auto m1 = RMatrix!(2, 2)(2, 4, 6, 8);
	auto m2 = m1.rref();
	auto expected1 = RMatrix!(2, 2).identity();
	assert(m2 == expected1, "RMatrix rref test failed");

	auto m3 = RMatrix!(3, 3)(2, 4, 6, 8, 10, 12, 14, 16, 18);
	auto m4 = m3.rref();
	auto expected2 = RMatrix!(3, 3)(1, 0, -1, 0, 1, 2, 0, 0, 0);
	assert(m4 == expected2, "RMatrix rref test failed");

	auto m5 = RMatrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m6 = m5.rref();
	auto expected3 = RMatrix!(3, 3).identity();
	assert(m6 == expected3, "RMatrix rref test failed");

	/*
	 0     4     5
     0     0     7
     0     3     1
	*/
	auto m7 = RMatrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m8 = m7.rref();
	auto expected4 = RMatrix!(3, 3)(0, 1, 0, 0, 0, 1, 0, 0, 0);
	assert(m8 == expected4, "RMatrix rref test failed");


	auto m9 = RMatrix!(3, 6)(0, 0, 2, 1, 0, 0, 0, 5, 1, 0, 1, 0, 6, 3, 1, 0, 0, 1);
	auto m10 = m9.rref();
	/*
	-0.0333   -0.1000    0.1667
	-0.1000    0.2000        0
	0.5000         0         0
	*/
	auto expected5 = RMatrix!(3, 6)(1, 0, 0, -1.0/30.0, -1.0/10.0, 1.0/6.0, 0, 1, 0, -1.0/10.0, 1.0/5.0, 0, 0, 0, 1.0, 1.0/2.0, 0, 0);
	assert(m10 == expected5, "RMatrix rref test failed");
}

@nogc @system unittest
{
	auto m1 = RMatrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m2 = m1.inverse();
	assert(!m2.isNull, "Matrix inversion failed, should not have returned null");
	auto expected = RMatrix!(3, 3)(-1.0/30.0, -1.0/10.0, 1.0/6.0, -1.0/10.0, 1.0/5.0, 0, 1.0/2.0, 0, 0);
	assert(m2.get == expected, "RMatrix inverse test failed");

	auto m3 = RMatrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m4 = m3.inverse();
	assert(m4.isNull, "RMatrix inverse test failed, should have been empty Nullable");
}

@nogc @system unittest
{
	auto vec = RVector!(3)(1.0, 1.0, 0.0);
	auto mag = vec.magnitude();
	auto expected = sqrt(cast(double)2);	
	assert(mag == expected, "RVector magnitude test failed");
}

@nogc @system unittest
{
	auto vec1 = RVector!(3)(1, 1, 0);
	auto vec2 = vec1.normalize();
	auto expected = RVector!(3)(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0);
	assert(vec2 == expected, "RVector normalize test failed");
}

@nogc @system unittest
{
	auto vec1 = RVector!(3)(1, 2, 3);
	auto vec2 = RVector!(3)(4, 5, 6);
	auto res = vec1.dot(vec2);
	double expected = 32;
	assert(res == expected, "RVector dot product test failed");
}

@nogc @system unittest
{
	auto vec1 = RVector!(3)(1, 2, 3);
	auto vec2 = RVector!(3)(4, 5, 6);
	auto vec3 = vec1.cross(vec2);
	auto expected = RVector!(3)(-3, 6, -3);
	assert(vec3 == expected, "RVector cross product failed");
}
+/

@safe @nogc unittest {
	immutable pos1 = Vector!(3, size_t)(0, 1, 2);
	immutable pos2 = Vector!(3, size_t)(2, 1, 0);
	immutable pos3 = Vector!(3, size_t)(0, 1, 2);

	assert(pos1.toHash != pos2.toHash);
	assert(pos1.toHash == pos3.toHash);

	assert(pos1 != pos2);
	assert(pos1 == pos3);
}

alias Vector(size_t l, T = double) = Matrix!(l, 1, T);

struct Matrix(size_t r, size_t c, T = double)
{
	alias ThisType = typeof(this);

@nogc
{
	T[r*c] mData;

	this(const T[] values)
	{
		mData[] = values[];
	}
	
	this(T[r*c] values)
	{
		mData[] = values;
	}

	this(immutable T[r*c] values...)
	{
		mData[] = values;
	}

	this(T init)
	{
		static if(isArray!T) {
			foreach(idx; 0..r*c) {
				mData[idx][] = init[];
			}
		} else {
			mData[] = init;
		}
	}

	static if(isArray!T) {
		this(ForeachType!T init)
		{
			foreach(idx; 0..r*c) {
				mData[idx][] = init;
			}
		}
	}

	@trusted bool opEquals(ref const Matrix!(r, c, T) rhs) const {
		bool equal = true;
		foreach(idx; 0..r*c) {
			equal = equal && (mData[idx] == rhs.mData[idx]);
		}
		return equal;
	}

	@trusted size_t toHash() const nothrow {
		return mData.hashOf;
	}

	@trusted ref T opIndex(size_t index)
	{
		assert(index < r*c);
		return mData[index];
	}

	@trusted T opIndex(size_t index) immutable
	{
		assert(index < r*c);
		return mData[index];
	}

	@trusted T opIndex(size_t index) const
	{
		assert(index < r*c);
		return mData[index];
	}

	@trusted void opIndexAssign(T element, size_t index)
	{
		static if(isArray!T) {
			mData[index][] = element[];
		} else {
			mData[index] = element;
		}
	}

	@trusted void opIndexOpAssign(string op)(T element, size_t index)
	{
		assert(index < r*c);
		static if(isArray!T) {
			mixin("mData[index][] "~op~"= element[];");
		} else {
			mixin("mData[index] "~op~"= element;");
		}

		
	}

	@trusted ref T opIndex(size_t row, size_t col)
	{
		assert(row < r);
		assert(col < c);

		size_t idx = row*c + col;
		return mData[idx];
	}
	
	@trusted T opIndex(size_t row, size_t col) immutable
	{
		assert(row < r);
		assert(col < c);

		size_t idx = row*c + col;
		return mData[idx];
	}

	@trusted T opIndex(size_t row, size_t col) const
	{
		assert(row < r);
		assert(col < c);

		size_t idx = row*c + col;
		return mData[idx];
	}
	

	@trusted void opIndexAssign(T element, size_t row, size_t col)
	{
		assert(row < r);
		assert(col < c);

		size_t idx = row*c + col;
		static if(isArray!T) {
			mData[idx][] = element[];
		} else {
			mData[idx] = element;
		}
		
	}

	@trusted void opIndexOpAssign(string op)(T element, size_t row, size_t col)
	{
		assert(row < r);
		assert(col < c);

		size_t idx = row*c + col;
		static if(isArray!T) {
			mixin("mData[idx][] "~op~"= element[];");
		} else {
			mixin("mData[idx] "~op~"= element;");
		}
		
	}

	@trusted T[] opIndex()
	{
		return mData[];
	}

	@trusted size_t opDollar(size_t pos)()
	{
		return mData.length;
	}

	@trusted T[] opSlice(size_t a, size_t b)
	{
		assert(a < r*c);
		assert(b <= r*c);

		return mData[a..b];
	}

	static if(isArray!T) {
		@trusted Matrix!(r, ic, T) opBinary(string op, size_t ir, size_t ic, rhsType)(ref Matrix!(ir, ic, rhsType) rhs)
		{
			pragma(inline, true);
			static if(op == "+" || op == "-")
			{
				static assert(ic == c, "Incompatible matricies");
				static assert(ir == r, "Incompatible matricies");
				auto res = ThisType(0);
				static if(isArray!T && isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx][]"~op~"rhs.mData[idx][];");
					}
				} else static if(isArray!T && !isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx][]"~op~"rhs.mData[idx];");
					}
				} else {
					mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
				}
				return res;
			}
			else static if(op == "*")
			{

				auto res = Matrix!(r, ic, rhsType)(0);
				static assert(c == ir, "Incompatible matricies for multiplication");
				
				foreach(i; 0..r) {
					foreach(j; 0..ic) {
						foreach(k; 0..ir) {
							static if(isArray!T && isArray!rhsType) {
								T tmp = mData[c*i + k][]*rhs.mData[k*ic + j][];
								T tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else static if(isArray!T && !isArray!rhsType) {
								T tmp = mData[c*i + k][]*rhs.mData[k*ic + j];
								T tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else {
								res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
							}
						}
					}
				}

				return res;
			}
			else static assert(0, "Operator not implimented");
		}

		@trusted Matrix!(r, ic, T) opBinary(string op, size_t ir, size_t ic, rhsType)(Matrix!(ir, ic, rhsType) rhs)
		{
			pragma(inline, true);
			static if(op == "+" || op == "-")
			{
				static assert(ic == c, "Incompatible matricies");
				static assert(ir == r, "Incompatible matricies");
				auto res = ThisType(0);
				static if(isArray!T && isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx][]"~op~"rhs.mData[idx][];");
					}
				} else static if(isArray!T && !isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx][]"~op~"rhs.mData[idx];");
					}
				} else {
					mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
				}
				return res;
			}
			else static if(op == "*")
			{

				auto res = Matrix!(r, ic, rhsType)(0);
				static assert(c == ir, "Incompatible matricies for multiplication");
				
				foreach(i; 0..r) {
					foreach(j; 0..ic) {
						foreach(k; 0..ir) {
							static if(isArray!T && isArray!rhsType) {
								T tmp = mData[c*i + k][]*rhs.mData[k*ic + j][];
								T tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else static if(isArray!T && !isArray!rhsType) {
								T tmp = mData[c*i + k][]*rhs.mData[k*ic + j];
								T tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else {
								res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
							}
						}
					}
				}

				return res;
			}
			else static assert(0, "Operator not implimented");
		}
	}
	static if(!isArray!T) {
		@trusted auto opBinary(string op, size_t ir, size_t ic, rhsType)(immutable Matrix!(ir, ic, rhsType) rhs) immutable
		{
			pragma(inline, true);
			static if(op == "+" || op == "-")
			{
				static assert(ic == c, "Incompatible matricies");
				static assert(ir == r, "Incompatible matricies");
				
				//auto res = Matrix!(r, ic, rhsType)(0);
				static if(isArray!rhsType) {
					static assert(false);
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx]"~op~"rhs.mData[idx][];");
					}
				} else {
					mixin("rhsType[r*ic] ret_data = (cast(immutable(rhsType)[])mData)[]"~op~"rhs.mData[];");
				}
				return Matrix!(r, ic, rhsType)(ret_data);
			}
			else static if(op == "*")
			{

				auto res = Matrix!(r, ic, rhsType)(0);
				static assert(c == ir, "Incompatible matricies for multiplication");
				
				foreach(i; 0..r) {
					foreach(j; 0..ic) {
						foreach(k; 0..ir) {
							static if(isArray!rhsType) {
								rhsType tmp = mData[c*i + k]*rhs.mData[k*ic + j][];
								rhsType tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else {
								res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
							}
						}
					}
				}

				return res;
			}
			else static assert(0, "Operator not implimented");
		}

		@trusted Matrix!(r, ic, rhsType) opBinary(string op, size_t ir, size_t ic, rhsType)(ref Matrix!(ir, ic, rhsType) rhs)
		{
			pragma(inline, true);
			static if(op == "+" || op == "-")
			{
				static assert(ic == c, "Incompatible matricies");
				static assert(ir == r, "Incompatible matricies");
				auto res = Matrix!(r, ic, rhsType)(0);
				static if(isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx]"~op~"rhs.mData[idx][];");
					}
				} else {
					mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
				}
				return res;
			}
			else static if(op == "*")
			{

				auto res = Matrix!(r, ic, rhsType)(0);
				static assert(c == ir, "Incompatible matricies for multiplication");
				
				foreach(i; 0..r) {
					foreach(j; 0..ic) {
						foreach(k; 0..ir) {
							static if(isArray!rhsType) {
								rhsType tmp = mData[c*i + k]*rhs.mData[k*ic + j][];
								rhsType tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else {
								res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
							}
						}
					}
				}

				return res;
			}
			else static assert(0, "Operator not implimented");
		}

		@trusted Matrix!(r, ic, T) opBinary(string op, size_t ir, size_t ic, rhsType)(Matrix!(ir, ic, rhsType) rhs)
		{
			pragma(inline, true);
			static if(op == "+" || op == "-")
			{
				static assert(ic == c, "Incompatible matricies");
				static assert(ir == r, "Incompatible matricies");
				auto res = Matrix!(r, ic, rhsType)(0);
				static if(isArray!rhsType) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx]"~op~"rhs.mData[idx][];");
					}
				} else {
					mixin("res.mData[] = mData[]"~op~"rhs.mData[];");
				}
				return res;
			}
			else static if(op == "*")
			{

				auto res = Matrix!(r, ic, rhsType)(0);
				static assert(c == ir, "Incompatible matricies for multiplication");
				
				foreach(i; 0..r) {
					foreach(j; 0..ic) {
						foreach(k; 0..ir) {
							static if(isArray!rhsType) {
								rhsType tmp = mData[c*i + k]*rhs.mData[k*ic + j][];
								rhsType tmp1 = res.mData[i*ic + j][] + tmp[];
								res.mData[i*ic + j][] = tmp1[];
							} else {
								res.mData[i*ic + j] += mData[c*i + k]*rhs.mData[k*ic + j];
							}
						}
					}
				}

				return res;
			}
			else static assert(0, "Operator not implimented");
		}
	}

	@trusted inout ThisType opBinary(string op)(T rhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = ThisType(0);
			static if(isArray!T) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = mData[idx][]"~op~"rhs[];");
				}
			} else {
				mixin("res.mData[] = mData[]"~op~"rhs;");
			}
			
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	static if(isArray!T) {
		@trusted inout ThisType opBinary(string op)(ForeachType!T rhs)
		{
			static if(op == "*" || op == "/")
			{
				auto res = ThisType(0);
				static if(isArray!T) {
					foreach(idx; 0..r*c) {
						mixin("res.mData[idx][] = mData[idx][]"~op~"rhs;");
					}
				} else {
					mixin("res.mData[] = mData[]"~op~"rhs;");
				}
				
				return res;
			}
			else static assert(0, "Operator not implemented");
		}
	}

	@trusted inout Matrix!(ir, c, T) opBinaryRight(string op, size_t ir, size_t ic, lhsType)(ref inout Matrix!(ir, ic, lhsType) lhs)
	{
		//static assert(is (T : lhsType));
		static if(op == "*")
		{
			static assert(ic == r, "Incompatible matricies. ic == "~ic.to!string~", r == "~r.to!string);
			auto res = Matrix!(ir, c, T)(0);
			for(int i = 0; i < ir; i++) {
				for(int j = 0; j < c; j++) {
					for(int k = 0; k < ic; k++) {
						static if(isArray!T && isArray!lhsType) {
							res.mData[i*c + j][] += lhs.mData[ic*i + k][]*mData[k*c + j][];
						} else static if(isArray!T && !isArray!lhsType) {
							res.mData[i*c + j][] += lhs.mData[ic*i + k]*mData[k*c + j][];
						} else {
							res.mData[i*c + j] += lhs.mData[ic*i + k]*mData[k*c + j];
						}
					}
				}
			}

			return res;
		}
		else static if((op == "+") || (op == "-"))
		{
			static assert((ic == c) && (ir == r), "Incompatible matricies. ic == "~ic.to!string~", r == "~r.to!string);
			auto res = Matrix!(ir, c, T)(0);
			static if(isArray!T && isArray!lhsType) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs.mData[idx][] "~op~" mData[idx][];");
				}
			} static if(isArray!T && !isArray!lhsType) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs.mData[idx] "~op~" mData[idx][];");
				}
			} else {
				mixin("res.mData[] = lhs.mData[] "~op~" mData[];");
			}
			

			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	@trusted auto opBinaryRight(string op, size_t ir, size_t ic, lhsType)(immutable Matrix!(ir, ic, lhsType) lhs) immutable
	{
		//static assert(is (T : lhsType));
		static if(op == "*")
		{
			static assert(ic == r, "Incompatible matricies. ic == "~ic.to!string~", r == "~r.to!string);
			auto res = Matrix!(ir, c, T)(0);
			for(int i = 0; i < ir; i++) {
				for(int j = 0; j < c; j++) {
					for(int k = 0; k < ic; k++) {
						static if(isArray!T && isArray!lhsType) {
							res.mData[i*c + j][] += lhs.mData[ic*i + k][]*mData[k*c + j][];
						} else static if(isArray!T && !isArray!lhsType) {
							res.mData[i*c + j][] += lhs.mData[ic*i + k]*mData[k*c + j][];
						} else {
							res.mData[i*c + j] += lhs.mData[ic*i + k]*mData[k*c + j];
						}
					}
				}
			}

			return res;
		}
		else static if((op == "+") || (op == "-"))
		{
			static assert((ic == c) && (ir == r), "Incompatible matricies. ic == "~ic.to!string~", r == "~r.to!string);
			auto res = Matrix!(ir, c, T)(0);
			static if(isArray!T && isArray!lhsType) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs.mData[idx][] "~op~" mData[idx][];");
				}
			} static if(isArray!T && !isArray!lhsType) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs.mData[idx] "~op~" mData[idx][];");
				}
			} else {
				mixin("res.mData[] = lhs.mData[] "~op~" mData[];");
			}
			

			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	@trusted inout ThisType opBinaryRight(string op)(T lhs)
	{
		static if(op == "*" || op == "/")
		{
			auto res = ThisType(0);
			static if(isArray!T) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs[]"~op~"mData[idx][];");
				}
			} else {
				mixin("res.mData[] = lhs"~op~"mData[];");
			}
			
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	@trusted auto opBinaryRight(string op, _T)(_T lhs) immutable
	{
		static if(op == "*" || op == "/")
		{
			alias RT = typeof(mixin("mData[0]"~op~"lhs"));
			alias RetType = Matrix!(r, c, RT);
			auto res = RetType(0);
			static if(isArray!T) {
				foreach(idx; 0..r*c) {
					mixin("res.mData[idx][] = lhs[]"~op~"mData[idx][];");
				}
			} else {
				mixin("res.mData[] = lhs"~op~"mData[];");
			}
			
			return res;
		}
		else static assert(0, "Operator not implemented");
	}

	@trusted bool opEquals(ref ThisType o)
	{
		for(int i = 0; i < r*c; i++)
		{
			static if(isArray!T) {
				foreach(sub_idx; 0..mData[i].length) {
					if(abs(mData[i][sub_idx] - o.mData[i][sub_idx]) > fpTol)
						return false;
				}
			} else {
				if(abs(mData[i] - o.mData[i]) > fpTol)
					return false;
			}
		}
		return true;
	}

	@trusted bool opEquals(immutable ThisType o) immutable
	{
		for(int i = 0; i < r*c; i++)
		{
			static if(isArray!T) {
				foreach(sub_idx; 0..mData[i].length) {
					if(abs(mData[i][sub_idx] - o.mData[i][sub_idx]) > fpTol)
						return false;
				}
			} else {
				if(abs(mData[i] - o.mData[i]) > fpTol)
					return false;
			}
		}
		return true;
	}

	static if(!isArray!T) {
		@trusted ThisType rref()
		{
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
							T[c] swapRow;
							swapRow[] = rref.mData[j*c..j*c+c];
							rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
							rref.mData[currRow*c..currRow*c+c] = swapRow[];
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
			return rref;
		}
	}

	@trusted inout Matrix!(c, r, T) transpose()
	{
		auto ret = Matrix!(c, r, T)(0);
		// Loop over new rows
		for(int i = 0; i < c; i++)
		{
			// Loop over new cols
			for(int j = 0; j < r; j++)
			{
				static if(isArray!T) {
					ret.mData[i*r+j][] = mData[j*c + i][];
				} else {
					ret.mData[i*r+j] = mData[j*c + i];
				}
				
			}
		}
		return ret;
	}

	// Things that only apply to squar matricies.
	static if(c == r)
	{
		@trusted static ThisType identity()
		{
			auto ident = ThisType(0);
			for(int i = 0; i < r*c; i += r+1)
			{
				ident.mData[i] = 1;
			}
			return ident;
		}

		static if(!isArray!T) {
			@trusted Nullable!ThisType inverse()
			{
				auto inv = identity();
				auto rref = ThisType(mData);

				int currRow = 0, currCol = 0;
				while((currRow < r) && (currCol < c))
				{
					if(abs(rref.mData[currRow*c + currCol]) > fpTol)
					{
						T[] row, invRow;
						row = rref.mData[currRow*c..currRow*c+c];
						invRow = inv.mData[currRow*c..currRow*c+c];
						// Normalize row so first entry is 1
						if(abs(rref.mData[currRow*c + currCol] - 1) > fpTol)
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

								if(abs(rref.mData[j*c + currCol]) > fpTol)
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
							if(abs(rref.mData[j*c + currCol]) > fpTol)
							{
								T[c] swapRow, invSwapRow;
								swapRow[] = rref.mData[j*c..j*c+c];
								invSwapRow[] = inv.mData[j*c..j*c+c];

								rref.mData[j*c..j*c+c] = rref.mData[currRow*c..currRow*c+c];
								inv.mData[j*c..j*c+c] = inv.mData[currRow*c..currRow*c+c];

								rref.mData[currRow*c..currRow*c+c] = swapRow[];
								inv.mData[currRow*c..currRow*c+c] = invSwapRow[];
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

				import std.algorithm : fold;
				auto allZeros = rref.mData[$-c..$].fold!((a, b) => (fabs(b) <= fpTol) && a)(true);

				if(allZeros)
				{
					return Nullable!ThisType();
				}
				else
				{
					return inv.nullable!ThisType;
				}
			}
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

	ThisType opUnary(string s)() immutable if (s == "-")
	{
		auto neg = ThisType(0);
		static if(isArray!T) {
			foreach(idx, ref arr; mData[]) {
				neg.mData[idx][] = -arr[];
			}
		} else {
			neg.mData[] = -mData[];
		}
		
		return neg;
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
			Vector!(3, T) cross(Vector!(3, T) rhs)
			{
				auto res = Vector!(3, T)(0);
				
				static if(isArray!T) {
					res.mData[0][] = mData[1][]*rhs.mData[2][] - mData[2][]*rhs.mData[1][];
					res.mData[1][] = mData[2][]*rhs.mData[0][] - mData[0][]*rhs.mData[2][];
					res.mData[2][] = mData[0][]*rhs.mData[1][] - mData[1][]*rhs.mData[0][];
				} else {
					res.mData[0] = mData[1]*rhs.mData[2] - mData[2]*rhs.mData[1];
					res.mData[1] = mData[2]*rhs.mData[0] - mData[0]*rhs.mData[2];
					res.mData[2] = mData[0]*rhs.mData[1] - mData[1]*rhs.mData[0];
				}
				return res;
			}
		}
		
		inout T dot(ref inout Vector!(r, T) rhs)
		{
			T res = 0;
			for(size_t i = 0; i < r; i++) {
				static if(isArray!T) {
					res[] += mData[i][]*rhs.mData[i][];
				} else {
					res += mData[i]*rhs.mData[i];
				}
			}
			
			return res;
		}

		inout T dot(Vector!(r, T) rhs)
		{
			T res = 0;
			for(size_t i = 0; i < r; i++) {
				static if(isArray!T) {
					res[] += mData[i][]*rhs.mData[i][];
				} else {
					res += mData[i]*rhs.mData[i];
				}
			}

			return res;
		}

		inout T dot(alias vector)()
		{
			T res = 0;
			for(size_t i = 0; i < r; i++) {
				static if(isArray!T) {
					res[] += mData[i][]*vector.mData[i][];
				} else {
					res += mData[i]*vector.mData[i];
				}
			}

			return res;
		}

		static if(isArray!T) {
			inout T dot(ref inout Vector!(r, ForeachType!T) rhs)
			{
				T res = 0;
				for(size_t i = 0; i < r; i++) {
					res[] += mData[i][]*rhs.mData[i];
				}
				return res;
			}
		}

		T magnitude() immutable
		{
			static if(isArray!T) {
				static assert(isNumeric!(ForeachType!T), "T is not a numeric type.");
				
				T res = 0;
				foreach(ref element; mData)
					res[] += element[]*element[];
			
				static if(isFloatingPoint!(ForeachType!T)) {
					foreach(ref mag_i ; res) {
						mag_i = sqrt(mag_i);
					}
				} else static if(isIntegral!(ForeachType!T)) {
					// cast to double if we are not a floating point type and round result.
					foreach(ref mag_i ; res) {
						mag_i = cast(ForeachType!T)round(sqrt(cast(double)mag_i));
					}
				} else {
					foreach(ref mag_i ; res) {
						mag_i = cast(ForeachType!T)sqrt(cast(double)mag_i);
					}
				}
				return res;

			} else {
				static assert(isNumeric!T, "T is not a numeric type.");

				T res = 0;
				foreach(ref element; mData)
					res += element*element;
				
				static if(isFloatingPoint!T) {
					res = sqrt(res);
				} else static if(isIntegral!T) {
					// cast to double if we are not a floating point type and round result.
					res = cast(T)round(sqrt(cast(double)res));
				} else {
					res = cast(T)sqrt(cast(double)res);
				}
				return res;

			}
		}

		T magnitude()
		{
			/+static assert(isNumeric!T, "T is not a numeric type.");

			T res = 0;
			foreach(ref element; mData)
				res += element*element;
			
			static if(isFloatingPoint!T) {
				res = sqrt(res);
			} else static if(isIntegral!T) {
				// cast to double if we are not a floating point type and round result.
				res = cast(T)round(sqrt(cast(double)res));
			} else {
				res = cast(T)sqrt(cast(double)res);
			}
			return res;+/
			static if(isArray!T) {
				static assert(isNumeric!(ForeachType!T), "T is not a numeric type.");
				
				T res = 0;
				foreach(ref element; mData)
					res[] += element[]*element[];
			
				static if(isFloatingPoint!(ForeachType!T)) {
					foreach(ref mag_i ; res) {
						mag_i = sqrt(mag_i);
					}
				} else static if(isIntegral!(ForeachType!T)) {
					// cast to double if we are not a floating point type and round result.
					foreach(ref mag_i ; res) {
						mag_i = cast(ForeachType!T)round(sqrt(cast(double)mag_i));
					}
				} else {
					foreach(ref mag_i ; res) {
						mag_i = cast(ForeachType!T)sqrt(cast(double)mag_i);
					}
				}
				return res;

			} else {
				static assert(isNumeric!T, "T is not a numeric type.");

				T res = 0;
				foreach(ref element; mData)
					res += element*element;
				
				static if(isFloatingPoint!T) {
					res = sqrt(res);
				} else static if(isIntegral!T) {
					// cast to double if we are not a floating point type and round result.
					res = cast(T)round(sqrt(cast(double)res));
				} else {
					res = cast(T)sqrt(cast(double)res);
				}
				return res;

			}

		}
		
		T magnitude_squared() immutable
		{
			T res = 0;
			static if(isArray!T) {
				static assert(isNumeric!(ForeachType!T), "T is not a numeric type.");

				foreach(ref element; mData)
					res[] += element[]*element[];
			} else {
				static assert(isNumeric!T, "T is not a numeric type.");

				foreach(ref element; mData)
					res += element*element;
			}
			return res;
		}

		T magnitude_squared()
		{
			T res = 0;
			static if(isArray!T) {
				static assert(isNumeric!(ForeachType!T), "T is not a numeric type.");

				foreach(ref element; mData)
					res[] += element[]*element[];
			} else {
				static assert(isNumeric!T, "T is not a numeric type.");

				foreach(ref element; mData)
					res += element*element;
			}
			return res;

		}

		ThisType normalize()
		{
			auto res = Vector!(r, T)(0);
			static if(isArray!T) {
				immutable T mag = 1/magnitude[];
				foreach(idx; 0..mag.length) {
					res.mData[idx][] = mData[idx][]*mag[];
				}
			} else {
				immutable mag = 1/magnitude;
				res.mData[] = mData[]*mag;
			}
			return res;
		}

		ThisType normalize() immutable
		{
			/+auto res = Vector!(r, T)(0);
			immutable mag = 1/magnitude;
			res.mData[] = mData[]*mag;
			return res;+/
			auto res = Vector!(r, T)(0);
			static if(isArray!T) {
				immutable T mag = 1/magnitude[];
				foreach(idx; 0..mag.length) {
					res.mData[idx][] = mData[idx][]*mag[];
				}
			} else {
				immutable mag = 1/magnitude;
				res.mData[] = mData[]*mag;
			}
			return res;
		}
	}

	enum size_t rows = r;
	enum size_t columns = c;

}
	string toString()
	{
		import std.format : format;
		import std.string : rightJustify;
		string matStr;
		
		for(int i = 0; i < r; i++)
		{
			matStr ~= "[";
			for(int j = 0; j < c; j++)
				matStr ~= " " ~ mData[i*c + j].format!"%s";
			
			matStr ~= " ]";
		}
		return matStr;
	}
}

// opEquals
@nogc @safe unittest
{
	auto m1 = Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	auto m2 = Matrix!(3, 2)(25, 28,
		57, 64,
		89, 100);
	assert(m1 == m2, "Matrix equality test failed");
}

// identity matrix generation test.
@nogc @safe unittest
{
	auto ident = Matrix!(3, 3).identity();
	auto testIdent = Matrix!(3, 3)(1, 0, 0, 0, 1, 0, 0, 0, 1);
	assert(ident == testIdent, "Matrix identity failed");
}

// op =
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 2)(1, 2,
		3, 4);
	auto m2 = Matrix!(2, 2)(0);
	m2 = m1;
	assert(m1 == m2, "Matrix assignment test failed");

	// Force failure
	m1.mData[0] = 3;
	assert(m1 != m2, "Matrix assignment test failed");
}

// opBinaryRight * (non square)
@nogc @safe unittest
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
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	auto m2 = Matrix!(2, 2)(5, 6, 7, 8);
	auto m3 = Matrix!(2, 2)(0);
	m3 = m1*m2;
	auto expected = Matrix!(2, 2)(19, 22, 43, 50);
	assert(m3 == expected, "Matrix opBinary * test failed");
}

// opBinary * (non-square)
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 3)(1, 2, 3, 4, 5, 6);
	auto m2 = Matrix!(3, 2)(5, 6, 7, 8,  9, 10);
	auto m3 = Matrix!(2, 2)(0);
	m3 = m1*m2;
	auto expected = Matrix!(2, 2)(46, 52, 109, 124);
	assert(m3 == expected, "Matrix opBinary * test failed");
}


// op +
@nogc @safe unittest
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
@nogc @safe unittest
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
@nogc @safe unittest
{
	auto m1 = Matrix!(3, 2)(2, 4,
							6, 8,
							10, 12);
	double scalar = 2;
	auto m2 = m1/scalar;
	auto expected = Matrix!(3, 2)(1, 2,
								  3, 4,
								  5, 6);
	assert(m2 == expected, "Matrix scalar division test failed");
}

// opUnary -
@nogc @safe unittest
{
	auto m1 = Matrix!(3, 3)(1, 2, 3, 4, 5, 6, 7, 8, 9);
	auto m2 = -m1;
	auto expected = Matrix!(3, 3)(-1, -2, -3, -4, -5, -6, -7, -8, -9);
	assert(m2 == expected, "Matrix negative test failed");
}

// op *=
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	m1 *= 2;
	auto expected = Matrix!(2, 2)(2, 4, 6, 8);
	assert(m1 == expected, "Matrix *= test failed");
}

// op +=
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 2)(1, 2, 3, 4);
	m1 += Matrix!(2, 2).identity();
	auto expected = Matrix!(2, 2)(2, 2, 3, 5);
	assert(m1 == expected, "Matrix += test failed");
}

// transpose
@nogc @safe unittest
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
@nogc @safe unittest
{
	auto m1 = Matrix!(2, 2)(2, 4, 6, 8);
	auto m2 = m1.rref();
	auto expected1 = Matrix!(2, 2).identity();
	assert(m2 == expected1, "Matrix rref test failed");

	auto m3 = Matrix!(3, 3)(2, 4, 6, 8, 10, 12, 14, 16, 18);
	auto m4 = m3.rref();
	auto expected2 = Matrix!(3, 3)(1, 0, -1, 0, 1, 2, 0, 0, 0);
	assert(m4 == expected2, "Matrix rref test failed");

	auto m5 = Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m6 = m5.rref();
	auto expected3 = Matrix!(3, 3).identity();
	assert(m6 == expected3, "Matrix rref test failed");

	/*
	 0     4     5
     0     0     7
     0     3     1
	*/
	auto m7 = Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m8 = m7.rref();
	auto expected4 = Matrix!(3, 3)(0, 1, 0, 0, 0, 1, 0, 0, 0);
	assert(m8 == expected4, "Matrix rref test failed");


	auto m9 = Matrix!(3, 6)(0, 0, 2, 1, 0, 0, 0, 5, 1, 0, 1, 0, 6, 3, 1, 0, 0, 1);
	auto m10 = m9.rref();
	/*
	-0.0333   -0.1000    0.1667
	-0.1000    0.2000        0
	0.5000         0         0
	*/
	auto expected5 = Matrix!(3, 6)(1, 0, 0, -1.0/30.0, -1.0/10.0, 1.0/6.0, 0, 1, 0, -1.0/10.0, 1.0/5.0, 0, 0, 0, 1.0, 1.0/2.0, 0, 0);
	assert(m10 == expected5, "Matrix rref test failed");
}

@nogc @safe unittest
{
	auto m1 = Matrix!(3, 3)(0, 0, 2, 0, 5, 1, 6, 3, 1);
	auto m2 = m1.inverse();
	auto expected = Matrix!(3, 3)(-1.0/30.0, -1.0/10.0, 1.0/6.0, -1.0/10.0, 1.0/5.0, 0, 1.0/2.0, 0, 0);
	//assert(m2 == expected, "Matrix inverse test failed");

	auto m3 = Matrix!(3, 3)(0, 4, 5, 0, 0, 7, 0, 3, 1);
	auto m4 = m3.inverse();
	assert(m4.isNull, "Matrix inverse test failed, should have been empty Nullable");
}

@nogc @safe unittest
{
	auto vec = Vector!(3)(1.0, 1.0, 0.0);
	auto mag = vec.magnitude();
	auto expected = sqrt(cast(double)2);	
	assert(mag == expected, "Vector magnitude test failed");
}

@nogc @safe unittest
{
	auto vec1 = Vector!(3)(1, 1, 0);
	auto vec2 = vec1.normalize();
	auto expected = Vector!(3)(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0);
	assert(vec2 == expected, "Vector normalize test failed");
}

@nogc @safe unittest
{
	auto vec1 = Vector!(3)(1, 2, 3);
	auto vec2 = Vector!(3)(4, 5, 6);
	auto res = vec1.dot(vec2);
	double expected = 32;
	assert(res == expected, "Vector dot product test failed");
}

@nogc @safe unittest
{
	auto vec1 = Vector!(3)(1, 2, 3);
	auto vec2 = Vector!(3)(4, 5, 6);
	auto vec3 = vec1.cross(vec2);
	auto expected = Vector!(3)(-3, 6, -3);
	assert(vec3 == expected, "Vector cross product failed");
}
