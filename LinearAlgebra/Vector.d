module LinearAlgebra.Vector;

import LinearAlgebra.Matrix;

import std.stdio;		
alias Vector(size_t l, T = real) = Matrix!(l, 1, T);

Vector!(3) cross(Vector!(3) lhs, Vector!(3) rhs)
{
	auto res = new Vector!(3);

	res.mData[0] = lhs.mData[1]*rhs.mData[2] - lhs.mData[2]*rhs.mData[1];
	res.mData[1] = lhs.mData[2]*rhs.mData[0] - lhs.mData[0]*rhs.mData[2];
	res.mData[2] = lhs.mData[0]*rhs.mData[1] - lhs.mData[1]*rhs.mData[0];

	return res;
}

unittest
{
	auto vec1 = new Vector!(3)(1, 2, 3);
	auto vec2 = new Vector!(3)(4, 5, 6);

	auto vec3 = vec1.cross(vec2);

	auto expected = new Vector!(3)(-3, 6, -3);

	assert(vec3 == expected);
}