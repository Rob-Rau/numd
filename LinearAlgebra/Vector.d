module LinearAlgebra.Vector;

import LinearAlgebra.Matrix;

import std.stdio;
import std.math;

class Vector(size_t l, T = real) : Matrix!(l, 1, T)
{

	this(in T[l] values...)
	{
		super(values);
	}
	
	this()
	{
		super();
	}
	
	this(Vector!(l, T) vec)
	{
		super(vec);
	}

	Vector!(3) cross(Vector!(3) rhs)
	{
		auto res = new Vector!(3);

		res.mData[0] = mData[1]*rhs.mData[2] - mData[2]*rhs.mData[1];
		res.mData[1] = mData[2]*rhs.mData[0] - mData[0]*rhs.mData[2];
		res.mData[2] = mData[0]*rhs.mData[1] - mData[1]*rhs.mData[0];

		return res;
	}

	T dot(Vector rhs)
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

	Vector normalize()
	{
		auto res = new Vector!(l, T);
		res.mData[] = mData[]/magnitude();
		return res;
	}
}

unittest
{
	auto vec = new Vector!(3)(1.0, 1.0, 0.0);
	auto mag = vec.magnitude();
	auto expected = sqrt(cast(real)2);
	
	assert(mag == expected, "Vector magnitude test failed");
}

unittest
{
	auto vec1 = new Vector!(3)(1, 1, 0);
	auto vec2 = vec1.normalize();
	auto expected = new Vector!(3)(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0);
	
	assert(vec2 == expected, "Vector normalize test failed");
}

unittest
{
	auto vec1 = new Vector!(3)(1, 2, 3);
	auto vec2 = new Vector!(3)(4, 5, 6);
	
	auto res = vec1.dot(vec2);
	
	real expected = 32;
	
	assert(res == expected, "Vector dot product test failed");
}

unittest
{
	auto vec1 = new Vector!(3)(1, 2, 3);
	auto vec2 = new Vector!(3)(4, 5, 6);
	
	auto vec3 = vec1.cross(vec2);
	
	auto expected = new Vector!(3)(-3, 6, -3);
	
	assert(vec3 == expected, "Vector cross product failed");
}