module main;

import std.stdio;

import LinearAlgebra.Matrix;
import LinearAlgebra.Vector;

void main(string[] args)
{
	auto m1 = new Matrix!(3, 2)(1, 2,
								3, 4,
								5, 6);
	auto m2 = new Matrix!(2, 2)(7, 8,
								9, 10);

	auto m3 = new Matrix!(3, 2)(7, 8,
								9, 10,
								11, 12);

	auto m4 = m1 * m2;

	auto m5 = m1 + m3;

	writeln("m1:");
	writeln(m1.ToString());

	writeln("m2:");
	writeln(m2.ToString());

	writeln("m3:");
	writeln(m3.ToString());

	writeln("m4 = m1 * m2:");
	writeln(m4.ToString());

	writeln("m5 = m1 + m3:");
	writeln(m5.ToString());

	if(m2 == m3)
		writeln("m2 == m3");

	if(m1 != m2)
		writeln("m1 != m2");

	if(m1 != m4)
		writeln("m1 != m4");
}

