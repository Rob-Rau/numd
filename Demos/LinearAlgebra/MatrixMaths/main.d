module main;

import std.stdio;

import numd.linearalgebra.matrix;

void main(string[] args)
{
	auto m1 = Matrix!(3, 2)(1, 2,
							3, 4,
							5, 6);

	auto m2 = Matrix!(2, 2)(7, 8,
							9, 10);

	auto m3 = Matrix!(3, 2)(7, 8,
							9, 10,
							11, 12);

	auto m4 = m1 * m2;

	auto m5 = m1 + m3;

	writeln("m1:");
	writeln(m1);

	writeln("m2:");
	writeln(m2);

	writeln("m3:");
	writeln(m3);

	writeln("m4 = m1 * m2:");
	writeln(m4);

	writeln("m5 = m1 + m3:");
	writeln(m5);
}

