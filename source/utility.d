module numd.utility;

import std.conv;

/++
+ This is the tolerance used when doing floating point comparisons.
+/
public double fpTol = 1.0e-16;

struct Meshgrid(T = double)
{
	T[][] X;
	T[][] Y;
}

Meshgrid!T meshgrid(T = double)(T[] x, T[] y)
{
	Meshgrid!T mesh = {X: new T[][](y.length, x.length), Y: new T[][](y.length, x.length)};

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.X[i][j] = x[j];

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.Y[i][j] = y[i];

	return mesh;
}

T[] linspace(T = double)(T start, T end, int points)
{
	//real h = (end-start)/(points-1);
	T h = ((end-start)/(points)).to!T;
	T[] x = new T[points];
	x[0] = start;
	for(int i = 1; i < points; i++)
	{
		x[i] = x[i-1]+h;
	}
	return x;
}