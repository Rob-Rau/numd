module Optimization.MatrixOps;

import std.stdio;

double[] IdentityMatrix(int size)
{
	double[] I = new double[size^^2];
	I[] = 0;
	for(int i = 0; i < size^^2; i += size+1)
	{
		I[i] = 1;
	}
	return I;
}

void PrintMatrix(double[] mat, int clen)
{
	int rlen = cast(int)mat.length/clen;
	for(int i = 0; i < clen; i++)
	{
		write("[ ");
		for(int j = 0; j < mat.length; j+=clen)
		{
			write(mat[i+j], " ");
		}
		write("]\n");
	}
}