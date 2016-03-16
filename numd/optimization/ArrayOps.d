module numd.optimization.ArrayOps;

import std.complex;
import std.math;
import std.random;
import std.stdio;

double[] MultiplyScalar(double[] arr, double scal)
{
	for(int i = 0; i < arr.length; i++)
	{
		arr[i] *= scal;
	}
	return arr;
}

Complex!double[] complex(double[] rhs)
{
	Complex!double[] lhs = new Complex!double[rhs.length];
	//lhs.length = rhs.length;
	foreach(int i, double val; rhs)
	{
		lhs[i] = std.complex.complex!(double)(val, 0.0);
	}
	return lhs;
}

double[] Real(Complex!double[] rhs)
{
	double[] lhs = new double[rhs.length];
	//lhs.length = rhs.length;
	foreach(int i, Complex!double val; rhs)
	{
		lhs[i] = val.re;
	}
	return lhs;
}

double[] abs(double[] rhs)
{
	double[] lhs = new double[rhs.length];
	//lhs.length = rhs.length;
	foreach(int i, double val; rhs)
	{
		lhs[i] = std.math.abs(val);
	}
	return lhs;
}

void WriteArrayCSV(ref File f, double[] x)
{
	for(int i = 0; i < x.length-1; i++)
		f.writef("%40.40f, ", x[i]);

	f.writef("%40.40f\n", x[$-1]);
}

double[] RandomVector(int size, double min, double max)
{
	double[] rnd = new double[size];
	auto gen = Random(unpredictableSeed);

	//Mt19937 gen;
	//gen.min!(double)(min);
	//gen.max!(double)(max);
	//gen.seed(unpredictableSeed);
	
	for(int i = 0; i < size; i++)
	{
		rnd[i] = uniform!("[]", double, double)(min, max, gen);
		//writeln(val);
		//val = gen.front;
		//gen.popFront();
	}
	
	return rnd;
}