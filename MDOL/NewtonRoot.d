module MDOL.NewtonRoot;

import std.complex;
import std.math;
import std.stdio;

import MDOL.RootFinder;

class NewtonRoot : RootFinder
{
	final override Complex!double doSolve(RootEquation eqn, Complex!double xk)
	{
		int iterations = 0;
		double xkNext;
		auto f = complex!double(1.0, 0);
		xk.im = h;
		while(abs(cast(double)f.re) > Tolerance)
		{
			f = eqn(xk);
			double fp = f.im/h;
			xkNext = xk.re - f.re/fp;
			xk.re = xkNext;
			//xk.im = h;
			//writeln(f.re);
			iterations++;
		}
		writeln("Newton root finder converged in ", iterations, " iterations.");
		return xk;
	}
}