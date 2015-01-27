module Optimization.SecantRoot;

import std.complex;
import std.math;
import std.stdio;

import Optimization.RootFinder;

class SecantRoot : RootFinder
{
	final override Complex!double doSolve(RootEquation eqn, Complex!double xk)
	{
		int iterations = 0;
		Complex!double xkNext;
		Complex!double xkPrev = xk - 1;
		auto f = complex!double(1.0, 0);
		
		while(abs(f.re) > Tolerance)
		{
			f = eqn(xk);
			auto fp = eqn(xkPrev);
			xkNext = xk - f*((xk - xkPrev)/(f - fp));
			xkPrev = xk;
			xk = xkNext;
			//writeln(f.re);
			iterations++;
		}
		//writeln("Secant root finder converged in ", iterations, " iterations.");
		return xk;
	}
}