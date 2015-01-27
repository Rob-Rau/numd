module MDOL.RootFinder;

import std.complex;

/*
iNterface RootEquation
{
	Complex!double Equation(Complex!double input);
}
*/

alias Complex!double delegate(Complex!double input) RootEquation;

abstract class RootFinder
{
	public Complex!double Solve(RootEquation eqn, Complex!double xk)
	{
		return doSolve(eqn, xk);
	}

	protected Complex!double doSolve(RootEquation, Complex!double xk);

	@property void Tolerance(double tolerance) { mTolerance = tolerance; }
	@property double Tolerance() { return mTolerance; }

	double mTolerance = 1.0e-10;
	double h = 10.0e-40;
}