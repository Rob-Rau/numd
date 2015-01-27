module MDOL.Derivative;

import std.complex;

import MDOL.ObjectiveFunction;

alias Complex!double delegate(Complex!double[] input) DerivativeEquation;

abstract class Derivative
{
	double Compute(ObjectiveFunction func, double[] point, int dimension);
	double Compute(DerivativeEquation func, double[] point, int dimension);

	@property void StepSize(double step) { h = step; }
	@property double StepSize() { return h; }

	private double h = 1.0e-10;
}