module Optimization.Derivative;

import std.complex;

import Optimization.ObjectiveFunction;

alias Complex!double delegate(Complex!double[] input) DerivativeEquation;

abstract class Derivative
{
	double Compute(ObjectiveFunction func, in double[] point, int dimension);
	double Compute(DerivativeEquation func, in double[] point, int dimension);

	@property void StepSize(double step) { h = step; }
	@property double StepSize() { return h; }

	private double h = 1.0e-10;
}