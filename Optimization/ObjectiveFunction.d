module MDOL.ObjectiveFunction;

import std.complex;
import std.math;
import std.stdio;

import MDOL.Gradient;
import MDOL.Derivative;
import MDOL.ComplexStep;
import MDOL.FiniteDifference;

abstract class ObjectiveFunction : IGradient
{
	Complex!double Compute(Complex!double[] designVar);

	Complex!double[] Constraint(Complex!double[] designVar);

	double[] Gradient(double[] point)
	{
		//writeln("Computing gradient");
		auto diff = cast(Derivative)Object.factory(DerivativeType);
		double[] grad = new double[point.length];

		if(diff is null)
		{
			writeln("Cannot allocate derivative machine");
			return grad;
		}
		diff.StepSize = StepSize;
		//writeln("alloc derivative");

		for(int i = 0; i < grad.length; i++)
		{
			//writeln(i);
			grad[i] = diff.Compute(this, point, i+1);
		}
		
		return grad;
	}

	double[][] ConstraintGradient(double[] point)
	{
		//writeln("Computing gradient");
		auto diff = cast(Derivative)Object.factory(DerivativeType);
		int i = 0;
		int j = 0;
		double[][] grad;// = new double[point.length];
		grad.length = Constraints;
		for(i = 0; i < Constraints; i++)
		{
			grad[i].length = point.length;
		}
		
		if(diff is null)
		{
			writeln("Cannot allocate derivative machine");
			return grad;
		}
		diff.StepSize = StepSize;

		//alias Complex!double delegate(Complex!double[] input) DerivativeEquation;
		Complex!double ConstraintDerivative(Complex!double[] input)
		{
			return this.Constraint(input)[i];
		}
		//writeln("alloc derivative");
		
		for(i = 0; i < grad.length; i++)
		{
			for(j = 0; j < grad[i].length; j++)
			{
			//writeln(i);
				grad[i][j] = diff.Compute(&ConstraintDerivative, point, j+1);
			}
		}
		
		return grad;
	}

	@property void Constraints(int constraints) { mConstraints = constraints; }
	@property int Constraints() { return mConstraints; }

	@property void StepSize(double step) { h = step; }
	@property double StepSize() { return h; }

	@property DerivativeType(string derType) { mDerivativeType = derType; }
	@property DerivativeType() { return mDerivativeType; }

	private string mDerivativeType = "MDOL.ComplexStep.ComplexStep";
	private double h = 1.0e-10;
	private int mConstraints = 0;
}