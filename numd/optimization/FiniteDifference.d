module numd.optimization.FiniteDifference;

import std.complex;
import std.stdio;

import numd.optimization.Derivative;
import numd.optimization.ObjectiveFunction;

class FiniteDifference : Derivative
{
	final override double Compute(ObjectiveFunction func, in double[] point, int dimension)
	{
		//writeln("In finite diff");
		Complex!double[] leftPoint, rightPoint;

		leftPoint.length = point.length;
		rightPoint.length = point.length;

		//writeln("going to complexify");
		//leftPoint = cast(Complex!double[])point;
		//rightPoint = cast(Complex!double[])point;

		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			//writeln("Complixify i: ", i);
			if(i == (dimension-1))
			{
				//writeln("Found dimension");
				leftPoint[i].re = point[i] + StepSize;
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i] - StepSize;
				rightPoint[i].im = 0;
			}
			else
			{
				leftPoint[i].re = point[i];
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i];
				rightPoint[i].im = 0;
			}
		}

		//writeln(leftPoint);
		//writeln("setting step");
		//leftPoint[dimension-1] += StepSize;
		//rightPoint[dimension-1] -= StepSize;

		//writeln("trying to compute");
		return (func.Compute(leftPoint).re - func.Compute(rightPoint).re)/(2*StepSize);
	}

	final override double Compute(DerivativeEquation func, in double[] point, int dimension)
	{
		//writeln("In finite diff");
		Complex!double[] leftPoint, rightPoint;
		
		leftPoint.length = point.length;
		rightPoint.length = point.length;
		
		//writeln("going to complexify");
		//leftPoint = cast(Complex!double[])point;
		//rightPoint = cast(Complex!double[])point;
		
		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			//writeln("Complixify i: ", i);
			if(i == (dimension-1))
			{
				//writeln("Found dimension");
				leftPoint[i].re = point[i] + StepSize;
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i] - StepSize;
				rightPoint[i].im = 0;
			}
			else
			{
				leftPoint[i].re = point[i];
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i];
				rightPoint[i].im = 0;
			}
		}
		
		//writeln(leftPoint);
		//writeln("setting step");
		//leftPoint[dimension-1] += StepSize;
		//rightPoint[dimension-1] -= StepSize;
		
		//writeln("trying to compute");
		return (func(leftPoint).re - func(rightPoint).re)/(2*StepSize);
	}
}