module numd.optimization.ComplexStep;

import std.complex;
import std.stdio;

import numd.optimization.Derivative;
import numd.optimization.ObjectiveFunction;

class ComplexStep : Derivative
{
	final override double Compute(ObjectiveFunction func, in double[] point, int dimension)
	{
		//writeln("In ComplexStep");
		Complex!double[] cPoint;
		cPoint.length = point.length;

		// Complexify the point
		for(int i = 0; i < cPoint.length; i++)
		{
			//writeln("Complixify i: ", i);
			if(i == (dimension-1))
			{
				//writeln("Found dimension: ", i+1);
				cPoint[i].re = point[i];
				cPoint[i].im = StepSize;
			}
			else
			{
				cPoint[i].re = point[i];
				cPoint[i].im = 0;
			}
		}
		//writeln("Computing derivative");
		Complex!double f = func.Compute(cPoint);
		//writeln("Done Computing derivative");
		return f.im/StepSize;
	}

	final override double Compute(DerivativeEquation func, in double[] point, int dimension)
	{
		//writeln("In ComplexStep");
		Complex!double[] cPoint;
		cPoint.length = point.length;
		
		// Complexify the point
		for(int i = 0; i < cPoint.length; i++)
		{
			//writeln("Complixify i: ", i);
			if(i == (dimension-1))
			{
				//writeln("Found dimension: ", i+1);
				cPoint[i].re = point[i];
				cPoint[i].im = StepSize;
			}
			else
			{
				cPoint[i].re = point[i];
				cPoint[i].im = 0;
			}
		}
		//writeln("Computing derivative");
		Complex!double f = func(cPoint);
		//writeln("Done Computing derivative");
		return f.im/StepSize;
	}
}