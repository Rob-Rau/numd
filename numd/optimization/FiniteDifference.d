module numd.optimization.FiniteDifference;

import std.complex;
import std.stdio;

import numd.optimization.Derivative;
import numd.optimization.ObjectiveFunction;

class FiniteDifference : Derivative
{
	final override double Compute(ObjectiveFunction func, in double[] point, int dimension)
	{
		Complex!double[] leftPoint, rightPoint;

		leftPoint.length = point.length;
		rightPoint.length = point.length;
		
		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			if(i == (dimension-1))
			{
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

		return (func.Compute(leftPoint).re - func.Compute(rightPoint).re)/(2*StepSize);
	}

	final override double Compute(DerivativeEquation func, in double[] point, int dimension)
	{
		Complex!double[] leftPoint, rightPoint;
		
		leftPoint.length = point.length;
		rightPoint.length = point.length;

		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			if(i == (dimension-1))
			{
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

		return (func(leftPoint).re - func(rightPoint).re)/(2*StepSize);
	}
}

class FiniteDifferenceEqualized : Derivative
{
	final override double Compute(ObjectiveFunction func, in double[] point, int dimension)
	{
		//writeln("In finite diff");
		Complex!double[] leftPoint, rightPoint;

		leftPoint.length = point.length;
		rightPoint.length = point.length;

		double addBack = StepSize/(cast(double)point.length - 1.0);


		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			if(i == (dimension-1))
			{
				leftPoint[i].re = point[i] + StepSize;
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i] - StepSize;
				rightPoint[i].im = 0;
			}
			else
			{
				leftPoint[i].re = point[i] - addBack;
				leftPoint[i].im = 0;
				rightPoint[i].re = point[i] + addBack;
				rightPoint[i].im = 0;
			}
		}

		return (func.Compute(leftPoint).re - func.Compute(rightPoint).re)/(2*StepSize);
	}

	final override double Compute(DerivativeEquation func, in double[] point, int dimension)
	{
		Complex!double[] leftPoint, rightPoint;
		
		leftPoint.length = point.length;
		rightPoint.length = point.length;
		double addBack = StepSize/(cast(double)point.length - 1.0);

		
		// Complexify the point
		for(int i = 0; i < point.length; i++)
		{
			if(i == (dimension-1))
			{
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

		return (func(leftPoint).re - func(rightPoint).re)/(2*StepSize);
	}
}