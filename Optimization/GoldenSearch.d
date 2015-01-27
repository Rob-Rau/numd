module Optimization.GoldenSearch;

import Optimization.Optimizer;
import Optimization.ArrayOps;
import Optimization.ObjectiveFunction;

import core.thread;
import std.array;
import std.algorithm;
import std.complex;
import std.stdio;
import std.math;

const double Tau = 0.618033988749895;

class GoldenSearch : Optimizer
{
	final override protected Result doOptimize(ObjectiveFunction objectiveFunction)
	{
		Result result;
		ulong numDesignVars = InitialGuess.length;
		version(X86)
		{
			uint iterations = 0;
		}
		else
		{
			ulong iterations = 0;
		}
		double minPoint = 1;
		double minPointOld = 1;
		double minPointErr = 1;
		double minVar;
		double minVarOld = InitialGuess[0];
		double intStart = InitialGuess[0];
		double tau = InitialGuess[0] + (Tau*Interval);
		double oneMinusTau = InitialGuess[0] + (1 - Tau)*Interval;
		double intEnd = InitialGuess[0] + Interval;

		Complex!double[] points;
		points.length = 2;

		points[0] = objectiveFunction.Compute(array([complex(oneMinusTau)]));
		points[1] = objectiveFunction.Compute(array([complex(tau)]));

		result.Error.length += 100;
		result.IterArr.length += 100;
		while(minPointErr > Tolerance)
		{
			// Find the smaller interior point
			if(points[0].re < points[1].re)
			{
				double newTau = oneMinusTau;
				double newOneMinusTau = (tau - intStart)*(1-Tau) + intStart;
				double newIntEnd = tau;
				minVar = oneMinusTau;
				intEnd = newIntEnd;
				tau = newTau;
				oneMinusTau = newOneMinusTau;
				points[1] = points[0];
				minPoint = points[0].re;
				points[0] = objectiveFunction.Compute(array([complex(oneMinusTau)]));
			}
			else
			{
				double newTau = Tau*(intEnd - oneMinusTau) + oneMinusTau;
				double newOneMinusTau = tau;
				double newIntStart = oneMinusTau;
				minVar = oneMinusTau;
				intStart = newIntStart;
				tau = newTau;
				oneMinusTau = newOneMinusTau;
				points[0] = points[1];
				minPoint = points[1].re;
				points[1] = objectiveFunction.Compute(array([complex(tau)]));
			}

			minPointErr = abs(minPoint - minPointOld)/(1+abs(minPoint)) + abs(minVar - minVarOld)/(1+abs(minVar));
			if(iterations >= result.Error.length)
			{
				result.Error.length += 100;
				result.IterArr.length += 100;
			}

			result.Error[iterations] = minPointErr;
			result.IterArr[iterations] = iterations;

			//writeln(minPointErr);
			minVarOld = minVar;
			minPointOld = minPoint;

			iterations++;
		}

		result.Iterations = iterations;
		result.ObjectiveFunctionValue = minPoint;
		result.DesignVariables = array([tau]);
		return result;
	}
}