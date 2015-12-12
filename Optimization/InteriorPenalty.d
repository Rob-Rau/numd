module Optimization.InteriorPenalty;

import Optimization.ArrayOps;
import Optimization.BfgsNewton;
import Optimization.BracketAndZoom;
import Optimization.Complex;
import Optimization.MatrixOps;
import Optimization.ObjectiveFunction;
import Optimization.Optimizer;

import core.thread;

import std.algorithm;
import std.complex;
import std.conv;
import std.math;
import std.stdio;

import cblas;
//import scid.bindings.lapack.dlapack;
//import scid.bindings.blas.dblas;

abstract class InteriorPenaltyFunction : ObjectiveFunction
{
	@property void ObjectiveFunc(ObjectiveFunction objectiveFunc) { mObectiveFunction = objectiveFunc; }
	@property ObjectiveFunction ObjectiveFunc() { return mObectiveFunction; }
	
	@property void Mu(double mu) { mMu = mu; }
	@property double Mu() { return mMu; }
	
	private ObjectiveFunction mObectiveFunction;
	private double mMu;
}

class LogarithmicPenalty : InteriorPenaltyFunction
{
	final override Complex!double Compute(Complex!double[] designVar)
	{
		Complex!double logSum = 0.0;
		for(int i = 0; i < ObjectiveFunc.Constraints; i++)
		{
			auto tmp = log(ObjectiveFunc.Constraint(designVar)[i]);
			//writeln("Here");

			if(isNaN(tmp.re) || isInfinity(tmp.re))
			{
				//tmp = complex(-std.math.log(abs(ObjectiveFunc.Constraint(designVar)[i].re))-10000, tmp.im);
				tmp = complex(log(0), tmp.im);
				//tmp = complex(10000, tmp.im);
			}

			logSum += tmp;
		}
		return ObjectiveFunc.Compute(designVar) - Mu*logSum;
	}
	
	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}
}

class InteriorPenalty : Optimizer
{

	final override protected Result doOptimize(ObjectiveFunction objectiveFunction)
	{
		version(X86)
		{
			uint iterations = 0;
			uint minorIterations = 0;
		}
		else
		{
			ulong iterations = 0;
			ulong minorIterations = 0;
		}
		double decrease = 0.1;
		Result result;
		Result innerRes;
		Result innerResLast;
		InteriorPenaltyFunction penalty = new LogarithmicPenalty;
		Optimizer optimizer = new BfgsNewton!9;
		optimizer.DebugMode = DebugMode;
		//optimizer.PointFilename = PointFilename;
		//optimizer.PointFilename = BfgsPoints.;
		//optimizer.ErrorFilename = ErrorFilename;
		optimizer.FileOutput = FileOutput;
		//auto f = File(PointFilename, "w");
		//auto ferr = File(ErrorFilename, "w");
		File f;
		if(FileOutput) f = File(PointFilename, "w");
		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");

		auto activeConstraint = minPos( Real(objectiveFunction.Constraint(complex(InitialGuess))))[0];
		double[] cGrad = objectiveFunction.ConstraintGradient(InitialGuess)[0];
		double[] cGradNorm = new double[cGrad.length];
		cGradNorm[] = nrm2(cast(int)cGrad.length, cast(double*)cGrad, 1)^^(-1)*cGrad[];
		double[] testPoint = new double[cGrad.length];
		testPoint[] = InitialGuess;
		if(DebugMode) writeln("Grad C = ", cGrad);

		double alpha = 1;

		if(FileOutput) WriteArrayCSV(f, InitialGuess);

		while(activeConstraint < 0)
		{
			testPoint[] = alpha*cGradNorm[] + InitialGuess[];
			activeConstraint = minPos( objectiveFunction.Constraint(testPoint.complex()).Real())[0];
			if(DebugMode) writeln("ActiveConstraint = ", activeConstraint, "\tGrad C = ", cGrad);
			alpha = 1.5*alpha;
		}

		InitialGuess[] = testPoint[];

		if(FileOutput) WriteArrayCSV(f, InitialGuess);

		optimizer.InitialGuess = new double[InitialGuess.length];
		optimizer.InitialGuess[] = InitialGuess;

		penalty.ObjectiveFunc = objectiveFunction;
		penalty.Mu = 10;

		optimizer.PointFilename = "BfgsPoints" ~ to!string(iterations) ~ ".csv";
		optimizer.ErrorFilename = "BfgsError" ~ to!string(iterations) ~ ".csv";

		if(DebugMode) writeln("Error: ", 100, "\t\tActive constraint: ", activeConstraint, "\t\tPoint: ", InitialGuess);
		innerResLast = optimizer.Optimize(penalty);
		if(FileOutput) WriteArrayCSV(f, innerResLast.DesignVariables);
		auto muLast = penalty.Mu;
		//innerResLast = innerRes;
		penalty.Mu = decrease*muLast;

		optimizer.InitialGuess[] = innerResLast.DesignVariables;
		//if(DebugMode) writeln("Made it here");
		innerRes = optimizer.Optimize(penalty);

		double error = (innerRes.ObjectiveFunctionValue - innerResLast.ObjectiveFunctionValue)/(penalty.Mu - muLast);
		activeConstraint = minPos( Real(objectiveFunction.Constraint(complex(innerRes.DesignVariables))))[0];

		if(DebugMode) writeln("Error: ", (error), "\t\tActive constraint: ", activeConstraint, "\t\tPoint: ", innerRes.DesignVariables);

		while( (abs(activeConstraint) > Tolerance) || (activeConstraint < 0.0))
		{
			optimizer.PointFilename = "BfgsPoints" ~ to!string(iterations+1) ~ ".csv";
			optimizer.ErrorFilename = "BfgsError" ~ to!string(iterations+1) ~ ".csv";

			innerResLast = innerRes;
			muLast = penalty.Mu;
			penalty.Mu = decrease*muLast;
			optimizer.InitialGuess[] = innerResLast.DesignVariables;
			innerRes = optimizer.Optimize(penalty);
			if(FileOutput) WriteArrayCSV(f, innerRes.DesignVariables);
			optimizer.InitialGuess[] = innerRes.DesignVariables;

			error = (innerRes.ObjectiveFunctionValue - innerResLast.ObjectiveFunctionValue)/(penalty.Mu - muLast);
			activeConstraint = minPos( Real(objectiveFunction.Constraint(complex(innerRes.DesignVariables))))[0];

			if(DebugMode) writeln("Error: ", error, "\t\tActive constraint: ", activeConstraint, "\t\tPoint: ", innerRes.DesignVariables);

			minorIterations += innerRes.Iterations;
			//if(DebugMode) writeln(iterations);
			iterations++;
		}

		if(DebugMode) writeln();
		double[] grad = penalty.Gradient(innerRes.DesignVariables);
		if(DebugMode) writeln("Gradient", grad, "\tNorm = ", nrm2(cast(int)grad.length, cast(double*)(grad), 1));

		result.ObjectiveFunctionValue = innerRes.ObjectiveFunctionValue;
		result.Iterations = iterations;
		result.MinorIterations = minorIterations;
		result.DesignVariables = innerRes.DesignVariables;

		return result;
	}
}