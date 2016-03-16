module numd.optimization.ConjugateGradient;

import numd.optimization.ArrayOps;
import numd.optimization.Optimizer;
import numd.optimization.ObjectiveFunction;
import numd.optimization.BracketAndZoom;

import core.thread;

import std.math;
import std.stdio;

//import scid.bindings.lapack.dlapack;
//import scid.bindings.blas.dblas;
import cblas;

class ConjugateGradient : Optimizer
{
	enum UpdateMethods
	{
		FletcherReeves,
		PolakRibiere,
		DaiYuan
	}

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
		auto lineSearch = new BracketAndZoom;
		lineSearch.DebugMode = DebugMode;
		Result lineResult;
		Result lineResultLast = lineResult;
		Result result;
		bool converged = false;
		double error = 1;
		double[] xk = InitialGuess;
		double[] gk = objectiveFunction.Gradient(InitialGuess);
		double[] pkLast = new double[gk.length];
		double[] gkLast = new double[gk.length];
		double[] xkLast = new double[gk.length];
		double[] pk = new double[gk.length];//gk;//new double[gk.length];
		double[] tmp = new double[gk.length];
		auto f = File(PointFilename, "w");
		auto ferr = File(ErrorFilename, "w");
		//auto f = File("ConjugatePoints.csv", "w");
		//auto ferr = File("ConjugateError.csv", "w");
		double beta;
		WriteArrayCSV(f, xk);
		//lineSearch.P.length = pk.length;
		pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);
		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		
		lineSearch.P = pk;
		lineSearch.InitialGuess = xk;
		lineResult = lineSearch.Optimize(objectiveFunction);
		minorIterations += lineResult.Iterations;
		lineResultLast = lineResult;
		xkLast[] = xk;
		pkLast[] = pk;
		gkLast[] = gk;
		xk = lineResult.DesignVariables;
		if(DebugMode) writeln("xk = ", xk);
		
		// Compute gradient at new point.
		gk[] = objectiveFunction.Gradient(xk);
		// Compute new direction at point.
		pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);
		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		lineSearch.AlphaInitial = (dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));

		double epa = 1.0e-6;
		double epr = 1.0e-6;
		double epg = 1.0e-3;

		//while(error >= 1.0e-6)
		while( !(abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) > (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && !(nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1) <= epg) )
		{
			if(DebugMode) writeln();
			if(FileOutput) WriteArrayCSV(f, xk);
			//writeln("Alpha = ", lineSearch.AlphaInitial);
			//lineSearch.AlphaInitial = 1;
			//lineSearch.AlphaInitial = lineSearch.AlphaInitial*(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//writeln("gk-1*pk-1 = ", dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1));
			//writeln("gk*pk = ", dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//writeln("AlphaNew = ", lineSearch.AlphaInitial);
			lineSearch.P = pk;
			
			lineSearch.InitialGuess[] = xk;
			
			lineResult = lineSearch.Optimize(objectiveFunction);
			minorIterations += lineResult.Iterations;
			//lineResultLast = lineResult;
			xkLast[] = xk;
			
			pkLast[] = pk;
			
			gkLast[] = gk;
			
			xk[] = lineResult.DesignVariables;
			
			// Compute gradient at new point.
			gk = objectiveFunction.Gradient(xk);
			// Compute new direction at point.
			//pk[] = gk;
			//if(iterations%1.0)
			if(iterations%5)
			//if(iterations%xk.length)
			{
				if(UpdateMethod == UpdateMethods.FletcherReeves)
				{
					beta = dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)gk, 1)/dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)gkLast, 1);
				}
				else if(UpdateMethod == UpdateMethods.PolakRibiere)
				{
					tmp[] = gk[] - gkLast[];
					beta = dot(cast(int)tmp.length, cast(double*)gk, 1, cast(double*)tmp, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)gk, 1);
				}
				else if(UpdateMethod == UpdateMethods.DaiYuan)
				{
					tmp[] = gk[] - gkLast[];
					beta = nrm2(cast(int)gk.length, cast(double*)gk, 1)^^2/dot(cast(int)tmp.length, cast(double*)tmp, 1, cast(double*)pkLast, 1);
				}
			}
			else
			{
				beta = 0;
			}
			pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[] + beta*pkLast[];
			//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
			
			tmp[] = xk[] - xkLast[];
			error = nrm2(cast(int)tmp.length, cast(double*)(tmp), 1)/(1 + nrm2(cast(int)xkLast.length, cast(double*)xkLast, 1)) + abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue)/(1+abs(lineResultLast.ObjectiveFunctionValue));
			lineResultLast = lineResult;
			lineSearch.AlphaInitial = 1;
			//lineSearch.AlphaInitial = lineSearch.AlphaInitial*(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//lineSearch.AlphaMax = lineSearch.AlphaInitial+10;
			if(DebugMode)
			{
				writeln("gkLast = ", gkLast, "\tpkLast = ", pkLast);
				writeln("gk = ", gk, "\tpk = ", pk);
				writeln("Error = ", error, "\txk = ", xk);
				Thread.sleep(dur!("msecs")(50));
			}

			iterations++;
			if(FileOutput) ferr.writefln("%d, %40.40f", iterations, error);
		}
		ferr.close();
		f.close();
		result = lineResult;
		result.Iterations = iterations;
		result.MinorIterations = minorIterations;
		return result;
	}

	@property void UpdateMethod(UpdateMethods updateMethod) { mUpdateMethod = updateMethod; }
	@property UpdateMethods UpdateMethod() { return mUpdateMethod; }

	private UpdateMethods mUpdateMethod = UpdateMethods.FletcherReeves;
}