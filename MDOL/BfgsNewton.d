module MDOL.BfgsNewton;

import MDOL.ArrayOps;
import MDOL.BracketAndZoom;
import MDOL.MatrixOps;
import MDOL.ObjectiveFunction;
import MDOL.Optimizer;

import core.thread;

import std.math;
import std.stdio;

import scid.bindings.blas.dblas;

class BfgsNewton : Optimizer
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
		double[] sk = new double[gk.length];
		double[] skLast = new double[gk.length];
		double[] yk = new double[gk.length];
		double[] ykLast = new double[gk.length];
		double[] Vk = new double[gk.length^^2];
		double[] VkLast = new double[gk.length^^2];
		double[] I = IdentityMatrix(cast(int)gk.length);
		//bool converged = false;
		File f;
		if(FileOutput) f = File(PointFilename, "w");
		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");
		double beta;

		double[] Vtmp1 = new double[Vk.length];
		double[] Vtmp2 = new double[Vk.length];
		double[] Vtmp3 = new double[Vk.length];

		Vk[] = I;
		VkLast[] = I;

		if(FileOutput) WriteArrayCSV(f, xk);
		//lineSearch.P.length = pk.length;
		//pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];

		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		Vtmp1[] = -Vk[];
		gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk[], 1, 0, cast(double*)pk, 1);
		pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);
		
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
		//pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);



		yk[] = gk[] - gkLast[];
		sk[] = xk[] - xkLast[];

		double skykDot = dot(cast(int)yk.length, cast(double*)yk, 1, cast(double*)sk, 1);
		double[] skyk = new double[Vk.length];
		double[] yksk = new double[Vk.length];
		double[] sksk = new double[Vk.length];

		gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
		gemm('n', 'n', cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);

		skyk[] *= skykDot^^(-1);
		Vtmp1[] = I[] - skyk[];
		gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)Vk, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);

		gemm('n', 'n', cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
		sksk[] *= skykDot^^(-1);

		Vtmp1[] = I[] - yksk[]*skykDot^^(-1);

		gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);

		Vk[] = Vtmp3[] + sksk[];
		VkLast[] = Vk;

		Vtmp1[] = -Vk[];
		gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
		pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);
		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		lineSearch.AlphaInitial = 1;//(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
		
		double epa = 1.0e-6;
		double epr = 1.0e-6;
		double epg = 1.0e-3;
		
		//while(error >= Tolerance)
		//while( !(abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) > (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && !(nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1) <= epg) )
		while(!converged)
		{

			if(DebugMode) writeln();
			if(FileOutput) WriteArrayCSV(f, xk);

			//gemm('n', 'n', 2, 2, 1, 1, cast(double*)v1, 2, cast(double*)v2, 1, 0, cast(double*)ans, 2);

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

			if(iterations%5)
			//if(true)
			//if(iterations%xk.length)
			{
				//beta = dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)gk, 1)/dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)gkLast, 1);

				yk[] = gk[] - gkLast[];
				sk[] = xk[] - xkLast[];

				skykDot = dot(cast(int)yk.length, cast(double*)yk, 1, cast(double*)sk, 1);

				
				gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
				gemm('n', 'n', cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);
				
				skyk[] *= skykDot^^(-1);
				Vtmp1[] = I[] - skyk[];
				gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)VkLast, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
				
				gemm('n', 'n', cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
				sksk[] *= skykDot^^(-1);

				Vtmp1[] = I[] - yksk[]*skykDot^^(-1);
				
				gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);
				
				Vk[] = Vtmp3[] + sksk[];
				VkLast[] = Vk;
			
			}
			else
			{
				Vk[] = I;
			}
			Vtmp1[] = -Vk[];
			gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
			pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];

			/*
			// Compute new direction at point.
			//pk[] = gk;
			//if(iterations%1.0)
			pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[] + beta*pkLast[];
			*/
			//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);

			if(ErrorConvergence)
			{
				if(error <= Tolerance)
					converged = true;
			}
			else
			{
				//writefln(" current = %20.20f \t previous = %20.20f", lineResult.ObjectiveFunctionValue, lineResultLast.ObjectiveFunctionValue);
				//writefln("%20.20f should be > %20.20f", abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue), (epa + epr*abs(lineResultLast.ObjectiveFunctionValue)));
				//writefln("%20.20f should be < %20.20f", nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1), epg);
				if((abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) < (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && (nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1) <= epg))
					converged = true;
			}

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
				writeln("Error = ", error, "\txk = ", lineResult.DesignVariables);
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
		//writeln("Error = ", error, "\txk = ", result.DesignVariables);
		return result;
	}

	@property void ErrorConvergence(bool errConv) { mErrorConverg = errConv; }
	@property bool ErrorConvergence() { return mErrorConverg; }

	private bool mErrorConverg = false;
}