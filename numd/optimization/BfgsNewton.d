module numd.optimization.BfgsNewton;

import numd.optimization.ArrayOps;
import numd.optimization.BracketAndZoom;
import numd.optimization.MatrixOps;
import numd.optimization.ObjectiveFunction;
import numd.optimization.Optimizer;

import numd.linearalgebra.matrix;

import core.thread;

import std.math;
import std.stdio;

class BfgsNewton(size_t dims) : Optimizer
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

		alias Vector!(dims, double) Vec;
		alias Matrix!(dims, dims, double) Mat;


		auto xk = Vec(InitialGuess[]);
		auto gk = Vec(objectiveFunction.Gradient(InitialGuess));
		auto pkLast = Vec(0);
		auto gkLast = Vec(0);
		auto xkLast = Vec(0);
		auto pk = Vec(0);
		auto tmp = Vec(0);
		auto sk = Vec(0);
		auto skLast = Vec(0);
		auto yk = Vec(0);
		auto ykLast = Vec(0);
		auto Vk = Mat.Identity();
		auto VkLast = Mat.Identity();
		auto I = Mat.Identity();

		File f;
		if(FileOutput) f = File(PointFilename, "w");
		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");
		double beta;

		auto Vtmp1 = Mat(0);
		auto Vtmp2 = Mat(0);
		auto Vtmp3 = Mat(0);

		if(FileOutput) WriteArrayCSV(f, xk);

		Vtmp1 = -Vk;

		pk = (Vtmp1*gk).normalize();
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);

		lineSearch.P = pk.getData();
		lineSearch.InitialGuess = xk.getData();
		lineResult = lineSearch.Optimize(objectiveFunction);
		minorIterations += lineResult.Iterations;
		lineResultLast = lineResult;

		xkLast = xk;
		pkLast = pk;
		gkLast = gk;

		//xk[] = lineResult.DesignVariables;
		xk = lineResult.DesignVariables;
		if(DebugMode) writeln("xk = ", xk);
		
		// Compute gradient at new point.
		//gk[] = objectiveFunction.Gradient(xk[]);
		gk = objectiveFunction.Gradient(xk.getData());

		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);


		yk = gk - gkLast;
		sk = xk - xkLast;

		auto skykDot = yk.dot(sk);

		auto skyk = sk*yk.transpose();
	
		auto yksk = yk*sk.transpose();

		skyk /= skykDot;

		Vtmp1 = I - skyk;

		Vtmp2 = Vtmp1*Vk;

		auto sksk = sk*sk.transpose();

		sksk /= skykDot;
		Vtmp1 = I - yksk/skykDot;

		Vtmp3 = Vtmp2*Vtmp1;

		Vk = Vtmp3 + sksk;
		VkLast = Vk;
		Vtmp1 = -Vk;

		pk = (Vtmp1*gk).normalize();

		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);

		lineSearch.AlphaInitial = 1;//(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
		
		double epa = 1.0e-6;
		double epr = 1.0e-6;
		double epg = 1.0e-3;

		while(!converged)
		{

			if(DebugMode) writeln();
			if(FileOutput) WriteArrayCSV(f, xk);

			lineSearch.P = pk.getData();

			lineSearch.InitialGuess = xk.getData();

			lineResult = lineSearch.Optimize(objectiveFunction);

			minorIterations += lineResult.Iterations;

			xkLast = xk;
			pkLast = pk;
			gkLast = gk;
			xk = lineResult.DesignVariables;

			// Compute gradient at new point.
			gk = objectiveFunction.Gradient(xk.getData());

			if(iterations%5)
			{

				yk = gk - gkLast;
				sk = xk - xkLast;

				skykDot = sk.dot(yk);

				skyk = sk*yk.transpose();

				yksk = yk*sk.transpose();

				skyk /= skykDot;

				Vtmp1 = I - skyk;

				Vtmp2 = Vtmp1*VkLast;

				sksk = sk*sk.transpose();
				sksk /= skykDot;

				Vtmp1 = I - yksk/skykDot;

				Vtmp3 = Vtmp2*Vtmp1;
				//Vtmp2.mult(Vtmp1, Vtmp3);
				Vk = Vtmp3 + sksk;
	
				VkLast = Vk;
			
			}
			else
			{
	
				Vk = Mat.Identity();
			}
	
			Vtmp1 = -Vk;

			pk = (Vtmp1*gk).normalize();

			if(ErrorConvergence)
			{
				if(error <= Tolerance)
					converged = true;
			}
			else
			{
				if((abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) < (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && (gkLast.magnitude() <= epg))
					converged = true;
			}


			tmp = xk - xkLast;
		
			error = tmp.magnitude()/(1 + xkLast.magnitude()) + abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue)/(1+abs(lineResultLast.ObjectiveFunctionValue));
			lineResultLast = lineResult;
			lineSearch.AlphaInitial = 1;

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

version(unittest)
{
	import core.memory;
	import numd.optimization.SecantRoot;
	import std.math;
	import std.complex;

	class Drag : ObjectiveFunction
	{
		this() { Constraints = 1; }
		
		final override Complex!double Compute(Complex!double[] designVar)
		{
			Complex!double Sref = designVar[1];		// m^2
			Complex!double Swet = 2.05*Sref;		// m^2
			Complex!double AR = designVar[0];
			Complex!double b = sqrt(AR*Sref);		// m
			Complex!double c = Sref/b;				// m
			Complex!double Re = (rho*V*c)/mu;
			Complex!double Cf = complex(0);
			
			if(c.re <= 0.85)
			{
				Cf = 1.328/(Re^^0.5);
			}
			else
			{
				Cf = 0.074/(Re^^0.2);
			}
			
			Complex!double Ww = WingWeight(AR, Sref);// N
			Wtot = Ww+W0;
			CL = (2*Wtot)/(rho*V^^2*Sref);
			
			Complex!double CD = 0.03062702/Sref + k*Cf*(Swet/Sref) + CL^^2/(PI*AR*e);
			
			//return CD;
			
			Complex!double LoverD = CL/CD;
			return Wtot/LoverD;
		}
		
		final override Complex!double[] Constraint(Complex!double[] designVar)
		{
			Complex!double[] c = new Complex!double[1];
			auto Vmin = 22;
			auto CLmax = 2;
			auto AR = designVar[0];
			auto Sref = designVar[1];
			Complex!double Ww = WingWeight(AR, Sref);// N
			Wtot = Ww+W0;
			
			c[0] = Sref - (2*Wtot)/(rho*Vmin^^2*CLmax);
			return c;
		}
		
		Complex!double WingWeight(Complex!double AR, Complex!double Sref)
		{
			auto rootFind = new SecantRoot;
			auto span = sqrt(Sref*AR);
			double Nult = 2.5;
			double tc = 0.12;
			auto phi = ((8.71e-5*Nult*span^^3)/(Sref*tc))^^2;
			auto b = -(90.84*Sref + phi*W0);
			auto c = 2062.98*Sref^^2 - phi*W0^^2;
			
			Complex!double WeightDelegate(Complex!double Ww)
			{
				return 45.42*Sref + 8.71e-5*((Nult*span^^3*sqrt(W0*(Ww+W0)))/(Sref*tc)) - Ww;
			}
			
			//return rootFind.Solve(&WeightDelegate, complex(1000.0));
			return (-b + sqrt(b^^2 - 4*c))/2;
		}
		
		private Complex!double CL;
		private Complex!double Wtot;	// N
		private double W0 = 4940;		// N
		private double rho = 1.23;		// kg/m^3
		private double mu = 17.9e-6;	// kg/(m s)
		private double V = 35;			// m/s
		private double k = 1.2;
		private double e = 0.96;
	}

	class Rosenbrock : ObjectiveFunction
	{
		final override Complex!double Compute(Complex!double[] designVars)
		{
			Complex!double sum = complex!(double, double)(0, 0);
			
			for(int i = 0; i < designVars.length-1; i++)
			{
				sum += 100.0*(designVars[i+1] - designVars[i]^^2.0)^^2.0 + (1.0-designVars[i])^^2.0;
				//sum += 100.0*(designVar[i+1] - designVar[i]*designVar[i])*(designVar[i+1] - designVar[i]*designVar[i]) + (1.0-designVar[i])*(1.0-designVar[i]);
			}
			
			/*
			foreach(Complex!double designVar; designVars)
			{
				sum += 100.0*(designVar - designVar^^2.0)^^2.0 + (1.0-designVar)^^2.0;
			}
			*/
			return sum;
		}
		
		final override Complex!double[] Constraint(Complex!double[] designVar)
		{
			Complex!double[] c = new Complex!double[1];
			c[0] = 0;
			return c;
		}
	}

	unittest
	{
		auto bfgs = new BfgsNewton!2;
		
		double[2] initialPoint = [10, 40];
		auto drag = new Drag;
		Result result;
		
		//bfgs.FileOutput = true;
		
		writeln();
		writeln("Starting drag optimization. Start point [35, 10]");
		bfgs.InitialGuess = [35, 10];
		
		writeln();
		//bfgs.DebugMode = true;
		bfgs.ErrorConvergence = true;
		bfgs.PointFilename = "BfgsDragPoints1.csv";
		bfgs.ErrorFilename = "BfgsDragError1.csv";
		result = bfgs.Optimize(drag);
		writeln("BFGS Quasi-Newton");
		writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
		writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
		writeln("\tConverged in ", result.Iterations, " iterations.");
		writeln("\tComputation time: ", result.ComputationTime, " usecs.");
		writeln("\tMinor iterations: ", result.MinorIterations);

		assert(abs(result.ObjectiveFunctionValue - 191.9016258727) < 1.0e-6, "BFGS Quasi-Newton test 2 failed");
	}

	unittest
	{
		auto bfgs = new BfgsNewton!8;
		
		auto rose = new Rosenbrock;
		Result result;
		
		//bfgs.FileOutput = true;
		
		writeln();
		writeln("Starting 8D Rosenbrock.");
		bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
		
		writeln();
		//bfgs.DebugMode = true;
		bfgs.PointFilename = "BfgsRose8DPoints.csv";
		bfgs.ErrorFilename = "BfgsRose8DError.csv";
		result = bfgs.Optimize(rose);
		writeln("BFGS Quasi-Newton");
		writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
		writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
		writeln("\tConverged in ", result.Iterations, " iterations.");
		writeln("\tComputation time: ", result.ComputationTime, " usecs.");
		writeln("\tMinor iterations: ", result.MinorIterations);

		assert(abs(result.ObjectiveFunctionValue) < 1.0e-6, "BFGS Quasi-Newton test 1 failed");
	}

	unittest
	{
		auto bfgs = new BfgsNewton!4;
		
		auto rose = new Rosenbrock;
		Result result;

		//bfgs.FileOutput = true;
		
		writeln();
		writeln("Starting 4D Rosenbrock.");
		//bfgs.DebugMode = true;
		bfgs.Tolerance = 1.0e-10;
		bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
		rose.DerivativeType = "Optimization.FiniteDifference.FiniteDifference";
		//bfgs.ErrorConvergence = true;
		
		writeln();
		//bfgs.DebugMode = true;
		bfgs.PointFilename = "BfgsRose4DPoints.csv";
		bfgs.ErrorFilename = "BfgsRose4DError.csv";
		result = bfgs.Optimize(rose);
		writeln("BFGS Quasi-Newton");
		writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
		writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
		writeln("\tConverged in ", result.Iterations, " iterations.");
		writeln("\tComputation time: ", result.ComputationTime, " usecs.");
		writeln("\tMinor iterations: ", result.MinorIterations);
		assert(true);
	}
}