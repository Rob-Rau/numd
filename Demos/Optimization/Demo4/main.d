module main;

import core.memory;
import core.thread;

import std.algorithm;
import std.array;
import std.complex;
import std.file;
import std.math;
import std.stdio;

import plplot;
import Plotting.FigPlot;
import Plotting.Line;
import Plotting.Color;

import MDOL.Complex;
import MDOL.InteriorPenalty;
import MDOL.NewtonRoot;
import MDOL.ObjectiveFunction;
import MDOL.Optimizer;
import MDOL.SecantRoot;
import MDOL.SQP;

import scid.bindings.blas.dblas;
import scid.matrix;
import scid.linalg;

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
		Complex!double Cf = 0.074/(Re^^0.2);
		
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

int main(string[] args)
{
	GC.disable();

	auto interiorPen = new InteriorPenalty;
	auto sqp = new SQP;
	auto drag = new Drag;
	Result result;
	double[2] initialGuess = [48, 52];

	writeln();
	writeln("Initial guess = ", initialGuess);

	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//interiorPen.DebugMode = true;
	//interiorPen.InitialGuess = [8, 5];
	interiorPen.InitialGuess = new double[initialGuess.length];
	interiorPen.InitialGuess[] = initialGuess;
	interiorPen.PointFilename = "InterPenPoints.csv";
	interiorPen.ErrorFilename = "InterPenError.csv";
	interiorPen.FileOutput = false;
	result = interiorPen.Optimize(drag);
	writeln();
	writeln("InteriorPenalty using logarithmic penalty function");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);


	writeln();
	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//sqp.DebugMode = true;
	//sqp.InitialGuess = [8, 5];
	sqp.InitialGuess = new double[initialGuess.length];
	sqp.InitialGuess[] = initialGuess;
	sqp.PointFilename = "SQPpoints.csv";
	sqp.ErrorFilename = "SQPerror.csv";
	sqp.FileOutput = false;
	result = sqp.Optimize(drag);
	writeln();
	writeln("SQP:");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln("-------------------------------------------------------------------------------------------------------------");


	writeln();
	initialGuess = [8, 32];
	writeln("Initial guess = ", initialGuess);

	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//interiorPen.DebugMode = true;
	//interiorPen.InitialGuess = [8, 5];
	interiorPen.InitialGuess = new double[initialGuess.length];
	interiorPen.InitialGuess[] = initialGuess;
	interiorPen.PointFilename = "InterPenPoints.csv";
	interiorPen.ErrorFilename = "InterPenError.csv";
	interiorPen.FileOutput = false;
	result = interiorPen.Optimize(drag);
	writeln();
	writeln("InteriorPenalty using logarithmic penalty function");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	
	writeln();
	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//sqp.DebugMode = true;
	//sqp.InitialGuess = [8, 5];
	sqp.InitialGuess = new double[initialGuess.length];
	sqp.InitialGuess[] = initialGuess;
	sqp.PointFilename = "SQPpoints.csv";
	sqp.ErrorFilename = "SQPerror.csv";
	sqp.FileOutput = false;
	result = sqp.Optimize(drag);
	writeln();
	writeln("SQP:");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln("-------------------------------------------------------------------------------------------------------------");


	writeln();
	initialGuess = [28, 10];
	writeln("Initial guess = ", initialGuess);

	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//interiorPen.DebugMode = true;
	//interiorPen.InitialGuess = [8, 5];
	interiorPen.InitialGuess = new double[initialGuess.length];
	interiorPen.InitialGuess[] = initialGuess;
	interiorPen.PointFilename = "InterPenPoints.csv";
	interiorPen.ErrorFilename = "InterPenError.csv";
	interiorPen.FileOutput = false;
	result = interiorPen.Optimize(drag);
	writeln();
	writeln("InteriorPenalty using logarithmic penalty function");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	
	writeln();
	drag.DerivativeType = "MDOL.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-3;
	//sqp.DebugMode = true;
	//sqp.InitialGuess = [8, 5];
	sqp.InitialGuess = new double[initialGuess.length];
	sqp.InitialGuess[] = initialGuess;
	sqp.PointFilename = "SQPpoints.csv";
	sqp.ErrorFilename = "SQPerror.csv";
	sqp.FileOutput = false;
	result = sqp.Optimize(drag);
	writeln();
	writeln("SQP:");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();

	GC.enable();
	return 0;
}



