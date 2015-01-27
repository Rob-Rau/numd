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

import Optimization.Optimizer;
//import Optimization.GoldenSearch;
//import Optimization.BracketAndZoom;
import Optimization.ObjectiveFunction;
//import Optimization.SecantRoot;
import Optimization.NewtonRoot;
//import Optimization.Derivative;
import Optimization.SteepestDescent;
import Optimization.ConjugateGradient;
import Optimization.BfgsNewton;

import scid.bindings.blas.dblas;

class Drag : ObjectiveFunction
{
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
		c[0] = 0;
		return c;
	}
	
	Complex!double WingWeight(Complex!double AR, Complex!double Sref)
	{
		auto rootFind = new NewtonRoot;
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
		
		//return rootFind.Solve(&WeightDelegate, complex(100.0));
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
	final override Complex!double Compute(Complex!double[] designVar)
	{
		Complex!double sum = 0;
		for(int i = 0; i < designVar.length-1; i++)
		{
			sum += 100*(designVar[i+1] - designVar[i]^^2)^^2 + (1-designVar[i])^^2;
		}
		return sum;
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}
}


void DragMinimize();
void DragMinimizeCGFletcherReeves();
void DragCrappyDerivatives();
void Rosenbrock2();
void Rosenbrock2CrappyDerivatives();
void Rosenbrock4();
void Rosenbrock8();
void Rosenbrock16();
void Rosenbrock32();
void Rosenbrock64();
void Rosenbrock128();
void Rosenbrock512();

int main(string[] args)
{
	GC.disable();

	DragMinimize1();
	writeln("================================================================");
	DragMinimize2();
	writeln("================================================================");
	DragMinimize3();
	writeln("================================================================");
	DragCrappyDerivatives();
	writeln("================================================================");
	Rosenbrock2();
	writeln("================================================================");
	Rosenbrock2CrappyDerivatives();
	writeln("================================================================");
	Rosenbrock4();
	writeln("================================================================");
	Rosenbrock8();
	/*
	writeln("================================================================");
	Rosenbrock16();
	writeln("================================================================");
	Rosenbrock32();
	writeln("================================================================");
	Rosenbrock64();
	writeln("================================================================");
	Rosenbrock128();
	writeln("================================================================");
	Rosenbrock512();
	writeln("================================================================");
	*/
	GC.enable();
	return 0;
}

void DragMinimize1()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	writeln();
	writeln("Starting drag optimization. Start point [50, 40]");
	conjGrad.InitialGuess = [50, 40];
	steepest.InitialGuess = [50, 40];
	bfgs.InitialGuess = [50, 40];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsDragPoints1.csv";
	bfgs.ErrorFilename = "BfgsDragError1.csv";
	result = bfgs.Optimize(drag);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradDragPoints1.csv";
	conjGrad.ErrorFilename = "ConjGradDragError1.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [50, 40];
	conjGrad.PointFilename = "ConjGradDragPointsFR1.csv";
	conjGrad.ErrorFilename = "ConjGradDragErrorFR1.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestDragPoints1.csv";
	steepest.ErrorFilename = "SteepestDragError1.csv";
	result = steepest.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void DragMinimize2()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting drag optimization. Start point [20, 40]");
	conjGrad.InitialGuess = [20, 40];
	steepest.InitialGuess = [20, 40];
	bfgs.InitialGuess = [20, 40];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsDragPoints2.csv";
	bfgs.ErrorFilename = "BfgsDragError2.csv";
	result = bfgs.Optimize(drag);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradDragPoints2.csv";
	conjGrad.ErrorFilename = "ConjGradDragError2.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [20, 40];
	conjGrad.PointFilename = "ConjGradDragPointsFR2.csv";
	conjGrad.ErrorFilename = "ConjGradDragErrorFR2.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestDragPoints2.csv";
	steepest.ErrorFilename = "SteepestDragError2.csv";
	result = steepest.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void DragMinimize3()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting drag optimization. Start point [40, 10]");
	conjGrad.InitialGuess = [40, 10];
	steepest.InitialGuess = [40, 10];
	bfgs.InitialGuess = [40, 10];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsDragPoints3.csv";
	bfgs.ErrorFilename = "BfgsDragError3.csv";
	result = bfgs.Optimize(drag);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradDragPoints3.csv";
	conjGrad.ErrorFilename = "ConjGradDragError3.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [40, 10];
	conjGrad.PointFilename = "ConjGradDragPointsFR3.csv";
	conjGrad.ErrorFilename = "ConjGradDragErrorFR3.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestDragPoints3.csv";
	steepest.ErrorFilename = "SteepestDragError3.csv";
	result = steepest.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void DragCrappyDerivatives()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	drag.DerivativeType = "Optimization.FiniteDifference.FiniteDifference";
	drag.StepSize = 1.0e-2;

	writeln();
	writeln("Starting drag optimization with crappy derivatives");
	conjGrad.InitialGuess = [50, 40];
	steepest.InitialGuess = [50, 40];
	bfgs.InitialGuess = [50, 40];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsDragPointsCD.csv";
	bfgs.ErrorFilename = "BfgsDragErrorCD.csv";
	result = bfgs.Optimize(drag);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradDragPointsCD.csv";
	conjGrad.ErrorFilename = "ConjGradDragErrorCD.csv";
	//conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [50, 40];
	conjGrad.PointFilename = "ConjGradDragPointsFR.csv";
	conjGrad.ErrorFilename = "ConjGradDragErrorFR.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(drag);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestDragPointsCD.csv";
	steepest.ErrorFilename = "SteepestDragErrorCD.csv";
	result = steepest.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock2()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;

	auto rose = new Rosenbrock;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	writeln();
	writeln("Starting 2D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0];
	steepest.InitialGuess = [-1.2, 1.0];

	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose2DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose2DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose2DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFR.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFR.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose2DPoints.csv";
	steepest.ErrorFilename = "SteepestRose2DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock2CrappyDerivatives()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	rose.DerivativeType = "Optimization.FiniteDifference.FiniteDifference";
	rose.StepSize = 1.0e-2;

	writeln();
	writeln("Starting 2D Rosenbrock with crappy derivatives.");
	bfgs.InitialGuess = [-1.2, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0];
	steepest.InitialGuess = [-1.2, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose2DPointsCD.csv";
	bfgs.ErrorFilename = "BfgsRose2DErrorCD.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose2DPointsCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose2DPointsCD.csv";
	steepest.ErrorFilename = "SteepestRose2DErrorCD.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock4()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	writeln();
	writeln("Starting 4D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose4DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose4DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose4DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose4DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose4DPoints.csv";
	steepest.ErrorFilename = "SteepestRose4DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock8()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	writeln();
	writeln("Starting 8D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose8DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose8DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose8DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose8DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose8DPoints.csv";
	steepest.ErrorFilename = "SteepestRose8DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock16()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;

	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;

	writeln();
	writeln("Starting 16D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose16DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose16DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose16DPoints.csv";
	steepest.ErrorFilename = "SteepestRose16DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock32()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting 32D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose16DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose16DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose16DPoints.csv";
	steepest.ErrorFilename = "SteepestRose16DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock64()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting 64D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose16DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose16DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose16DPoints.csv";
	steepest.ErrorFilename = "SteepestRose16DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock128()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting 128D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose16DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose16DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose16DPoints.csv";
	steepest.ErrorFilename = "SteepestRose16DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock512()
{
	auto steepest = new SteepestDescent;
	auto conjGrad = new ConjugateGradient;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	conjGrad.FileOutput = true;
	steepest.FileOutput = true;
	
	writeln();
	writeln("Starting 512D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	steepest.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	//writeln("BFGS optimal Point = ", result.DesignVariables, "\t\t\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("BFGS Quasi-Newton");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.PointFilename = "ConjGradRose16DPoints.csv";
	conjGrad.ErrorFilename = "ConjGradRose16DError.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.PolakRibiere;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Polak-Ribiere update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//conjGrad.DebugMode = true;
	conjGrad.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	conjGrad.PointFilename = "ConjGradRose2DPointsFRCD.csv";
	conjGrad.ErrorFilename = "ConjGradRose2DErrorFRCD.csv";
	conjGrad.UpdateMethod = ConjugateGradient.UpdateMethods.FletcherReeves;
	result = conjGrad.Optimize(rose);
	//writeln("Conjugate gradient using Fletcher-Reeves update optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Conjugate gradient using Fletcher-Reeves update");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	steepest.PointFilename = "SteepestRose16DPoints.csv";
	steepest.ErrorFilename = "SteepestRose16DError.csv";
	result = steepest.Optimize(rose);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tFunction value = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, "usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Steepest decent");
	writeln("\tOptimal Point = ", result.DesignVariables);
	writeln("\tDrag = ", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}
