module main;

import core.memory;
import core.thread;

import std.algorithm;
import std.array;
import std.complex;
import std.file;
import std.math;
import std.stdio;

/*import plplot;
import Plotting.FigPlot;
import Plotting.Line;
import Plotting.Color;
*/
import Optimization.BfgsNewton;
import Optimization.Complex;
import Optimization.NewtonRoot;
import Optimization.ObjectiveFunction;
import Optimization.Optimizer;
import Optimization.ParticleSwarmOptimizer;
import Optimization.SecantRoot;

import cblas;
//import scid.bindings.blas.dblas;
//import scid.matrix;
//import scid.linalg;

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

int main(string[] args)
{
	GC.disable();
	/*
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	//double[2] initialPoint = [10, 40];
	auto rose = new Rosenbrock;
	auto drag = new Drag;
	Result result;

	pso.Bounds = [[-1, 1], [-1, 1]];
	pso.Particles = 20;
	pso.Dimensions = 2;
	//pso.DebugMode = true;


	bfgs.InitialGuess = [35, 15];
	bfgs.PointFilename = "BfgsDragPoints.csv";
	bfgs.ErrorFilename = "BfgsDragError.csv";
	bfgs.FileOutput = true;
	//bfgs.DebugMode = true;
	result = bfgs.Optimize(drag);
	writeln();
	writeln("bfgs");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tOptimal Point Value = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "  usecs.");
	//writeln("\tMinor iterations: ", result.MinorIterations);

	pso.Particles = 2;
	pso.DebugMode = true;
	pso.Bounds = [[10, 60], [10, 60]];
	pso.PointFilename = "ParticleSwarmPoints";
	pso.ErrorFilename = "ParticleSwarmError.csv";
	pso.FileOutput = true;
	result = pso.Optimize(drag);
	writeln();
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, "  usecs.");
	*/
	DragMinimize1();
	writeln("================================================================");
	DragMinimize2();
	writeln("================================================================");
	Rosenbrock2();
	writeln("================================================================");
	Rosenbrock4();
	writeln("================================================================");
	Rosenbrock8();
	writeln("================================================================");
	/*
	Rosenbrock16();
	writeln("================================================================");
	Rosenbrock32();
	writeln("================================================================");
	Rosenbrock64();
	writeln("================================================================");
	Rosenbrock128();
	writeln("================================================================");
	*/
	//Rosenbrock512();
	//writeln("================================================================");


	GC.enable();
	return 0;
}

void DragMinimize1()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;
	
	//bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting drag optimization. Start point [35, 10]");
	bfgs.InitialGuess = [35, 10];
	
	writeln();
	bfgs.DebugMode = false;
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

	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[10, 60], [10, 60]];
	pso.PointFilename = "ParticleSwarmDrag1Points";
	pso.ErrorFilename = "ParticleSwarmDrag1Error.csv";
	//pso.FileOutput = true;
	pso.Dimensions = 2;
	writeln();
	//steepest.DebugMode = true;
	result = pso.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, " usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void DragMinimize2()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	double[2] initialPoint = [10, 40];
	auto drag = new Drag;
	Result result;
	
	//bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting drag optimization. Start point [15, 40]");
	bfgs.Tolerance = 1.0e-10;
	bfgs.InitialGuess = [15, 40];
	
	writeln();
	bfgs.ErrorConvergence = true;
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsDragPoints2.csv";
	bfgs.ErrorFilename = "BfgsDragError2.csv";
	result = bfgs.Optimize(drag);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[10, 60], [10, 60]];
	//pso.FileOutput = true;
	pso.Dimensions = 2;
	pso.PointFilename = "ParticleSwarmDrag2Points";
	pso.ErrorFilename = "ParticleSwarmDrag2Error.csv";
	result = pso.Optimize(drag);
	//writeln("Steepest decent optimal Point = ", result.DesignVariables, "\tDrag = ", result.ObjectiveFunctionValue, "\tConverged in ", result.Iterations, " iterations.\tComputation time: ", result.ComputationTime, " usecs.\t Minor iterations: ", result.MinorIterations);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock2()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	//bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 2D Rosenbrock.");
	bfgs.Tolerance = 1.0e-10;
	bfgs.InitialGuess = [-1.2, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose2DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose2DError.csv";
	result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Tolerance = 1.0e-10;
	pso.Bounds = [[-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 2;
	pso.PointFilename = "ParticleSwarmRose2Doints";
	pso.ErrorFilename = "ParticleSwarmRose2DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock4()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	//bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 4D Rosenbrock.");
	//bfgs.DebugMode = true;
	bfgs.Tolerance = 1.0e-10;
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0];
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
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Tolerance = 1.0e-10;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 4;
	pso.PointFilename = "ParticleSwarmRose4Doints";
	pso.ErrorFilename = "ParticleSwarmRose4DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock8()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
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

	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Tolerance = 1.0e-10;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 8;
	pso.PointFilename = "ParticleSwarmRose8Doints";
	pso.ErrorFilename = "ParticleSwarmRose8DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock16()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 16D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose16DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose16DError.csv";
	result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 16;
	pso.PointFilename = "ParticleSwarmRose16Doints";
	pso.ErrorFilename = "ParticleSwarmRose16DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock32()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 32D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose32DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose32DError.csv";
	result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 32;
	pso.PointFilename = "ParticleSwarmRose32Doints";
	pso.ErrorFilename = "ParticleSwarmRose32DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock64()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 64D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose64DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose64DError.csv";
	//result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	//pso.FileOutput = true;
	pso.Dimensions = 64;
	pso.PointFilename = "ParticleSwarmRose64Doints";
	pso.ErrorFilename = "ParticleSwarmRose64DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock128()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 128D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose128DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose128DError.csv";
	//result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	
	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	pso.FileOutput = true;
	pso.Dimensions = 128;
	pso.PointFilename = "ParticleSwarmRose128Doints";
	pso.ErrorFilename = "ParticleSwarmRose128DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}

void Rosenbrock512()
{
	auto pso = new ParticleSwarmOptimizer;
	auto bfgs = new BfgsNewton;
	
	auto rose = new Rosenbrock;
	Result result;
	
	bfgs.FileOutput = true;
	
	writeln();
	writeln("Starting 512D Rosenbrock.");
	bfgs.InitialGuess = [-1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0, -1.2, 1.0, 1.0, 1.0];
	
	writeln();
	//bfgs.DebugMode = true;
	bfgs.PointFilename = "BfgsRose512DPoints.csv";
	bfgs.ErrorFilename = "BfgsRose512DError.csv";
	//result = bfgs.Optimize(rose);
	writeln("BFGS Quasi-Newton");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);

	writeln();
	//steepest.DebugMode = true;
	pso.Particles = 40;
	//pso.DebugMode = true;
	pso.Bounds = [[-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10], [-10, 10]];
	pso.FileOutput = true;
	pso.Dimensions = 512;
	pso.PointFilename = "ParticleSwarmRose512Doints";
	pso.ErrorFilename = "ParticleSwarmRose512DError.csv";
	result = pso.Optimize(rose);
	writeln("Particle Swarm Optimizer");
	writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
	writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
	writeln("\tConverged in ", result.Iterations, " iterations.");
	writeln("\tComputation time: ", result.ComputationTime, " usecs.");
	writeln("\tMinor iterations: ", result.MinorIterations);
	writeln();
}




