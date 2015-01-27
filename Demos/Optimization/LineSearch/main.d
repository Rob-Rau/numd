module main;

import core.thread;
import std.stdio;
import std.math;
import std.array;
import std.complex;

import plplot;
import Plotting.FigPlot;
import Plotting.Line;
import Plotting.Color;

import Optimization.Optimizer;
import Optimization.GoldenSearch;
import Optimization.BracketAndZoom;
import Optimization.ObjectiveFunction;

void TestLine2D();
void TestCountour();

void test1();
void test2();
void test3();
void test4();

class Drag : ObjectiveFunction
{
	final override Complex!double Compute(Complex!double[] designVar)
	{
		double rho = 1.23;		// kg/m^3
		double mu = 17.9e-6;	// kg/(m s)
		double V = 35;			// m/s
		double Sref = 11.8;		// m^2
		double Swet = 2.05*Sref;// m^2
		double k = 1.2;
		double CL = 0.3;
		double e = 0.96;
		Complex!double AR = designVar[0];
		Complex!double b = sqrt(AR*Sref);
		Complex!double c = Sref/b;
		Complex!double Re = (rho*V*c)/mu;
		Complex!double Cf = 0.074/(Re^^0.2);

		Complex!double CD = k*Cf*(Swet/Sref) + CL^^2/(PI*AR*e);

		return CD;
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}

}

class TestFunc : ObjectiveFunction
{
	final override Complex!double Compute(Complex!double[] designVar)
	{
		return (designVar[0] - 3)*designVar[0]^^3*(designVar[0]-6)^^4;
	}
	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}

}

int main(char[][] args)
{
	BracketAndZoom bazOptimizer = new BracketAndZoom;
	Result bazResult;
	auto drag = new Drag;

	bazOptimizer.InitialGuess = 30;
	bazOptimizer.AlphaMax = 10;
	bazOptimizer.Tolerance = 1e-5;
	bazOptimizer.P = [-1];
	/*
	writeln();
	test1();
	writeln();
	test2();
	writeln();
	test3();
	writeln();
	test4();
	*/
	bazResult = bazOptimizer.Optimize(drag);
	writeln("Bracket and zoom converged in ", bazResult.Iterations, " iterations.\tOptimized aspect ratio = ", bazResult.DesignVariables[0], "\tOptimized drag coefficient = ", bazResult.ObjectiveFunctionValue);
	//writeln("Average bracket and zoom iteration time: ", bazIterTime, " usecs/iteration");
	return 0;
}

void test1()
{
	auto goldenOptimizer = new GoldenSearch;
	auto bazOptimizer = new BracketAndZoom;
	auto drag = new Drag;
	auto test = new TestFunc;
	double[] bazIters;
	double[] goldenIters;
	Result goldenResult;
	Result bazResult;
	double goldenCompTime = 0;
	double bazCompTime = 0;
	double bazIterTime = 0;
	double goldenIterTime = 0;
	int runTimes = 10;
	Figure fig = new Figure;
	Line2D plot1 = new Line2D;
	

	version(Win32)
	{
		fig.Driver = "svg";
		fig.Filename = "Convergance_1_30.svg";
	}
	else
	{
		fig.Driver = "epsqt";
		fig.Filename = "Convergance_1_30.eps";
	}
	
	goldenOptimizer.InitialGuess = 1;
	goldenOptimizer.Interval = 30;
	goldenOptimizer.Tolerance = 1e-3;
	goldenResult = goldenOptimizer.Optimize(drag);
	writeln("Starting Golden section search algo");
	for(int i = 0; i < runTimes; i++)
	{
		goldenResult = goldenOptimizer.Optimize(drag);
		goldenCompTime += cast(double)goldenResult.ComputationTime;
	}
	goldenCompTime /= runTimes;
	goldenIterTime = cast(double)goldenResult.ComputationTime/cast(double)goldenResult.Iterations;
	
	bazOptimizer.InitialGuess = 1;
	bazOptimizer.AlphaMax = 30;
	bazOptimizer.Tolerance = 1e-4;
	
	writeln("Starting Bracket and Zoom algo");
	for(int i = 0; i < runTimes; i++)
	{
		bazResult = bazOptimizer.Optimize(drag);
		bazCompTime += cast(double)bazResult.ComputationTime;
	}
	bazCompTime /= runTimes;
	bazIterTime = cast(double)bazResult.ComputationTime/cast(double)bazResult.Iterations;
	
	goldenIters.length = goldenResult.Iterations;
	bazIters.length = bazResult.Iterations;
	
	for(int i = 0; i < bazIters.length; i++)
	{
		bazIters[i] = i;
	}
	
	for(int i = 0; i < goldenIters.length; i++)
	{
		goldenIters[i] = i;
	}
	
	writeln();
	writeln("Golden section search converged in ", goldenResult.Iterations, " iterations with an average computation time of ", goldenCompTime, " usecs.\tOptimized aspect ratio = ", goldenResult.DesignVariables[0], "\tOptimized drag coefficient = ", goldenResult.ObjectiveFunctionValue);
	writeln("Average golden search iteration time: ", goldenIterTime, " usecs/iteration"); 
	writeln();
	writeln("Bracket and zoom converged in ", bazResult.Iterations, " iterations with an average computation time of ", bazCompTime, " usecs.\tOptimized aspect ratio = ", bazResult.DesignVariables[0], "\tOptimized drag coefficient = ", bazResult.ObjectiveFunctionValue);
	writeln("Average bracket and zoom iteration time: ", bazIterTime, " usecs/iteration");
	writeln();
	
	plot1.AddDataSet(bazResult.IterArr, bazResult.Error);
	plot1.AddDataSet(goldenResult.IterArr, goldenResult.Error);
	plot1.AddLineColor(Blue);
	plot1.AddLineColor(Red);
	plot1.Label("Iterations", "Error");
	plot1.LogY();
	fig.AddPlot(plot1);
	fig.Show();	
}

void test2()
{
	auto goldenOptimizer = new GoldenSearch;
	auto bazOptimizer = new BracketAndZoom;
	auto drag = new Drag;
	auto test = new TestFunc;
	double[] bazIters;
	double[] goldenIters;
	Result goldenResult;
	Result bazResult;
	double goldenCompTime = 0;
	double bazCompTime = 0;
	double bazIterTime = 0;
	double goldenIterTime = 0;
	int runTimes = 10;
	Figure fig = new Figure;
	Line2D plot1 = new Line2D;
	
	
	version(Win32)
	{
		fig.Driver = "svg";
		fig.Filename = "Convergance_30_10.svg";
	}
	else
	{
		fig.Driver = "epsqt";
		fig.Filename = "Convergance_30_10.eps";
	}
	
	goldenOptimizer.InitialGuess = 30;
	goldenOptimizer.Interval = 10;
	goldenOptimizer.Tolerance = 1e-3;
	goldenResult = goldenOptimizer.Optimize(drag);
	writeln("Starting Golden section search algo");
	for(int i = 0; i < runTimes; i++)
	{
		goldenResult = goldenOptimizer.Optimize(drag);
		goldenCompTime += cast(double)goldenResult.ComputationTime;
	}
	goldenCompTime /= runTimes;
	goldenIterTime = cast(double)goldenResult.ComputationTime/cast(double)goldenResult.Iterations;
	
	bazOptimizer.InitialGuess = 30;
	bazOptimizer.AlphaMax = 10;
	bazOptimizer.Tolerance = 1e-4;
	
	writeln("Starting Bracket and Zoom algo");
	for(int i = 0; i < runTimes; i++)
	{
		bazResult = bazOptimizer.Optimize(drag);
		bazCompTime += cast(double)bazResult.ComputationTime;
	}
	bazCompTime /= runTimes;
	bazIterTime = cast(double)bazResult.ComputationTime/cast(double)bazResult.Iterations;
	
	goldenIters.length = goldenResult.Iterations;
	bazIters.length = bazResult.Iterations;
	
	for(int i = 0; i < bazIters.length; i++)
	{
		bazIters[i] = i;
	}
	
	for(int i = 0; i < goldenIters.length; i++)
	{
		goldenIters[i] = i;
	}
	
	writeln();
	writeln("Golden section search converged in ", goldenResult.Iterations, " iterations with an average computation time of ", goldenCompTime, " usecs.\tOptimized aspect ratio = ", goldenResult.DesignVariables[0], "\tOptimized drag coefficient = ", goldenResult.ObjectiveFunctionValue);
	writeln("Average golden search iteration time: ", goldenIterTime, " usecs/iteration"); 
	writeln();
	writeln("Bracket and zoom converged in ", bazResult.Iterations, " iterations with an average computation time of ", bazCompTime, " usecs.\tOptimized aspect ratio = ", bazResult.DesignVariables[0], "\tOptimized drag coefficient = ", bazResult.ObjectiveFunctionValue);
	writeln("Average bracket and zoom iteration time: ", bazIterTime, " usecs/iteration");
	writeln();
	
	plot1.AddDataSet(bazResult.IterArr, bazResult.Error);
	plot1.AddDataSet(goldenResult.IterArr, goldenResult.Error);
	plot1.AddLineColor(Blue);
	plot1.AddLineColor(Red);
	plot1.Label("Iterations", "Error");
	plot1.LogY();
	fig.AddPlot(plot1);
	fig.Show();
}

void test3()
{
	auto goldenOptimizer = new GoldenSearch;
	auto bazOptimizer = new BracketAndZoom;
	auto drag = new Drag;
	auto test = new TestFunc;
	double[] bazIters;
	double[] goldenIters;
	Result goldenResult;
	Result bazResult;
	double goldenCompTime = 0;
	double bazCompTime = 0;
	double bazIterTime = 0;
	double goldenIterTime = 0;
	int runTimes = 10;
	Figure fig = new Figure;
	Line2D plot1 = new Line2D;
	
	version(Win32)
	{
		fig.Driver = "svg";
		fig.Filename = "Convergance_40_90.svg";
	}
	else
	{
		fig.Driver = "epsqt";
		fig.Filename = "Convergance_40_90.eps";
	}
	
	goldenOptimizer.InitialGuess = 40;
	goldenOptimizer.Interval = 90;
	goldenOptimizer.Tolerance = 1e-3;
	goldenResult = goldenOptimizer.Optimize(drag);
	writeln("Starting Golden section search algo");
	for(int i = 0; i < runTimes; i++)
	{
		goldenResult = goldenOptimizer.Optimize(drag);
		goldenCompTime += cast(double)goldenResult.ComputationTime;
	}
	goldenCompTime /= runTimes;
	goldenIterTime = cast(double)goldenResult.ComputationTime/cast(double)goldenResult.Iterations;
	
	bazOptimizer.InitialGuess = 40;
	bazOptimizer.AlphaMax = 90;
	bazOptimizer.Tolerance = 1e-4;
	
	writeln("Starting Bracket and Zoom algo");
	for(int i = 0; i < runTimes; i++)
	{
		bazResult = bazOptimizer.Optimize(drag);
		bazCompTime += cast(double)bazResult.ComputationTime;
	}
	bazCompTime /= runTimes;
	bazIterTime = cast(double)bazResult.ComputationTime/cast(double)bazResult.Iterations;
	
	goldenIters.length = goldenResult.Iterations;
	bazIters.length = bazResult.Iterations;
	
	for(int i = 0; i < bazIters.length; i++)
	{
		bazIters[i] = i;
	}
	
	for(int i = 0; i < goldenIters.length; i++)
	{
		goldenIters[i] = i;
	}
	
	writeln();
	writeln("Golden section search converged in ", goldenResult.Iterations, " iterations with an average computation time of ", goldenCompTime, " usecs.\tOptimized aspect ratio = ", goldenResult.DesignVariables[0], "\tOptimized drag coefficient = ", goldenResult.ObjectiveFunctionValue);
	writeln("Average golden search iteration time: ", goldenIterTime, " usecs/iteration"); 
	writeln();
	writeln("Bracket and zoom converged in ", bazResult.Iterations, " iterations with an average computation time of ", bazCompTime, " usecs.\tOptimized aspect ratio = ", bazResult.DesignVariables[0], "\tOptimized drag coefficient = ", bazResult.ObjectiveFunctionValue);
	writeln("Average bracket and zoom iteration time: ", bazIterTime, " usecs/iteration");
	writeln();
	
	plot1.AddDataSet(bazResult.IterArr, bazResult.Error);
	plot1.AddDataSet(goldenResult.IterArr, goldenResult.Error);
	plot1.AddLineColor(Blue);
	plot1.AddLineColor(Red);
	plot1.Label("Iterations", "Error");
	plot1.LogY();
	fig.AddPlot(plot1);
	fig.Show();	
}

void test4()
{
	auto goldenOptimizer = new GoldenSearch;
	auto bazOptimizer = new BracketAndZoom;
	auto drag = new Drag;
	auto test = new TestFunc;
	double[] bazIters;
	double[] goldenIters;
	Result goldenResult;
	Result bazResult;
	double goldenCompTime = 0;
	double bazCompTime = 0;
	double bazIterTime = 0;
	double goldenIterTime = 0;
	int runTimes = 10;
	Figure fig = new Figure;
	Line2D plot1 = new Line2D;
	
	
	version(Win32)
	{
		fig.Driver = "svg";
		fig.Filename = "Convergance_20_90.svg";
	}
	else
	{
		fig.Driver = "epsqt";
		fig.Filename = "Convergance_20_90.eps";
	}
	
	goldenOptimizer.InitialGuess = 20;
	goldenOptimizer.Interval = 90;
	goldenOptimizer.Tolerance = 1e-3;
	goldenResult = goldenOptimizer.Optimize(drag);
	writeln("Starting Golden section search algo");
	for(int i = 0; i < runTimes; i++)
	{
		goldenResult = goldenOptimizer.Optimize(drag);
		goldenCompTime += cast(double)goldenResult.ComputationTime;
	}
	goldenCompTime /= runTimes;
	goldenIterTime = cast(double)goldenResult.ComputationTime/cast(double)goldenResult.Iterations;
	
	bazOptimizer.InitialGuess = 20;
	bazOptimizer.AlphaMax = 90;
	bazOptimizer.Tolerance = 1e-4;
	
	writeln("Starting Bracket and Zoom algo");
	for(int i = 0; i < runTimes; i++)
	{
		bazResult = bazOptimizer.Optimize(drag);
		bazCompTime += cast(double)bazResult.ComputationTime;
	}
	bazCompTime /= runTimes;
	bazIterTime = cast(double)bazResult.ComputationTime/cast(double)bazResult.Iterations;
	
	goldenIters.length = goldenResult.Iterations;
	bazIters.length = bazResult.Iterations;
	
	for(int i = 0; i < bazIters.length; i++)
	{
		bazIters[i] = i;
	}
	
	for(int i = 0; i < goldenIters.length; i++)
	{
		goldenIters[i] = i;
	}
	
	writeln();
	writeln("Golden section search converged in ", goldenResult.Iterations, " iterations with an average computation time of ", goldenCompTime, " usecs.\tOptimized aspect ratio = ", goldenResult.DesignVariables[0], "\tOptimized drag coefficient = ", goldenResult.ObjectiveFunctionValue);
	writeln("Average golden search iteration time: ", goldenIterTime, " usecs/iteration"); 
	writeln();
	writeln("Bracket and zoom converged in ", bazResult.Iterations, " iterations with an average computation time of ", bazCompTime, " usecs.\tOptimized aspect ratio = ", bazResult.DesignVariables[0], "\tOptimized drag coefficient = ", bazResult.ObjectiveFunctionValue);
	writeln("Average bracket and zoom iteration time: ", bazIterTime, " usecs/iteration");
	writeln();
	
	plot1.AddDataSet(bazResult.IterArr, bazResult.Error);
	plot1.AddDataSet(goldenResult.IterArr, goldenResult.Error);
	plot1.AddLineColor(Blue);
	plot1.AddLineColor(Red);
	plot1.Label("Iterations", "Error");
	plot1.LogY();
	fig.AddPlot(plot1);
	fig.Show();
}