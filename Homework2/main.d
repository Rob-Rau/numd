module main;

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

import MDOL.Optimizer;
import MDOL.GoldenSearch;
import MDOL.BracketAndZoom;
import MDOL.ObjectiveFunction;
import MDOL.SecantRoot;
import MDOL.NewtonRoot;
import MDOL.Derivative;

double PlotDerivativeError();

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
		
		return CD;
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

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}

	Complex!double Drag(Complex!double AR, Complex!double Sref)
	{
		Complex!double CD = Compute(array([AR, Sref]));
		Complex!double LoverD = CL/CD;
		return Wtot/LoverD;
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
	auto drag = new Drag;

	//writeln(drag.WingWeight(complex(5.0), complex(20.0)));

	double bestStep = PlotDerivativeError();


	drag.StepSize = 1.0e-42;
	writeln();
	writeln("About to compute gradient with the complex step method...");
	double[] grad = drag.Gradient(array([10.0, 20.0]));
	writefln("The gradient is %(%20.20f, %).", grad);

	writeln();
	drag.DerivativeType = "FiniteDifference.FiniteDifference";
	drag.StepSize = bestStep;
	writeln("About to compute gradient with the finite difference method (step size = ", drag.StepSize, ")...");
	grad = drag.Gradient(array([10.0, 20.0]));
	writefln("The gradient is %(%20.20f, %).", grad);
	writeln();
	//writeln(drag.WingWeight(complex(5.0), complex(20.0)));

	return 0;
}

double PlotDerivativeError()
{
	Figure fig = new Figure;
	Line2D plot1 = new Line2D;
	auto drag = new Drag;
	drag.StepSize = 1.0e-42;
	double[] baseline = drag.Gradient(array([10.0, 20.0]));
	//writeln("Completed baseline");

	double[] h, err, err1;
	h.length = 30;
	err.length = 30;
	err1.length = 30;

	h[0] = 1.0e-1;
	for(int i = 1; i < h.length; i++)
	{
		h[i] = h[i-1]/10;
	}
	//writeln("made h array");
	drag.DerivativeType = "FiniteDifference.FiniteDifference";

	for(int i = 0; i < h.length; i++)
	{
		drag.StepSize = h[i];

		drag.DerivativeType = "FiniteDifference.FiniteDifference";
		double[] grad = drag.Gradient(array([10.0, 20.0]));
		err[i] = abs((baseline[0] - grad[0])/baseline[0]);

		drag.DerivativeType = "ComplexStep.ComplexStep";
		grad = drag.Gradient(array([10.0, 20.0]));
		err1[i] = abs((baseline[0] - grad[0])/baseline[0]);
	}

	auto minErr = minPos(err)[0];
	int minIndex;

	for(int i = 0; i < err.length; i++)
	{
		if(err[i] == minErr)
		{
			minIndex = i;
			break;
		}
	}
	writeln();
	writeln("Minimum error = ", minErr, " with a step size of ", h[minIndex]);

	auto f = File("error.csv", "w");
	for(int i = 0; i < err.length; i++)
	{
		f.writefln("%40.40f, %40.40f, %40.40f", h[i], err[i], err1[i]);
	}
	/*
	plot1.LogLog();
	plot1.AddDataSet(h, err);
	plot1.AddDataSet(h, err1);
	plot1.AddLineColor(Blue);
	plot1.AddLineColor(Red);
	fig.AddPlot(plot1);
	fig.Show();
	*/
	return h[minIndex];
}












