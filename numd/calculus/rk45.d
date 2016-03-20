module numd.calculus.rk45;

import std.conv;
import std.math;

import numd.calculus.integration;
import numd.linearalgebra.matrix;

unittest
{
	immutable size_t dims = 3;

	//@nogc Vector!dimensions ode(size_t dimensions)(double t, Vector!dimensions x)
	@nogc void ode(size_t dimensions)(ref Vector!dimensions xdot, double t, Vector!dimensions x)
	{
		xdot = [x[1], x[2], -2*x[0]*x[2] - 1 + x[1]^^2];
		//Vector!dimensions xdot = [x[1], x[2], -2*x[0]*x[2] - 1 + x[1]^^2];
		//return xdot;
	}

	Vector!dims x0 = [0, 0, 1.3118];

	auto res = rk45!(ode!dims, dims)(0, 2, x0, 0.005);
}

@nogc IntegrationResults!dimensions rk45(alias func, size_t dimensions, Args...)(double t0, double te, Vector!dimensions x0, double timestep, Args args)
{
	import core.stdc.stdio : printf;

	int points = cast(int)(((te - t0)/timestep).ceil);
	auto results = IntegrationResults!dimensions(points);

	results.X[0][] = x0.mData;

	results.T[0] = t0;
	for(int i = 1; i < results.T.length; i++)
	{
		results.T[i] = results.T[i-1]+timestep;
	}

	alias vec = Vector!dimensions;

	auto k1 = vec(0);
	auto k2 = vec(0);
	auto k3 = vec(0);
	auto k4 = vec(0);
	auto x = vec(0);

	auto tmp = vec(0);

	for(int i = 0; i < (points - 1); i++)
	{
		x.mData[] = results.X[i];

		func(k1, results.T[i], x, args);
		//k1 *= timestep;
		//k1 = timestep*func(results.T[i], x);

		tmp = x + (timestep/2.0)*k1;
		//k2 = timestep*func(results.T[i] + timestep/2, x + k1/2);
		func(k2, results.T[i] + timestep/2.0, tmp, args);
		//k2 *= timestep;

		tmp = x + (timestep/2.0)*k2;
		//k3 = timestep*func(results.T[i] + timestep/2, x + k2/2);
		func(k3, results.T[i] + timestep/2.0, tmp, args);
		//k3 *= timestep;

		tmp = x + timestep*k3;
		//k4 = timestep*func(results.T[i] + timestep, x + k3);
		func(k4, results.T[i] + timestep, tmp, args);
		//k4 *= timestep;

		tmp = x + (timestep/6.0)*(k2 + 2*k2 + 2*k3 + k4);
		results.X[i+1][] = tmp.mData;

/+
		func(k1, results.T[i], x, args);
		k1 *= timestep;

		func(k2, results.T[i] + timestep/2, x + k1/2, args);
		k2 *= timestep;

		func(k3, results.T[i] + timestep/2, x + k2/2, args);
		k3 *= timestep;

		func(k4, results.T[i] + timestep, x + k3, args);
		k4 *= timestep;

		results.X[i+1][] = (x + k2/6 + (k2 + k3)/3 + k4/6).mData;
+/
		//printf("step\n");
	}

	return results;
}
