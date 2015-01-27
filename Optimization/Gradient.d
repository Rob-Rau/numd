module Optimization.Gradient;

interface IGradient
{
	double[] Gradient(double[] point);
}