module Optimization.Gradient;

interface IGradient
{
	double[] Gradient(in double[] point);
}