module Optimization.Complex;

import std.complex;
import std.math;

Complex!double log(Complex!double x)
{
	return complex!double(std.math.log(x.re), 1/x.re*x.im);
}