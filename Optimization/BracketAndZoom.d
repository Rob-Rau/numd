module MDOL.BracketAndZoom;

import MDOL.Optimizer;
import MDOL.ArrayOps;
import MDOL.ObjectiveFunction;

import core.thread;
import std.array;
import std.algorithm;
import std.complex;
import std.stdio;
import std.math;
import scid.bindings.blas.dblas;

class BracketAndZoom : Optimizer
{

	final override protected Result doOptimize(ObjectiveFunction objectiveFunction)
	{
		Result result;
		version(X86)
		{
			uint iterations = 0;
		}
		else
		{
			ulong iterations = 0;
		}
		bool converged = false;
		double minPoint = 1;
		double minPointOld = 1;
		double minPointErr = 1;
		double minVar;
		double minVarOld = InitialGuess[0];
		double alpha = AlphaInitial;
		double alphaPrev = 0;
		double alphaMax = AlphaMax;
		double[] xk = InitialGuess;
		double[] xkTest = new double[xk.length];
		double[] xHigh = new double[xk.length];
		//writeln("About to cast");
		Complex!double phiInitial = objectiveFunction.Compute(complex(xk));
		//writeln("Finished cast");
		Complex!double phiAlphaPrev = phiInitial;
		double[] g = objectiveFunction.Gradient(xk);
		double[] gInitial = g;
		//double p = -g/abs(g); // we want to head in the opposite direction of the gradient. just direction, no magnitude

		//writeln(DebugMode);
		xkTest[] = xk;
		if(DebugMode) writeln();
		//if(DebugMode) writeln("Pk = ", P);

		//writeln(complex(xk));
		// Zoom function. It's nested in the bracket function to minimize the amount of plumbing between the two functions
		Result Zoom(ObjectiveFunction objectiveFunction, Complex!double initialPoint, Complex!double phiAlphaHigh, Complex!double phiAlphaLow, double alphaHigh, double alphaLow)//, ulong currentIteration, ref Result result)
		{
			if(DebugMode) writeln("Entering zoom");
			while(!converged)
			{
				//double[] xHigh;
				xHigh[] = xk;
				axpy(cast(int)xHigh.length, alphaHigh, cast(double*)P, 1, cast(double*)xHigh, 1);
				//double[] xLow;
				//xLow[] = xk;
				//axpy(cast(int)xLow.length, alphaLow, cast(double*)P, 1, cast(double*)xLow, 1);
				double[] gHigh = objectiveFunction.Gradient(xHigh);
				//double[] gLow = objectiveFunction.Gradient(xLow);
				double a = (phiAlphaHigh.re - alphaHigh*dot(cast(int)P.length, cast(double*)abs(P), 1, cast(double*)gHigh, 1) - phiAlphaLow.re + alphaLow*dot(cast(int)P.length, cast(double*)abs(P), 1, cast(double*)gHigh, 1))/(2*alphaLow*alphaHigh - alphaLow^^2 - alphaHigh^^2);
				double b = dot(cast(int)P.length, cast(double*)abs(P), 1, cast(double*)gHigh, 1) - 2*a*alphaHigh;
				// This is the derivative of the quadratic set to 0 and solved for.
				//double alphaTrial = -b/(2*a);
				//alphaPrev = alpha;
				double alphaTrial = (alphaHigh + alphaLow)/2;
				//writeln();
				//writeln(complex(xk));
				if(DebugMode) writeln("Alpha high = ", alphaHigh, "\t\tAlpha low = ", alphaLow, "\t\tAlpha trial = ", alphaTrial);
				if(DebugMode) writeln("Ø high = ", phiAlphaHigh.re, "\tØ low = ", phiAlphaLow.re);
				//writeln("Xk = ", xkTest);
				//writeln("Gradient = ", g);
				//Thread.sleep(dur!("msecs")( 500 ));

				xkTest[] = xk;
				axpy(cast(int)xk.length, alphaTrial, cast(double*)P, 1, cast(double*)xkTest, 1);

				//writeln("Pk = ", P);
				auto phiTrialPoint = objectiveFunction.Compute(complex(xkTest));
				g = objectiveFunction.Gradient(xkTest);
				//writeln("Gradient = ", g);
				minPoint = phiTrialPoint.re;
				minVar = alphaTrial;


				//minPointErr = abs(minPoint - minPointOld)/(1+abs(minPoint)) + abs(minVar - minVarOld)/(1+abs(minVar));
				minPointErr = abs(phiTrialPoint.re - phiAlphaLow.re)/(1+abs(phiTrialPoint.re)) + abs(alphaTrial - alphaLow)/(1+abs(alphaTrial));
				if(DebugMode) writeln("Error = ", minPointErr, "\t\t\tminPoint = ", minPoint, "\t\tminPointOld = ", minPointOld);
				assert(!isNaN(minPointErr));

				if(iterations >= result.Error.length)
				{
					result.Error.length += 100;
					result.IterArr.length += 100;
				}
				
				result.Error[iterations] = minPointErr;
				result.IterArr[iterations] = iterations;

				if(minPointErr < Tolerance)
				{
					AlphaInitial = alphaTrial;
					result.Iterations = iterations;
					result.DesignVariables = xkTest;
					result.ObjectiveFunctionValue = phiTrialPoint.re;
					return result;
				}

				//writeln(minPointErr);
				minVarOld = minVar;
				minPointOld = minPoint;
				
				//if( (phiTrialPoint.re > (initialPoint.re + mu1*alphaTrial*p*g)) || (phiTrialPoint.re > phiAlphaLow.re))
				if( (phiTrialPoint.re > (initialPoint.re + mu1*alphaTrial*dot(cast(int)P.length, cast(double*)P, 1, cast(double*)gInitial, 1))) || (phiTrialPoint.re > phiAlphaLow.re))
				{
					//writeln("Here");
					alphaHigh = alphaTrial;
					phiAlphaHigh = phiTrialPoint;
				}
				else
				{
					if(abs(dot(cast(int)P.length, cast(double*)P, 1, cast(double*)g, 1)) <= Tolerance*abs(dot(cast(int)P.length, cast(double*)P, 1, cast(double*)gInitial, 1)))
					{
						AlphaInitial = alphaTrial;
						result.Iterations = iterations;
						//result.DesignVariables = array([xk.re + alphaTrial*p]);
						result.DesignVariables = xkTest;
						result.ObjectiveFunctionValue = phiTrialPoint.re;
						return result;
					}
					else if(dot(cast(int)P.length, cast(double*)P, 1, cast(double*)g, 1)*(alphaHigh - alphaLow) >= 0)
					{
						alphaHigh = alphaLow;
						phiAlphaHigh = phiAlphaLow;
					}
					alphaLow = alphaTrial;
					phiAlphaLow = phiTrialPoint;
				}
				iterations++;
			}
			
			return result;
		}
		// End zoom function

		result.Error.length += 100;
		result.IterArr.length += 100;

		while(!converged)
		{
			if(alpha > alphaMax)
				alpha = alphaMax;
			xkTest[] = xk;
			axpy(cast(int)xk.length, alpha, cast(double*)P, 1, cast(double*)xkTest, 1);
			//writeln();
			//if(DebugMode) writeln("Xk = ", xkTest);
			//Complex!double phiAlphaI = objectiveFunction.Compute(array([xk + complex(alpha*p, h)]));
			Complex!double phiAlphaI = objectiveFunction.Compute(complex(xkTest));

			if(DebugMode) writeln("Alpha = ", alpha, "\tØ = ", phiAlphaI.re);
			if(!isInfinity(phiAlphaI.re))
			{
				g = objectiveFunction.Gradient(xkTest);
				//writeln("Gradient = ", g);
				minPoint = phiAlphaI.re;
				minVar = alpha;

				minPointErr = abs(minPoint - minPointOld)/(1+abs(minPoint)) + abs(minVar - minVarOld)/(1+abs(minVar));
				if(iterations >= result.Error.length)
				{
					result.Error.length += 100;
					result.IterArr.length += 100;
				}
				
				result.Error[iterations] = minPointErr;
				result.IterArr[iterations] = iterations;

				if(minPointErr < Tolerance)
				{
					AlphaInitial = alpha;
					result.Iterations = iterations;
					//result.DesignVariables = array([xk.re + alpha*p]);
					result.DesignVariables = xkTest;
					result.ObjectiveFunctionValue = phiAlphaI.re;
					return result;
				}

				//writeln(minPointErr);
				minVarOld = minVar;
				minPointOld = minPoint;

				//if( (phiAlphaI.re > (phiInitial.re + mu1*alpha*g*p)) || ((phiAlphaI.re > phiAlphaPrev.re) && iterations > 0) )
				if( (phiAlphaI.re > (phiInitial.re + mu1*alpha*dot(cast(int)P.length, cast(double*)P, 1, cast(double*)g, 1))) || ((phiAlphaI.re > phiAlphaPrev.re) && iterations > 0) )
				{
					AlphaInitial = alpha;
					if(DebugMode) writeln("Alpha = ", alpha, "\t alphaPrev = ", alphaPrev, "\tØ = ", phiAlphaI.re);
					result = Zoom(objectiveFunction, phiInitial, phiAlphaI, phiAlphaPrev, alpha, alphaPrev);
					//result.Iterations += iterations;
					return result;
				}

				//if(abs((phiAlphaI.im/h)*p) <= Tolerance*abs((phiInitial.im)/h*p))
				if(abs(dot(cast(int)P.length, cast(double*)P, 1, cast(double*)g, 1)) <= Tolerance*abs(dot(cast(int)P.length, cast(double*)P, 1, cast(double*)gInitial, 1)))
				{
					AlphaInitial = alpha;
					result.Iterations = iterations;
					//result.DesignVariables = array([(xk + complex(alpha*p, h)).re]);
					result.DesignVariables = xkTest;
					result.ObjectiveFunctionValue = phiAlphaI.re;
					return result;
				}
				//else if( (phiAlphaI.im/h)*p >= 0)
				else if( dot(cast(int)P.length, cast(double*)P, 1, cast(double*)g, 1) >= 0)
				{
					AlphaInitial = alpha;
					//writeln("alphaPrev = ", alphaPrev, "\tAlpha = ", alpha);
					//writeln();
					result = Zoom(objectiveFunction, phiInitial, phiAlphaPrev, phiAlphaI, alphaPrev, alpha);
					//result.Iterations = iterations;
					return result;
				}
				else
				{
					//alphaPrev = alpha;
					if(DebugMode) writeln("alpha = ", alpha, "\tPhi = ", phiAlphaI.re);
					alpha = (alpha + alphaMax)/2;

					if(alpha > alphaMax)
						alpha = alphaMax;

					if(abs(alpha - alphaMax) < Tolerance)
					{
						AlphaInitial = alpha;
						result.Iterations = iterations;
						//result.DesignVariables = array([(xk + complex(alpha*p, h)).re]);
						result.DesignVariables = xkTest;
						result.ObjectiveFunctionValue = phiAlphaI.re;
						return result;
					}
				}
			}
			else
			{
				alpha = 0.5*alpha;
			}
			iterations++;
		}

		return result;
	}
	
	@property void AlphaInitial(double alpha) { mAlphaInitial = alpha; }
	@property double AlphaInitial() { return mAlphaInitial; }

	@property void AlphaMax(double alpha) { mAlphaMax = alpha; }
	@property double AlphaMax() { return mAlphaMax; }

	@property void P(double[] p) { mP = p; }
	@property double[] P() { return mP; }

	private double[] mP = [0];
	private double mAlphaMax = 10;
	private double mAlphaInitial = 1;
	private double mu1 = 10e-2;
}