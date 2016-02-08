module Optimization.SQP;

import Optimization.ArrayOps;
import Optimization.BracketAndZoom;
import Optimization.Complex;
import Optimization.InteriorPenalty;
import Optimization.MatrixOps;
import Optimization.ObjectiveFunction;
import Optimization.Optimizer;

import core.thread;

import std.algorithm;
import std.complex;
import std.math;
import std.stdio;

//import scid.bindings.lapack.dlapack;
//import scid.bindings.blas.dblas;
import cblas;
import scid.matrix;
import scid.linalg;

class MeritFunction : InteriorPenaltyFunction
{
	final override Complex!double Compute(Complex!double[] designVar)
	{
		Complex!double cNorm = complex(0, 0);
		Complex!double[] c = ObjectiveFunc.Constraint(designVar);
		for(int i = 0; i < c.length; i++)
		{
			cNorm += abs(c[i]);
		}
		return ObjectiveFunc.Compute(designVar) + (1/Mu)*cNorm;
	}
	
	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}
}

class SQP : Optimizer
{
	
	final override protected Result doOptimize(ObjectiveFunction objectiveFunction)
	{
		version(X86)
		{
			uint iterations = 0;
			uint minorIterations = 0;
		}
		else
		{
			ulong iterations = 0;
			ulong minorIterations = 0;
		}
		ObjectiveFunc = objectiveFunction;
		//auto lineSearch = new BracketAndZoom;
		//lineSearch.DebugMode = DebugMode;
		Result lineResult;
		Result lineResultLast = lineResult;
		Result result;
		MeritFunction Ø = new MeritFunction;
		Ø.ObjectiveFunc = objectiveFunction;
		int constraints = ObjectiveFunc.Constraints;
		bool converged = false;
		double error = 1;
		double delta = 0.01;
		double gamma;
		double[] xk = InitialGuess;
		double[] gk = objectiveFunction.Gradient(InitialGuess);
		double[] pkLast = new double[gk.length];
		double[] gkLast = new double[gk.length];
		double[] xkLast = new double[gk.length];
		double[] pk = new double[gk.length];//gk;//new double[gk.length];
		double[] tmp = new double[gk.length];
		double[] sk = new double[gk.length];
		double[] skLast = new double[gk.length];
		double[] yk = new double[gk.length];
		double[] rk = new double[gk.length];
		double[] rkTmp = new double[gk.length^^2];
		double[] ykLast = new double[gk.length];
		double[] Bk = new double[gk.length^^2];
		double[] BkLast = new double[gk.length^^2];
		double[] Bktmp = new double[gk.length^^2];
		double[] lambdaK = new double[constraints];
		double[] A = new double[xk.length*constraints];
		double[] AT = new double[xk.length*constraints];
		double[] I = IdentityMatrix(cast(int)gk.length);
		double[] BA0 = new double[(xk.length+constraints)^^2];
		double[] gc = new double[xk.length+constraints];
		double[] pl = new double[xk.length+constraints];
		double[] c = new double[constraints];
		double[] skTmp = new double[sk.length];
		double[] skTmp1 = new double[sk.length];
		double µ = 0;
		double alpha = 1;
		double D = 0;
		double tau = 0.5;
		double theta = 1;
		double kkt1 = 1;
		double kkt2 = 1;
		//auto f = File(PointFilename, "w");
		//auto ferr = File(ErrorFilename, "w");
		File f;
		if(FileOutput) f = File(PointFilename, "w");
		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");
		Bk[] = I;

		c[] = ObjectiveFunc.Constraint(xk.complex()).Real();
		//if(DebugMode) writeln("Found initial constraint values. c = ", c);

		while( (abs(kkt1) > Tolerance) || (abs(kkt2) > Tolerance) )
		{
			//c[] = ObjectiveFunc.Constraint(xk.complex()).Real();

			if(FileOutput) WriteArrayCSV(f, xk);

			gc[0..xk.length] = -gk[];
			gc[xk.length..$] = -c[];
			//if(DebugMode) writeln("Spliced gc vector. gc = ", gc, "\tgk = ", gk);

			BA0 = MakeBAMatrix(Bk, Ak(xk), AkT(xk));
			//if(DebugMode) writeln("Made BA matrix");
			//if(DebugMode) PrintMatrix(BA0, cast(int)xk.length+constraints);
			MatrixView!double mat;
			
			mat.array = BA0;
			mat.rows = xk.length+constraints;
			
			pl = solve(mat, gc);


			pk[] = pl[0..xk.length];
			lambdaK[] = pl[xk.length..$];

			gamma = minPos!("a > b")(lambdaK)[0];

			if(!(1/µ >= delta + gamma))
			{
				µ = 1/(gamma + 2*delta);
			}

			alpha = 1;

			Ø.Mu = µ;

			tmp[] = xk[] + alpha*pk[];
			D = dot(cast(int)pk.length, cast(double*)Ø.Gradient(xk), 1, cast(double*)pk, 1);
			while( Ø.Compute(tmp.complex()).re > (Ø.Compute(xk.complex()).re + Eta*alpha*D) )
			{
				alpha = tau*alpha;
				tmp[] = xk[] + alpha*pk[];
				minorIterations++;
				//D = dot(cast(int)pk.length, cast(double*)Ø.Gradient(xk), 1, cast(double*)pk, 1);
			}

			xkLast[] = xk[];

			sk[] = alpha*pk[];
			xk[] = xkLast[] + alpha*pk[];
			yk = UpdateYk(xkLast, xk, lambdaK);

			double skykDot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)yk, 1);

			gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)Bk, cast(int)sk.length, 0, cast(double*)skTmp, 1);
			//gemm('n', 'n', 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)Bk, cast(int)sk.length, 0, cast(double*)skTmp, 1);
			double skBkdot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)skTmp, 1);

			if( skykDot >= 0.2*skBkdot)
			{
				theta = 1;
			}
			else if(skykDot < 0.2*skBkdot)
			{
				theta = (0.8*skBkdot)/(skBkdot - skykDot);
			}
			//CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans
			gemv(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)Bk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)skTmp, 1);
			//gemv('n', cast(int)sk.length, cast(int)sk.length, 1, cast(double*)Bk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)skTmp, 1);
			rk[] = theta*yk[] + (1 - theta)*skTmp[];

			BkLast[] = Bk[];

			gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)BkLast, cast(int)sk.length, 0, cast(double*)skTmp1, 1);
			
			gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)skTmp.length, cast(int)skTmp.length, 1, 1, cast(double*)skTmp, cast(int)skTmp.length, cast(double*)skTmp1, 1, 0, cast(double*)Bktmp, cast(int)skTmp.length);
			
			gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)rk.length, cast(int)rk.length, 1, 1, cast(double*)rk, cast(int)rk.length, cast(double*)rk, 1, 0, cast(double*)rkTmp, cast(int)rk.length);
			/*
			gemm('n', 'n', 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)BkLast, cast(int)sk.length, 0, cast(double*)skTmp1, 1);

			gemm('n', 'n', cast(int)skTmp.length, cast(int)skTmp.length, 1, 1, cast(double*)skTmp, cast(int)skTmp.length, cast(double*)skTmp1, 1, 0, cast(double*)Bktmp, cast(int)skTmp.length);

			gemm('n', 'n', cast(int)rk.length, cast(int)rk.length, 1, 1, cast(double*)rk, cast(int)rk.length, cast(double*)rk, 1, 0, cast(double*)rkTmp, cast(int)rk.length);
			*/
			double skrkDot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)rk, 1);


			//if((iterations%4) || (iterations == 0) )
			if(true)
			{
				Bk[] = BkLast[] - (1/skBkdot)*Bktmp[] + (1/skrkDot)*rkTmp[];
			}
			else
			{
				//if(DebugMode) writeln("Here");
				Bk[] = I;
			}

			kkt1 = Kkt1(xk, lambdaK);

			c[] = ObjectiveFunc.Constraint(xk.complex()).Real();
			gk[] = objectiveFunction.Gradient(xk);

			kkt2 = 0;
			for(int i = 0; i < constraints; i++)
			{
				kkt2 += c[i];
			}

			if(DebugMode)
			{
				writeln("gk = ", gk, "\t\t\tpk = ", pk);
				writeln("f = ", ObjectiveFunc.Compute(xk.complex()), "\tKKT1 = ", kkt1, "\t\t\t\tKKT2 = ", kkt2, "\t\t\t\txk = ", xk);
				writeln("------------------------------------------------------------------------------------------------------------");
				writeln();
				Thread.sleep(dur!("msecs")(150));
			}
			
			iterations++;

		}

		if(DebugMode) writeln("Done!");

		result.DesignVariables = xk[];
		result.Iterations = iterations;
		result.ObjectiveFunctionValue = ObjectiveFunc.Compute(xk.complex()).re;
		result.MinorIterations = minorIterations;

		return result;
	}

	private double Kkt1(double[] xk, double[] lambda)
	{
		double[] gk = ObjectiveFunc.Gradient(xk);
		double[][] cGrad = ObjectiveFunc.ConstraintGradient(xk);
		double[] cSum = new double[xk.length];
		double[] yk = new double[xk.length];
		double gSum = 0;
		cSum[] = 0;
		
		for(int i = 0; i < gk.length; i++)
		{
			for(int j = 0; j < cGrad.length; j++)
			{
				cSum[i] += lambda[j]*cGrad[j][i];
			}
		}

		for(int i = 0; i < gk.length; i++)
		{
			gSum += (gk[i] - cSum[i]);
		}

		return gSum;
	}

	private double[] UpdateYk(double[] xk, double[] xkNext, double[] lambda)
	{
		double[] gk = ObjectiveFunc.Gradient(xk);
		double[] gkNext = ObjectiveFunc.Gradient(xkNext);
		double[][] cGrad = ObjectiveFunc.ConstraintGradient(xk);
		double[][] cGradNext = ObjectiveFunc.ConstraintGradient(xkNext);
		double[] cSum = new double[xk.length];
		double[] cSumNext = new double[xk.length];
		double[] yk = new double[xk.length];

		cSum[] = 0;
		cSumNext[] = 0;

		for(int i = 0; i < gk.length; i++)
		{
			for(int j = 0; j < cGrad.length; j++)
			{
				cSum[i] += lambda[j]*cGrad[j][i];
				cSumNext[i] += lambda[j]*cGradNext[j][i];
			}
		}

		yk[] = (gkNext[] - cSumNext[]) - (gk[] - cSum[]);

		return yk;
	}

	private double[] MakeBAMatrix(double[] B, double[] A, double[] AT)
	{
		int constraints = ObjectiveFunc.Constraints;
		int dims = cast(int)A.length/constraints;
		int csize = dims+constraints;
		int size = csize^^2;
		double[] tmpBA = new double[size];
		tmpBA[] = 0;
		int j = 0;
		int j1 = 0;
		
		//writeln("Dims = ", dims, " csize = ", csize, " size = ", size);
		for(int i = 0; i < size; i+=csize)
		{
			if(j < dims)
			{
				//writeln("if... i = ", i, " j = ", j);
				tmpBA[i..i+(dims)] = B[j*dims..(j*dims)+(dims)];
				//writeln("if first");
				tmpBA[i+dims..i+dims+(constraints)] = A[j*constraints..(j*constraints)+(constraints)];
				//writeln("if second");
			}
			else
			{
				//writeln("else... i = ", i, " j1 = ", j1);
				tmpBA[i..i+(dims)] = -AT[j1*dims..(j1*dims)+(dims)];
				j1++;
			}
			j++;
		}
		//tmpBA
		return tmpBA;
	}

	private double[] Ak(double[] xk)
	{
		double[] A = new double[xk.length*ObjectiveFunc.Constraints];
		double[][] tmpA = ObjectiveFunc.ConstraintGradient(xk);
		int cl = ObjectiveFunc.Constraints;
		int i1 = 0;
		//if(DebugMode) writeln("Ak: About to compute. cl = ", cl, " xk.length = ", xk.length);
		for(int i = 0; i < xk.length; i+=cl, i1++)
		{
			//if(DebugMode) writeln("Ak: i = ", i, "\t tmpA[0..$][i] = ", tmpA[0..$][i]);
			//A[i*xk.length..(i*xk.length)+(cl-1)] = tmpA[0..$][i];
			for(int j = 0; j < cl; j++)
			{
				A[i+j] = tmpA[j][i];
			}
		}

		return A;
	}

	private double[] AkT(double[] xk)
	{
		double[] A = new double[xk.length*ObjectiveFunc.Constraints];
		double[][] tmpA = ObjectiveFunc.ConstraintGradient(xk);
		int cl = ObjectiveFunc.Constraints;
		for(int i = 0; i < cl; i++)
		{
			A[i*cl..(i*cl)+(xk.length)] = tmpA[i][0..$];
		}
		
		return A;
	}

	@property void ObjectiveFunc(ObjectiveFunction objectiveFunction) { mObjectiveFunction = objectiveFunction; }
	@property ObjectiveFunction ObjectiveFunc() { return mObjectiveFunction; }

	@property void Eta(double eta) { mEta = eta; }
	@property double Eta() { return mEta; }

	private ObjectiveFunction mObjectiveFunction;
	private double[] mLambdaInitial;
	private double mEta = 0.1;
}