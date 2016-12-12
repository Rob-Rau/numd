module numd.optimization.SQP;

import numd.optimization.ArrayOps;
import numd.optimization.BracketAndZoom;
import numd.optimization.Complex;
import numd.optimization.InteriorPenalty;
import numd.optimization.MatrixOps;
import numd.optimization.ObjectiveFunction;
import numd.optimization.Optimizer;

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
		//return ObjectiveFunc.Compute(designVar) + (1.0/Mu)*cNorm;
		return (1.0/Mu)*cNorm;
	}
	
	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];
		c[0] = 0;
		return c;
	}

	override void UpdateActiveSet(double[] x)
	{

	}
}

class SQP : Optimizer
{
	int id;
	bool stop;
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
		MeritFunction meritFunc = new MeritFunction;
		meritFunc.ObjectiveFunc = objectiveFunction;
		int constraints = ObjectiveFunc.Constraints;
		bool converged = false;
		double error = 1;
		double delta = 0.01;
		double gamma;
		double[] xk = new double[InitialGuess.length];
		xk[] = InitialGuess[];
		double[] gk = new double[xk.length];
		gk[] = objectiveFunction.Gradient(InitialGuess)[];
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

		double[] I = IdentityMatrix(cast(int)gk.length);

		double[] lambdaK = new double[constraints];
		double[] A = new double[xk.length*constraints];
		double[] AT = new double[xk.length*constraints];
		double[] BA0 = new double[(xk.length+constraints)^^2];
		double[] gc = new double[xk.length+constraints];
		double[] pl = new double[xk.length+constraints];
		double[] c = new double[constraints];
		
		double[] skTmp = new double[sk.length];
		double[] skTmp1 = new double[sk.length];
		double[][] cGrad = ObjectiveFunc.ConstraintGradient(xk);
		double[][] cGradLast = ObjectiveFunc.ConstraintGradient(xk);

		//double mu = 1.0e-8;
		double mu = 1.0;
		//double mu = 0;
		double alpha = 1;
		double D = 0;
		double tau = 0.5;
		double theta = 1;
		double kkt1 = 1;
		double kkt2 = 1;
		//auto f = File(PointFilename, "w");
		//auto ferr = File(ErrorFilename, "w");
		File f;
		if(FileOutput && (id == 0)) f = File(PointFilename, "a");
		File ferr;
		if(FileOutput && (id == 0)) ferr = File(ErrorFilename, "a");
		Bk[] = I;

		c[] = ObjectiveFunc.Constraint(xk.complex()).Real();
		//if(DebugMode) writeln("Found initial constraint values. c = ", c);

		double funcVal;
		double funcValLast;

		while( (((abs(kkt1) > Tolerance) || (abs(kkt2) > Tolerance)) && !stop) || (iterations <= 1))
		{
			//c[] = ObjectiveFunc.Constraint(xk.complex()).Real();

			if(FileOutput && (id == 0)) WriteArrayCSV(f, xk);

			gc[0..xk.length] = -gk[];
			gc[xk.length..$] = -c[];
			//if(DebugMode) writeln("Spliced gc vector. gc = ", gc, "\tgk = ", gk);

			BA0 = MakeBAMatrix(Bk, Ak(cGrad), AkT(cGrad));

			//if(DebugMode) writeln("Made BA matrix");
			//if(DebugMode) PrintMatrix(BA0, cast(int)xk.length+constraints);
			MatrixView!double mat;
			
			mat.array = BA0;
			mat.rows = xk.length+constraints;
			
			pl = solve(mat, gc);


			pkLast[] = pk[];
			pk[] = pl[0..xk.length];
			lambdaK[] = pl[xk.length..$];

			gamma = minPos!("a > b")(lambdaK)[0];

			if(!(1.0/mu >= (delta + gamma)))
			{
				mu = 1.0/(gamma + 2.0*delta);
			}

			if(id == 0) writeln("mu = ", mu);
			
			alpha = 1.0;
			//alpha = 0.05;

			//if(id == 0) writeln("mu = ", mu);
			meritFunc.Mu = mu;

			//tmp[] = xk[] + alpha*pk[];
			auto mg = meritFunc.Gradient(xk);
			tmp[] = gk[] + mg[];
			//D = dot(cast(int)pk.length, cast(double*)tmp, 1, cast(double*)pk, 1);
			tmp[] *= pk[];
			import std.algorithm : sum;
			D = tmp.sum;
			//D = dot(cast(int)pk.length, cast(double*)meritFunc.Gradient(xk), 1, cast(double*)pk, 1);
			funcValLast = funcVal;
			funcVal = objectiveFunction.Compute(xk.complex).re;
			auto mfxk = funcVal + meritFunc.Compute(xk.complex()).re;
			minorIterations = 0;
			do
			{
				tmp[] = xk[] + alpha*pk[];
				objectiveFunction.UpdateActiveSet(tmp);
				alpha = tau*alpha;
				//auto mfTmp = meritFunc.Compute(xk.complex()).re;
				if(id == 0) writeln("minor iteration = ", minorIterations);
				//if(id == 0) writeln("minor iteration = ", minorIterations, "; mfTmp = ", mfTmp, "; Eta*alpha*D = ", Eta*alpha*D, "; mfTmp + Eta*alpha*D = ", mfTmp + Eta*alpha*D);
				minorIterations++;

				if(minorIterations > 15)
				{
					if(id == 0) writeln("minorIterations > 15, breaking out");
					stop = true;
					break;
				}

				//D = dot(cast(int)pk.length, cast(double*)meritFunc.Gradient(xk), 1, cast(double*)pk, 1);
			} while( ((objectiveFunction.Compute(tmp.complex).re + meritFunc.Compute(tmp.complex()).re) > (mfxk + Eta*alpha*D)) && !stop );

			if(stop)
			{
				break;
			}

			alpha /= tau;
			
			xkLast[] = xk[];

			sk[] = alpha*pk[];
			xk[] = xkLast[] + alpha*pk[];

			cGradLast[] = cGrad[];
			cGrad = ObjectiveFunc.ConstraintGradient(xk);
			c[] = ObjectiveFunc.Constraint(xk.complex()).Real()[];
			gkLast[] = gk[];
			gk[] = objectiveFunction.Gradient(xk)[];

			objectiveFunction.UpdateActiveSet(xk);
			objectiveFunction.Constraint(xk.complex);
			constraints = objectiveFunction.Constraints;
			lambdaK.length = constraints;
			A.length = xk.length*constraints;
			AT.length = xk.length*constraints;
			BA0.length = (xk.length+constraints)^^2;
			gc.length = xk.length+constraints;
			pl.length = xk.length+constraints;
			c.length = constraints;
			//if(id == 0) writeln("minor iterations done; alpha = ", alpha, " xk = ", xk);

			yk = UpdateYk(gkLast, gk, cGradLast, cGrad, lambdaK);

			//double skykDot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)yk, 1);
			tmp[] = sk[]*yk[];
			double skykDot = tmp.sum;

			gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)Bk.ptr, cast(int)sk.length, 0, cast(double*)skTmp, 1);
			
			//double skBkdot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)skTmp, 1);
			tmp[] = sk[]*skTmp[];
			double skBkdot = tmp.sum;

			if( skykDot >= 0.2*skBkdot)
			{
				theta = 1;
			}
			else if(skykDot < 0.2*skBkdot)
			{
				theta = (0.8*skBkdot)/(skBkdot - skykDot);
			}

			gemv(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)Bk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)skTmp, 1);

			rk[] = theta*yk[] + (1 - theta)*skTmp[];

			BkLast[] = Bk[];

			gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, 1, cast(int)sk.length, cast(int)sk.length, 1, cast(double*)sk, 1, cast(double*)BkLast, cast(int)sk.length, 0, cast(double*)skTmp1, 1);

			gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)skTmp.length, cast(int)skTmp.length, 1, 1, cast(double*)skTmp, cast(int)skTmp.length, cast(double*)skTmp1, 1, 0, cast(double*)Bktmp, cast(int)skTmp.length);

			gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)rk.length, cast(int)rk.length, 1, 1, cast(double*)rk, cast(int)rk.length, cast(double*)rk, 1, 0, cast(double*)rkTmp, cast(int)rk.length);

			double skrkDot = dot(cast(int)sk.length, cast(double*)sk, 1, cast(double*)rk, 1);


			//if((iterations%4) || (iterations == 0) )
			if(false)
			{
				//if(DebugMode) writeln("Here");
				Bk[] = I;
			}
			else //if(true)
			{
				Bk[] = BkLast[] - (1/skBkdot)*Bktmp[] + (1/skrkDot)*rkTmp[];
			}
			

			kkt1 = Kkt1(gk, cGrad, lambdaK);

			kkt2 = 0;
			for(int i = 0; i < constraints; i++)
			{
				kkt2 += c[i];
			}

			if(DebugMode && (id == 0))
			{
				writeln("gk = ", gkLast, "\t\t\tpk = ", pkLast);
				writeln("f = ", funcVal, "\tKKT1 = ", kkt1, "\t\t\t\tKKT2 = ", kkt2, "\t\t\t\txk = ", xkLast, "\t\t\t\t%error = ", (funcVal - funcValLast)/funcVal);
				writeln("------------------------------------------------------------------------------------------------------------");
				writeln();
				Thread.sleep(dur!("msecs")(150));
			}
			
			iterations++;

		}

		if(DebugMode) writeln("Done!");

		result.DesignVariables.length = xk.length;
		result.DesignVariables[] = xk[];
		result.Iterations = iterations;
		result.ObjectiveFunctionValue = ObjectiveFunc.Compute(xk.complex()).re;
		result.MinorIterations = minorIterations;

		return result;
	}

	private double Kkt1(double[] gk, double[][] cGrad, double[] lambda)
	{
		double[] cSum = new double[gk.length];
		double[] yk = new double[gk.length];
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

	private double[] UpdateYk(double[] gk, double[] gkNext, double[][] cGrad, double[][] cGradNext, double[] lambda)
	{
		double[] cSum = new double[gk.length];
		double[] cSumNext = new double[gk.length];
		double[] yk = new double[gk.length];

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

	private double[] Ak(double[][] cGrad)
	{
		double[] A = new double[cGrad[0].length*ObjectiveFunc.Constraints];
		int cl = ObjectiveFunc.Constraints;
		int i1 = 0;
		for(int i = 0; i < cGrad[0].length; i+=cl, i1++)
		{
			for(int j = 0; j < cl; j++)
			{
				A[i+j] = cGrad[j][i];
			}
		}

		return A;
	}

	private double[] AkT(double[][] cGrad)
	{
		double[] A = new double[cGrad[0].length*ObjectiveFunc.Constraints];
		int cl = ObjectiveFunc.Constraints;
		for(int i = 0; i < cl; i++)
		{
			A[i*cl..(i*cl)+(cGrad[0].length)] = cGrad[i][0..$];
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

version(unittest)
{
	class Drag : ObjectiveFunction
	{
		import numd.optimization.SecantRoot;

		this() { Constraints = 1; }

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

		final override void UpdateActiveSet(double[] designVars)
		{

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

	unittest
	{
		auto sqp = new SQP;
		auto drag = new Drag;
		Result result;
		//double[2] initialGuess = [48, 52];
		double[2] initialGuess = [30, 30];

		drag.DerivativeType = "numd.optimization.FiniteDifference.FiniteDifference";
		drag.StepSize = 1.0e-3;
		//sqp.DebugMode = true;
		//sqp.InitialGuess = [8, 5];
		sqp.InitialGuess = new double[initialGuess.length];
		sqp.InitialGuess[] = initialGuess;
		sqp.PointFilename = "SQPpoints.csv";
		sqp.ErrorFilename = "SQPerror.csv";
		sqp.FileOutput = false;
		sqp.Eta = 0.1;
		writeln("Starting optimization");
		result = sqp.Optimize(drag);
		writeln();
		writeln("SQP:");
		writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
		writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
		writeln("\tConverged in ", result.Iterations, " iterations.");
		writeln("\tComputation time: ", result.ComputationTime, " usecs.");
		writeln("\tMinor iterations: ", result.MinorIterations);
	}
}