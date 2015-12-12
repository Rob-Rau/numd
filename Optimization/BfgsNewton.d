module Optimization.BfgsNewton;

import Optimization.ArrayOps;
import Optimization.BracketAndZoom;
import Optimization.MatrixOps;
import Optimization.ObjectiveFunction;
import Optimization.Optimizer;

import LinearAlgebra.Matrix;
//import LinearAlgebra.Vector;

import core.thread;

import std.math;
import std.stdio;

//import scid.bindings.lapack.dlapack;
//import scid.bindings.blas.dblas;
import cblas;

class BfgsNewton(size_t dims) : Optimizer
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
		auto lineSearch = new BracketAndZoom;
		lineSearch.DebugMode = DebugMode;
		Result lineResult;
		Result lineResultLast = lineResult;
		Result result;
		bool converged = false;
		double error = 1;
		double new_error = 1;

		alias Vector!(dims, double) Vec;
		alias Matrix!(dims, dims, double) Mat;

		double[] xk = InitialGuess;
		auto vec_xk = Vec(InitialGuess[]);
		double[] gk = objectiveFunction.Gradient(InitialGuess);
		auto vec_gk = Vec(objectiveFunction.Gradient(InitialGuess));
		double[] pkLast = new double[gk.length];
		auto vec_pkLast = Vec(0);
		double[] gkLast = new double[gk.length];
		auto vec_gkLast = Vec(0);
		double[] xkLast = new double[gk.length];
		auto vec_xkLast = Vec(0);
		double[] pk = new double[gk.length];//gk;//new double[gk.length];
		auto vec_pk = Vec(0);
		double[] tmp = new double[gk.length];
		auto vec_tmp = Vec(0);
		double[] sk = new double[gk.length];
		auto vec_sk = Vec(0);
		double[] skLast = new double[gk.length];
		auto vec_skLast = Vec(0);
		double[] yk = new double[gk.length];
		auto vec_yk = Vec(0);
		double[] ykLast = new double[gk.length];
		auto vec_ykLast = Vec(0);
		double[] Vk = new double[gk.length^^2];
		auto mat_Vk = Mat.Identity();//new Mat();
		double[] VkLast = new double[gk.length^^2];
		auto mat_VkLast = Mat.Identity();//new Mat();
		double[] I = IdentityMatrix(cast(int)gk.length);
		auto mat_I = Mat.Identity();

		//bool converged = false;
		File f;
		if(FileOutput) f = File(PointFilename, "w");
		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");
		double beta;

		double[] Vtmp1 = new double[Vk.length];
		double[] Vtmp2 = new double[Vk.length];
		double[] Vtmp3 = new double[Vk.length];
		auto mat_Vtmp1 = Mat(0);
		auto mat_Vtmp2 = Mat(0);
		auto mat_Vtmp3 = Mat(0);

		Vk[] = I;
		VkLast[] = I;

		if(FileOutput) WriteArrayCSV(f, xk);
		//lineSearch.P.length = pk.length;
		//pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];

		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		Vtmp1[] = -Vk[];
		mat_Vtmp1 = -mat_Vk;
		//gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk[], 1, 0, cast(double*)pk, 1);
		//gemv(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk[], 1, 0, cast(double*)pk, 1);
		gemv(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk[], 1, 0, cast(double*)pk, 1);
		pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];
		vec_pk = (mat_Vtmp1*vec_gk).normalize();
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);
		/*
		writeln("gk = ", gk, "\tpk = ", pk);
		writeln("vec_gk = ", vec_gk.ToString());
		writeln("vec_pk = ", vec_pk.ToString());
		*/
		//lineSearch.P = pk;
		//lineSearch.InitialGuess = xk;
		lineSearch.P = vec_pk.mData;
		lineSearch.InitialGuess = vec_xk.mData;
		lineResult = lineSearch.Optimize(objectiveFunction);
		minorIterations += lineResult.Iterations;
		lineResultLast = lineResult;

		//vec_xkLast.swap(vec_xk);
		//vec_pkLast.swap(vec_pk);
		//vec_gkLast.swap(vec_gk);
		vec_xkLast = vec_xk;
		vec_pkLast = vec_pk;
		vec_gkLast = vec_gk;
		xkLast[] = xk;
		pkLast[] = pk;
		gkLast[] = gk;
		xk = lineResult.DesignVariables;
		vec_xk[] = lineResult.DesignVariables;
		if(DebugMode) writeln("xk = ", xk);
		
		// Compute gradient at new point.
		gk[] = objectiveFunction.Gradient(xk);
		vec_gk[] = objectiveFunction.Gradient(vec_xk[]);
		//vec_gk = gk;
		// Compute new direction at point.
		//pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);


		vec_yk = vec_gk - vec_gkLast;
		vec_sk = vec_xk - vec_xkLast;
		yk[] = gk[] - gkLast[];
		sk[] = xk[] - xkLast[];

		double skykDot = dot(cast(int)yk.length, cast(double*)yk, 1, cast(double*)sk, 1);
		auto vec_skykDot = vec_yk.dot(vec_sk);
		double[] skyk = new double[Vk.length];
		double[] yksk = new double[Vk.length];
		double[] sksk = new double[Vk.length];

		//CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans
		//gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
		//gemm(Order.RowMajor, Transpose.NoTrans, Transpose.NoTrans, cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
		gemm(Order.RowMajor, Transpose.NoTrans, Transpose.NoTrans, cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, 1, cast(double*)yk, cast(int)yk.length, 0, cast(double*)skyk, cast(int)sk.length);
		auto mat_skyk = vec_sk*vec_yk.transpose();
		if(DebugMode) writeln("there");
		gemm(Order.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)1, cast(double*)sk, cast(int)sk.length, 0, cast(double*)yksk, cast(int)yk.length);
		//gemm(Order.ColMajor, Transpose.NoTrans, Transpose.NoTrans, cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);
		auto mat_yksk = vec_yk*vec_sk.transpose();
		//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
		//gemm('n', 'n', cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);
		/*
		auto tmp_skyk = Mat(skyk);
		auto tmp_yksk = Mat(yksk);

		writeln(sk);
		writeln(yk);

		writeln("skyk");
		writeln(tmp_skyk.ToString());
		writeln(mat_skyk.ToString());

		writeln("yksk");
		writeln(tmp_yksk.ToString());
		writeln(mat_yksk.ToString());
		*/
		if(DebugMode) writeln("here");

		skyk[] *= skykDot^^(-1);
		mat_skyk *= vec_skykDot^^(-1);
		Vtmp1[] = I[] - skyk[];
		mat_Vtmp1 = mat_I - mat_skyk;
		gemm(Order.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)Vk, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
		mat_Vtmp2 = mat_Vtmp1*mat_Vk;
		gemm(Order.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
		auto mat_sksk = vec_sk*vec_sk.transpose();
		//gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)Vk, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
		//gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
		//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)Vk, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
		//gemm('n', 'n', cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
		sksk[] *= skykDot^^(-1);
		mat_sksk *= vec_skykDot^^(-1);
		/*
		auto tmp_mat = Mat(sksk);
		writeln("sksk");
		writeln(tmp_mat.ToString());
		writeln(mat_sksk.ToString());
		*/

		mat_Vtmp1 = mat_I - mat_yksk/vec_skykDot;
		Vtmp1[] = I[] - yksk[]*skykDot^^(-1);

		//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);
		gemm(Order.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);
		//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);
		mat_Vtmp3 = mat_Vtmp2*mat_Vtmp1;

		/*
		auto tmp_Vk = Mat(Vk);
		auto tmp_Vtmp1 = Mat(Vtmp1);
		auto tmp_Vtmp2 = Mat(Vtmp2);
		auto tmp_Vtmp3 = Mat(Vtmp3);

		writeln("Vk");
		writeln(tmp_Vk.ToString());
		writeln(mat_Vk.ToString());
		
		writeln("Vtmp1");
		writeln(tmp_Vtmp1.ToString());
		writeln(mat_Vtmp1.ToString());
		
		writeln("Vtmp2");
		writeln(tmp_Vtmp2.ToString());
		writeln(mat_Vtmp2.ToString());
		
		writeln("Vtmp3");
		writeln(tmp_Vtmp3.ToString());
		writeln(mat_Vtmp3.ToString());
		*/

		mat_Vk = mat_Vtmp3 + mat_sksk;
		Vk[] = Vtmp3[] + sksk[];
		VkLast[] = Vk;
		//mat_VkLast.swap(mat_Vk);
		mat_VkLast = mat_Vk;

		Vtmp1[] = -Vk[];
		mat_Vtmp1 = -mat_Vk;
		//gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
		gemv(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
		//gemv(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);

		vec_pk = (mat_Vtmp1*vec_gk).normalize();
		//vec_pk = vec_pk.normalize();
		pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];
		if(DebugMode) writeln("gk = ", gk, "\tpk = ", pk);

		//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);
		lineSearch.AlphaInitial = 1;//(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
		
		double epa = 1.0e-6;
		double epr = 1.0e-6;
		double epg = 1.0e-3;
		/*
		auto tmp_xk = Vec(xk);
		auto tmp_gk = Vec(gk);
		auto tmp_pk = Vec(pk);

		writeln("xk");
		writeln(tmp_xk.ToString());
		writeln(vec_xk.ToString());

		writeln("gk");
		writeln(tmp_gk.ToString());
		writeln(vec_gk.ToString());

		writeln("pk");
		writeln(tmp_pk.ToString());
		writeln(vec_pk.ToString());
		*/
		//while(error >= Tolerance)
		//while( !(abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) > (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && !(nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1) <= epg) )
		while(!converged)
		{

			if(DebugMode) writeln();
			if(FileOutput) WriteArrayCSV(f, xk);

			//gemm('n', 'n', 2, 2, 1, 1, cast(double*)v1, 2, cast(double*)v2, 1, 0, cast(double*)ans, 2);

			//writeln("Alpha = ", lineSearch.AlphaInitial);
			//lineSearch.AlphaInitial = 1;
			//lineSearch.AlphaInitial = lineSearch.AlphaInitial*(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//writeln("gk-1*pk-1 = ", dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1));
			//writeln("gk*pk = ", dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//writeln("AlphaNew = ", lineSearch.AlphaInitial);
			//lineSearch.P = pk;
			if(DebugMode) writeln("Assign pk");
			lineSearch.P = vec_pk.mData;

			if(DebugMode) writeln("Assign xk");
			//lineSearch.InitialGuess[] = xk;
			lineSearch.InitialGuess = vec_xk.mData;
			if(DebugMode) writeln("opt");
			lineResult = lineSearch.Optimize(objectiveFunction);
			if(DebugMode) writeln("line opt");
			minorIterations += lineResult.Iterations;
			//lineResultLast = lineResult;
			vec_xkLast = vec_xk;
			xkLast[] = xk;

			vec_pkLast = vec_pk;
			pkLast[] = pk;

			vec_gkLast = vec_gk;
			gkLast[] = gk;
			if(DebugMode) writeln("getting lineResult");
			xk[] = lineResult.DesignVariables;
			vec_xk[] = lineResult.DesignVariables;
			// Compute gradient at new point.
			gk = objectiveFunction.Gradient(xk);
			//vec_gk = gk;
			vec_gk[] = objectiveFunction.Gradient(vec_xk);

			if(iterations%5)
			//if(true)
			//if(iterations%xk.length)
			{
				//beta = dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)gk, 1)/dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)gkLast, 1);

				yk[] = gk[] - gkLast[];
				vec_yk = vec_gk - vec_gkLast;
				sk[] = xk[] - xkLast[];
				vec_sk = vec_xk - vec_xkLast;

				skykDot = dot(cast(int)yk.length, cast(double*)yk, 1, cast(double*)sk, 1);
				vec_skykDot = vec_sk.dot(vec_yk);

				gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)1, cast(double*)yk, cast(int)yk.length, 0, cast(double*)skyk, cast(int)sk.length);
				mat_skyk = vec_sk*vec_yk.transpose();
				gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)1, cast(double*)sk, cast(int)sk.length, 0, cast(double*)yksk, cast(int)yk.length);
				mat_yksk = vec_yk*vec_sk.transpose();
				//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
				//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);
				//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)yk, 1, 0, cast(double*)skyk, cast(int)sk.length);
				//gemm('n', 'n', cast(int)yk.length, cast(int)sk.length, 1, 1, cast(double*)yk, cast(int)yk.length, cast(double*)sk, 1, 0, cast(double*)yksk, cast(int)yk.length);
				
				skyk[] *= skykDot^^(-1);
				mat_skyk /= vec_skykDot;
				Vtmp1[] = I[] - skyk[];
				mat_Vtmp1 = mat_I - mat_skyk;
				//writeln("here");
				gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)VkLast, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
				mat_Vtmp2 = mat_Vtmp1*mat_VkLast;
				//writeln("there");
				gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, 1, cast(double*)sk, cast(int)sk.length, 0, cast(double*)sksk, cast(int)sk.length);
				mat_sksk = vec_sk*vec_sk.transpose();
				//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)VkLast, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
				//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
				//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp1, cast(int)sk.length, cast(double*)VkLast, cast(int)sk.length, 0, cast(double*)Vtmp2, cast(int)sk.length);
				//gemm('n', 'n', cast(int)sk.length, cast(int)sk.length, 1, 1, cast(double*)sk, cast(int)sk.length, cast(double*)sk, 1, 0, cast(double*)sksk, cast(int)sk.length);
				sksk[] *= skykDot^^(-1);
				mat_sksk /= vec_skykDot;

				Vtmp1[] = I[] - yksk[]*skykDot^^(-1);
				mat_Vtmp1 = mat_I - mat_yksk/vec_skykDot;

				//gemm('n', 'n', cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);
				gemm(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);				
				mat_Vtmp3 = mat_Vtmp2*mat_Vtmp1;
				//gemm(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, CBLAS_TRANSPOSE.NoTrans, cast(int)sk.length, cast(int)yk.length, cast(int)yk.length, 1, cast(double*)Vtmp2, cast(int)sk.length, cast(double*)Vtmp1, cast(int)sk.length, 0, cast(double*)Vtmp3, cast(int)sk.length);				
				Vk[] = Vtmp3[] + sksk[];
				mat_Vk = mat_Vtmp3 + mat_sksk;
				VkLast[] = Vk;
				mat_VkLast = mat_Vk;
			
			}
			else
			{
				Vk[] = I;
				mat_Vk = Mat.Identity();
			}
			Vtmp1[] = -Vk[];
			mat_Vtmp1 = -mat_Vk;
			//CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans
			gemv(CBLAS_ORDER.RowMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
			vec_pk = (mat_Vtmp1*vec_gk).normalize();
			//gemv(CBLAS_ORDER.ColMajor, CBLAS_TRANSPOSE.NoTrans, cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
			//gemv('n', cast(int)gk.length, cast(int)gk.length, 1, cast(double*)Vtmp1, cast(int)gk.length, cast(double*)gk, 1, 0, cast(double*)pk, 1);
			pk[] = nrm2(cast(int)pk.length, cast(double*)pk, 1)^^(-1)*pk[];

			/*
			// Compute new direction at point.
			//pk[] = gk;
			//if(iterations%1.0)
			pk[] = -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1)*gk[] + beta*pkLast[];
			*/
			//scal(cast(int)pk.length, -nrm2(cast(int)gk.length, cast(double*)gk, 1)^^(-1), cast(double*)pk, 1);

			if(ErrorConvergence)
			{
				if(error <= Tolerance)
					converged = true;
			}
			else
			{
				//writefln(" current = %20.20f \t previous = %20.20f", lineResult.ObjectiveFunctionValue, lineResultLast.ObjectiveFunctionValue);
				//writefln("%20.20f should be > %20.20f", abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue), (epa + epr*abs(lineResultLast.ObjectiveFunctionValue)));
				//writefln("%20.20f should be < %20.20f", nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1), epg);
				//if((abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) < (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && (nrm2(cast(int)gkLast.length, cast(double*)gkLast, 1) <= epg))
				if((abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue) < (epa + epr*abs(lineResultLast.ObjectiveFunctionValue))) && (vec_gkLast.magnitude() <= epg))
					converged = true;
			}

			tmp[] = xk[] - xkLast[];
			vec_tmp = vec_xk - vec_xkLast;
			error = nrm2(cast(int)tmp.length, cast(double*)(tmp), 1)/(1 + nrm2(cast(int)xkLast.length, cast(double*)xkLast, 1)) + abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue)/(1+abs(lineResultLast.ObjectiveFunctionValue));
			error = vec_tmp.magnitude()/(1 + vec_xkLast.magnitude()) + abs(lineResult.ObjectiveFunctionValue - lineResultLast.ObjectiveFunctionValue)/(1+abs(lineResultLast.ObjectiveFunctionValue));
			lineResultLast = lineResult;
			lineSearch.AlphaInitial = 1;
			//lineSearch.AlphaInitial = lineSearch.AlphaInitial*(dot(cast(int)gkLast.length, cast(double*)gkLast, 1, cast(double*)pkLast, 1)/dot(cast(int)gk.length, cast(double*)gk, 1, cast(double*)pk, 1));
			//lineSearch.AlphaMax = lineSearch.AlphaInitial+10;
			if(DebugMode)
			{
				writeln("gkLast = ", gkLast, "\tpkLast = ", pkLast);
				writeln("gk = ", gk, "\tpk = ", pk);
				writeln("Error = ", error, "\txk = ", lineResult.DesignVariables);
				Thread.sleep(dur!("msecs")(50));
			}
			
			iterations++;
			if(FileOutput) ferr.writefln("%d, %40.40f", iterations, error);

		}
		ferr.close();
		f.close();
		result = lineResult;
		result.Iterations = iterations;
		result.MinorIterations = minorIterations;
		/*
		writeln("error");
		writeln(error);
		writeln(new_error);

		auto tmp_xk = Vec(xk);
		auto tmp_gk = Vec(gk);
		auto tmp_pk = Vec(pk);
		
		writeln("xk");
		writeln(tmp_xk.ToString());
		writeln(vec_xk.ToString());
		
		writeln("gk");
		writeln(tmp_gk.ToString());
		writeln(vec_gk.ToString());
		
		writeln("pk");
		writeln(tmp_pk.ToString());
		writeln(vec_pk.ToString());
		*/
		//writeln("Error = ", error, "\txk = ", result.DesignVariables);
		return result;
	}

	@property void ErrorConvergence(bool errConv) { mErrorConverg = errConv; }
	@property bool ErrorConvergence() { return mErrorConverg; }

	private bool mErrorConverg = false;
}

version(unittest)
{
	import core.memory;
	import Optimization.SecantRoot;
	import std.math;
	import std.complex;

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

	unittest
	{
		auto bfgs = new BfgsNewton!2;
		
		double[2] initialPoint = [10, 40];
		auto drag = new Drag;
		Result result;
		
		//bfgs.FileOutput = true;
		
		writeln();
		writeln("Starting drag optimization. Start point [35, 10]");
		bfgs.InitialGuess = [35, 10];
		
		writeln();
		//bfgs.DebugMode = true;
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

		assert(abs(result.ObjectiveFunctionValue - 191.9016258727) < 1.0e-6, "BFGS Quasi-Newton test 2 failed");
	}

	unittest
	{
		auto bfgs = new BfgsNewton!8;
		
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

		assert(abs(result.ObjectiveFunctionValue) < 1.0e-6, "BFGS Quasi-Newton test 1 failed");
	}

	unittest
	{
		auto bfgs = new BfgsNewton!4;
		
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
		assert(true);
	}
}