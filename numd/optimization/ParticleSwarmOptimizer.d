module Optimization.ParticleSwarmOptimizer;

import Optimization.ArrayOps;
import Optimization.Complex;
import Optimization.Optimizer;
import Optimization.ObjectiveFunction;
import Optimization.BracketAndZoom;

import core.thread;

import std.algorithm;
import std.conv;
import std.math;
import std.stdio;

//import scid.bindings.lapack.dlapack;
//import scid.bindings.blas.dblas;	
import cblas;

class ParticleSwarmOptimizer : Optimizer
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
		Result result;
		double[][] position = new double[][](Particles, Dimensions);
		double[][] positionTmp = new double[][](Particles, Dimensions);
		double[][] velocity = new double[][](Particles, Dimensions);
		double[] F = new double[Particles];
		uint globalBestIndex;
		double globalBestValue;
		double globalBestValueLast;
		double[] globalBestPosition = new double[Dimensions];
		double[] globalBestPositionLast = new double[Dimensions];
		double[][] localBestPosition = new double[][](Particles, Dimensions);
		double[] localBestPosAve = new double[Dimensions];
		double[] localBestValue = new double[Particles];
		double localBestValueAve;
		double[] upperBound = new double[Dimensions];
		double[] lowerBound = new double[Dimensions];
		double[] tmp = new double[Dimensions];
		double[] errorTmp = new double[Particles];
		int zeroCount = 0;
		double error = 1;
		double errorLast = 1;
		double dt = 0.5;
		//double w = 0.729;
		double w = 0.8;
		//double c1 = 1.49445;
		double c1 = 0.7;
		//double c2 = 1.49445;
		double c2 = 0.9;
		File[] f = new File[Particles];
		string[] FileNames = new string[Particles];
		if(FileOutput)
		{
			for(int i = 0; i < Particles; i++)
			{
				FileNames[i] = PointFilename ~ "Particle" ~ to!string(i) ~ ".csv";
				f[i] = File(FileNames[i], "w");
			}
		}

		File ferr;
		if(FileOutput) ferr = File(ErrorFilename, "w");

		for(int i = 0; i < Dimensions; i++)
		{
			//if(DebugMode) writeln(i);
			lowerBound[i] = Bounds[i][0];
			upperBound[i] = Bounds[i][1];
		}

		for(int i = 0; i < Particles; i++)
		{
			//if(DebugMode) writeln(RandomVector(Dimensions, 0, 1)[]);
			position[i][] = (lowerBound[] - upperBound[])*RandomVector(Dimensions, 0, 1)[] + upperBound[];
			//if(DebugMode) writeln(position[i][]);
			localBestPosition[i][] = position[i][];
			F[i] = objectiveFunction.Compute(position[i].complex()).re;
			localBestValue[i] = F[i];
			velocity[i][] = RandomVector(Dimensions, -1.0, 1.0)[];
		}

		double[] minArr = minPos(F);
		globalBestIndex = cast(uint)(F.length - minArr.length);
		globalBestValue = F[globalBestIndex];
		globalBestPosition[] = position[globalBestIndex][];

		if(DebugMode) writeln("Starting main loop");
		while( (iterations < 10000) && (error > Tolerance) )
		//while(iterations < 10000)
		{
			if(FileOutput)
			{
				for(int i = 0; i < Particles; i++)
				{
					WriteArrayCSV(f[i], position[i][]);
				}
			}
			/*
			if(DebugMode)
			{
				writeln();
				for(int i = 0; i < Particles; i++)
					writeln(position[i][]);
				writeln();
			}
			*/
			for(int i = 0; i < Particles; i++)
			{
				velocity[i][] = w*velocity[i][] + (dt^^-1)*c1*RandomVector(Dimensions, 0, 1)[]*(localBestPosition[i][] - position[i][]) + (dt^^-1)*c2*RandomVector(Dimensions, 0, 1)[]*(globalBestPosition[] - position[i][]);
				positionTmp[i][] = position[i][] + dt*velocity[i][];

				for(int j = 0; j < Dimensions; j++)
				{
					if((positionTmp[i][j] > Bounds[j][1]) || (positionTmp[i][j] < Bounds[j][0]))
					{
						velocity[i][j] = -0.5*velocity[i][j];
						//if(DebugMode) writeln(positionTmp[i][j], " -> ", position[i][j] + dt*velocity[i][j]);
					}
				}

				position[i][] = position[i][] + dt*velocity[i][];

				F[i] = objectiveFunction.Compute(position[i].complex()).re;

				if(F[i] < localBestValue[i])
				{
					localBestValue[i] = F[i];
					localBestPosition[i][] = position[i][];
				}
			}

			minArr = minPos(F);
			auto globalBestIndexTmp = cast(uint)(F.length - minArr.length);

			globalBestPositionLast[] = globalBestPosition[];
			globalBestValueLast = globalBestValue;

			if(F[globalBestIndexTmp] < globalBestValue)
			{
				globalBestIndex = globalBestIndexTmp;
				globalBestValue = F[globalBestIndex];
				globalBestPosition[] = position[globalBestIndex][];
			}

			//for(int i = 0; i < Particles; i++)
			//{
			//tmp[] = globalBestPosition[] - position[globalBestIndexTmp][];
			//tmp[] = position[globalBestIndexTmp][] - globalBestPositionLast[];
			//tmp[] = globalBestPosition[] - globalBestPositionLast[];

			localBestPosAve[] = 0.0;

			for(int i = 0; i < Dimensions; i++)
			{
				for(int j = 0; j < Particles; j++)
				{
					localBestPosAve[i] += localBestPosition[j][i];
					//localBestPosAve[i] += position[j][i];
				}
				localBestPosAve[i] /= Particles;
			}
			localBestValueAve = 0;
			for(int i = 0; i < Particles; i++)
			{
				localBestValueAve += localBestValue[i];
				//localBestValueAve += F[i];
			}
			localBestValueAve /= Particles;

			//if(DebugMode) writefln("local best position average = [%20.20f, %20.20f]", localBestPosAve[0], localBestPosAve[1]);	
			tmp[] = globalBestPosition[] - localBestPosAve[];
			//tmp[] = localBestPosAve[] - globalBestPosition[];
			//error = nrm2(cast(int)tmp.length, cast(double*)(tmp), 1)/(1 + nrm2(cast(int)globalBestPosition.length, cast(double*)globalBestPosition, 1)) + abs(localBestValueAve - globalBestValue)/(1+abs(globalBestValue));
			error = nrm2(cast(int)tmp.length, cast(double*)(tmp), 1)/(1 + nrm2(cast(int)localBestPosAve.length, cast(double*)localBestPosAve, 1)) + abs(globalBestValue - localBestValueAve)/(1+abs(localBestValueAve));

			if(FileOutput) ferr.writefln("%d, %40.40f", iterations, error);
			//if(error < 2*Tolerance)
			//	w *= 0.9;

			if(error == errorLast)
			{
				zeroCount++;
				//writeln("zero");
			}

			if(zeroCount > 50)
			{
				zeroCount = 0;
				w = w*0.5;
				if(DebugMode) writeln("w = ", w);
			}

			errorLast = error;
			//error = nrm2(cast(int)tmp.length, cast(double*)(tmp), 1)/(1 + nrm2(cast(int)globalBestPositionLast.length, cast(double*)globalBestPositionLast, 1)) + abs(F[globalBestIndexTmp]-globalBestValueLast)/(1+abs(globalBestValueLast));
			//}
			//error = nrm2(cast(int)errorTmp.length, cast(double*)error, 1);
			//error = 
			//if(DebugMode) writefln("%20.20f", error);
			iterations++;
		}

		if(FileOutput)
		{
			for(int i = 0; i < Particles; i++)
			{
				f[i].close();
			}
		}
		ferr.close();

		result.Iterations = iterations;
		result.DesignVariables = globalBestPosition[];
		result.ObjectiveFunctionValue = globalBestValue;

		return result;
	}

	@property void Particles(int particles) { mParticles = particles; }
	@property int Particles() { return mParticles; }

	@property void Dimensions(int dimensions) { mDimensions = dimensions; }
	@property int Dimensions() { return mDimensions; }

	@property void Bounds(double[][] bounds) { mBounds = bounds; }
	@property double[][] Bounds() { return mBounds; }

	private double[][] mBounds;
	private int mDimensions;
	private int mParticles;
}