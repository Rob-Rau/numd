module Optimization.Optimizer;

import std.complex;
import std.math;
import std.stdio;
import core.time;
import Optimization.ObjectiveFunction;

struct Result
{
	version(X86)
	{
		uint Iterations;
		uint MinorIterations;
	}
	else
	{
		ulong Iterations;
		ulong MinorIterations;
	}
	double ObjectiveFunctionValue;
	double[] DesignVariables;
	double[] IterArr;
	double[] Error;
	long ComputationTime;
}

abstract class Optimizer
{
	public Result Optimize(ObjectiveFunction objectiveFunction)
	{
		Result result;
		auto Clock = new TickDuration;
		auto startTime = Clock.currSystemTick;

		result = doOptimize(objectiveFunction);

		auto endTime = Clock.currSystemTick;
		auto elapsedTime = endTime - startTime;
		result.ComputationTime = elapsedTime.usecs;

		return result;
	}

	@property void InitialGuess(double[] initGuess) { mInitialGuess = initGuess; }
	@property double[] InitialGuess() { return mInitialGuess; }
	@property void InitialGuess(double initGuess)
	{
		mInitialGuess.length = 1;
		mInitialGuess[0] = initGuess;
	}

	@property void Interval(double interval) { mInterval = interval; }
	@property double Interval() { return mInterval; }

	@property void Tolerance(double tolerance) { mTolerance = tolerance; }
	@property double Tolerance() { return mTolerance; }

	@property void DebugMode(bool debugMode) { mDebugMode = debugMode; }
	@property bool DebugMode() { return mDebugMode; }

	@property void FileOutput(bool fileOutput) { mFileOutput = fileOutput; }
	@property bool FileOutput() { return mFileOutput; }

	@property void PointFilename(string filename) { mPointFilename = filename; }
	@property string PointFilename() { return mPointFilename; }

	@property void ErrorFilename(string filename) { mErrorFilename = filename; }
	@property string ErrorFilename() { return mErrorFilename; }

	protected Result doOptimize(ObjectiveFunction objectiveFunction);

	// Complex step size
	protected string mPointFilename;
	protected string mErrorFilename;

	protected const double h = 10e-10;
	protected double[] IntermediateSln;
	protected double[] IntermediateOFV;

	private bool mFileOutput = false;
	private double[] mInitialGuess = [0];
	private double mInterval = 1;
	private double mTolerance = 1.0e-6;
	private bool mDebugMode = false;
}