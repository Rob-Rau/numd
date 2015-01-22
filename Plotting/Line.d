module Plotting.Line;

import plplot;
import Plotting.FigPlot;
import Plotting.Color;

import std.algorithm;
import std.math;
import std.stdio;

class Line2D : Plot
{
public:
	this(){}

	void AddDataSet(double[] x, double[] y, double[] z = [])
	{
		mX.length++;
		mY.length++;
		mColors.length++;
		mX[$-1] = x;
		mY[$-1] = y;
		mColors[$-1] = cast(Color)Black;
		mColors[$-1].Red = 100;
	}

	void AddLineColor(Color color)
	{
	    //mColors.length++;
	    mColors[mColorPushIndex++] = color;
	}
	
	void AddLineColor(const Color color)
	{
	    AddLineColor(cast(Color)color);
	}

    void AddLineColor(Color color, int lineIndex)
	{
	    assert(mColors.length > lineIndex);
	    mColors[lineIndex] = color;
	}
	
	void AddLineColor(const Color color, int lineIndex)
	{
	    AddLineColor(cast(Color)color, lineIndex);
	}

	void LogX()
	{
		mAxis = 10;
	}

	void LogY()
	{
		mAxis = 20;
	}

	void LogLog()
	{
		mAxis = 30;
	}
	
	protected override void Draw()
	{
	    double xMin, xMax;
	    double yMin, yMax;
	    double[] xMinTmp, xMaxTmp;
	    double[] yMinTmp, yMaxTmp;

        xMinTmp = new double[mX.length];
        xMaxTmp = new double[mX.length];
        yMinTmp = new double[mY.length];
        yMaxTmp = new double[mY.length];

		double[][] tmp;
		double[][] tmp1;
		if(mAxis == 20)
		{
			tmp.length = mY.length;
			for(int i = 0; i < mY.length; i++)
			{
				tmp[i].length = mY[i].length;
				//writeln(i);
				for(int j = 0; j < mY[i].length; j++)
				{
					//writeln(mY[i][j]);
					if(mY[i][j] == 0)
					{
						mY[i][j] = mY[i][j-1];
					}
					tmp[i][j] = log(cast(real)mY[i][j]);
				}
			}
		}
		else if(mAxis == 30)
		{
			//writeln("LogLog");
			tmp.length = mY.length;
			tmp1.length = mX.length;
			for(int i = 0; i < mY.length; i++)
			{
				tmp[i].length = mY[i].length;
				tmp1[i].length = mX[i].length;
				//writeln(i);
				for(int j = 0; j < mY[i].length; j++)
				{
					//writeln(mY[i][j]);
					if(mY[i][j] == 0)
					{
						mY[i][j] = mY[i][j-1];
					}
					if(mX[i][j] == 0)
					{
						mX[i][j] = mX[i][j-1];
					}
					tmp[i][j] = log(cast(real)mY[i][j]);
					tmp1[i][j] = log(cast(real)mX[i][j]);
				}
			}
		}

		//writeln("past log");

        for(int i = 0; i < mX.length; i++)
        {
            xMinTmp[i] = minPos(mX[i])[0];
            xMaxTmp[i] = minPos!("a > b")(mX[i])[0];

			if(mAxis == 20)
			{
	            yMinTmp[i] = minPos(tmp[i])[0];
				yMaxTmp[i] = minPos!("a > b")(tmp[i])[0];
			}
			else if(mAxis == 30)
			{
				yMinTmp[i] = minPos(tmp[i])[0];
				yMaxTmp[i] = minPos!("a > b")(tmp[i])[0];
				xMinTmp[i] = minPos(tmp1[i])[0];
				xMaxTmp[i] = minPos!("a > b")(tmp1[i])[0];
			}
			else
			{
				yMinTmp[i] = minPos(mY[i])[0];
				yMaxTmp[i] = minPos!("a > b")(mY[i])[0];
			}
        }
        xMin = minPos(xMinTmp)[0];
        xMax = minPos!("a > b")(xMaxTmp)[0];
		yMin = minPos(yMinTmp)[0];
		yMax = minPos!("a > b")(yMaxTmp)[0];


		//writeln("yMin = ", yMin, "\tyMax = ", yMax);
		//writeln("xMin = ", xMin, "\txMax = ", xMax);

        plscol0(15, 0, 0, 0);
        plcol0(15);
		plenv(xMin, xMax, yMin, yMax, 0, mAxis );
		plcol0(15);
		plbox("al", 1.0e-1, 1, "v", 0.1, 0);
	    // Create a labelled box to hold the plot.
		plcol0(15);
		pllab(LabelX, LabelY, Title);


		for(int i = 0; i < mX.length; i++)
		{
		    plscol0(15, mColors[i].Red, mColors[i].Green, mColors[i].Blue);
		    plcol0(15);

			if(mAxis == 20)
			{
				plline(mX[i], tmp[i]);
			}
			else if(mAxis == 30)
			{
				plline(tmp1[i], tmp[i]);
			}
			else
			{
				plline(mX[i], mY[i]);
			}
		}
	}

private:
	int mColorPushIndex = 0;
	double[][] mX;
	double[][] mY;
	Color[] mColors;
	int mAxis = 0;
}
