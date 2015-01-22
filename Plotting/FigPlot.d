module Plotting.FigPlot;

import Plotting.Color;
import plplot;

import std.process;

abstract class Plot
{
public:
	this(){}

	final void Label(string xLabel, string yLabel)
	{
		mLabelX = xLabel;
		mLabelY = yLabel;
	}

	@property void LabelX(string label) { mLabelX = label; }
	@property string LabelX() { return mLabelX; }

	@property void LabelY(string label) { mLabelY = label; }
	@property string LabelY() { return mLabelY; }

	@property void Title(string title) { mTitle = title; }
	@property string Title() { return mTitle; }

	protected void Draw();

private:
	string mLabelX;
	string mLabelY;
	string mTitle;
}

class Figure
{
public:
	this(){}

	void AddPlot(Plot plot)
	{
		mPlots.length++;
		mPlots[$-1] = plot;
	}

	void AddPlot(Plot plot, int rowPos, int colPos)
	{
        int index = Columns*rowPos + colPos;
        if(index >= mPlots.length)
        {
            mPlots.length = mPlots.length + ((index + 1) - mPlots.length);
        }
        mPlots[index] = plot;
	}

	void Layout(int rows, int columns)
	{
        Rows = rows;
        Columns = columns;
	}

	void Show()
	{
		version(Win32)
		{
			char[][] cArgs = [cast(char[])"./plplotApp", cast(char[])"-o", cast(char[])Filename];
		}
		version(linux)
		{
			char[][] cArgs;
			//char[][] cArgs = [cast(char[])"./plplotApp"];
			if(Driver != "xwin")
			{
				cArgs = [cast(char[])"./plplotApp", cast(char[])"-o", cast(char[])Filename];
			}
			else
			{
				cArgs = [cast(char[])"./plplotApp"];
			}
		}
		version(OSX)
		{
			char[][] cArgs;
			if(Driver != "xwin")
			{
				cArgs = [cast(char[])"./plplotApp", cast(char[])"-o", cast(char[])Filename];
			}
			else
			{
				cArgs = [cast(char[])"./plplotApp"];
			}
		}

		plparseopts( cArgs, PL_PARSE_FULL );

	    plscolbg(BackgroundColor.Red, BackgroundColor.Green, BackgroundColor.Blue);

        plstart(Driver, Columns, Rows);

		for(int i = 0; i < mPlots.length; i++)
		{
			mPlots[i].Draw();
		}

		plend();
	}

	@property void BackgroundColor(Color bgColor) { mBackgroundColor = bgColor; }
	@property Color BackgroundColor() { return mBackgroundColor; }

    @property void Driver(string driver) { mDriver = driver; }
    @property string Driver() { return mDriver; }

	@property void Filename(string filename) { mFilename = filename; }
	@property string Filename() { return mFilename; }

    @property void Rows(int rows) { mRows = rows; }
    @property int Rows() { return mRows; }

    @property void Columns(int cols) { mColumns = cols; }
    @property int Columns() { return mColumns; }

private:
	Plot[] mPlots;
	Color mBackgroundColor = new Color(255, 255, 255);
	version(linux)
	{
		string mDriver = "xwin";
	}
	version(OSX)
	{
		string mDriver = "xwin";
	}
	version(Win32)
	{
		string mDriver = "svg";
	}
	string mFilename = "output.svg";
	int mRows = 1;
	int mColumns = 1;
}
