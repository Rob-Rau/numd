module Plotting.Color;

import plplot;

class Color
{
public:
    this(){}
    this(int r, int g, int b)
    {
        mRGB[0] = r;
        mRGB[1] = g;
        mRGB[2] = b;
    }

    this(int[3] rgb)
    {
        mRGB = rgb;
    }

    @property void Red(int r) { mRGB[0] = r; }
    @property int Red() { return mRGB[0]; }

    @property void Green(int g) { mRGB[1] = g; }
    @property int Green() { return mRGB[1]; }

    @property void Blue(int b) { mRGB[2] = b; }
    @property int Blue() { return mRGB[2]; }

private:
    int[3] mRGB;
}

export immutable Color White = cast(immutable)new Color(255, 255, 255);
const Color Black = new Color(0, 0, 0);
const Color Red = new Color(255, 0, 0);
const Color Green = new Color(0, 255, 0);
const Color Blue = new Color(0, 0, 255);
const Color Yellow = new Color(255, 255, 0);
const Color Purple = new Color(255, 0, 255);
const Color Orange = new Color(255, 127, 0);