module numd.utility;

struct Meshgrid(T = double)
{
	T[][] X;
	T[][] Y;
}

Meshgrid meshgrid(T = double)(T[] x, T[] y)
{
	Meshgrid mesh = {X: new T[][](y.length, x.length), Y: new T[][](y.length, x.length)};

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.X[i][j] = x[j];

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.Y[i][j] = y[i];

	return mesh;
}

T[] linspace(T = double)(double start, double end, int points)
{
	//real h = (end-start)/(points-1);
	T h = (end-start)/(points);
	T[] x = new T[points];
	x[0] = start;
	for(int i = 1; i < points; i++)
	{
		x[i] = x[i-1]+h;
	}
	return x;
}