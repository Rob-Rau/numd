module numd.utility;

struct Meshgrid
{
	real[][] X;
	real[][] Y;
}

Meshgrid meshgrid(real[] x, real[] y)
{
	Meshgrid mesh = {X: new real[][](y.length, x.length), Y: new real[][](y.length, x.length)};

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.X[i][j] = x[j];

	for(int i = 0; i < y.length; i++)
		for(int j = 0; j < x.length; j++)
			mesh.Y[i][j] = y[i];

	return mesh;
}

real[] linspace(real start, real end, int points)
{
	real h = (end-start)/(points-1);
	real[] x = new real[points];
	x[0] = start;
	for(int i = 1; i < points; i++)
	{
		x[i] = x[i-1]+h;
	}
	return x;
}