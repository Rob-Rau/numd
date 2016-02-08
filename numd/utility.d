module numd.utility;

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