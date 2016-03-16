module numd.optimization.Gradient;

interface IGradient
{
	double[] Gradient(in double[] point);
}