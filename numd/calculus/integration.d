module numd.calculus.integration;

struct IntegrationResults(size_t dimensions)
{
	@nogc
	{
		import std.experimental.allocator.mallocator : Mallocator;
		import core.stdc.stdio : printf;

		// Temporal points
		double[] T = null;
		// Spatial points
		double[dimensions][] X = null;

		private int* referenceCount = null;

		this(this)
		{
			if(referenceCount !is null)
			{
				(*referenceCount)++;
			}
			else
			{
				printf("Error IntegrationResults referenceCount somehow null in postblit\n");
			}
		}

		this(size_t n)
		{
			T = cast(double[])Mallocator.instance.allocate(n*double.sizeof);
			X = cast(double[dimensions][])Mallocator.instance.allocate(dimensions*n*double.sizeof);

			T[] = 0;
			foreach(ref subArr; X)
			{
				subArr[] = 0;
			}

			referenceCount = cast(int*)Mallocator.instance.allocate(int.sizeof);
			(*referenceCount) = 1;
		}

		~this()
		{
			if(referenceCount !is null)
			{
				if((*referenceCount) == 1)
				{
					Mallocator.instance.deallocate(T);
					Mallocator.instance.deallocate(X);
					Mallocator.instance.deallocate(referenceCount[0..int.sizeof]);
				}
				else
				{
					(*referenceCount)--;
				}
			}
			else
			{
				printf("Error IntegrationResults referenceCount somehow null in destructor\n");
			}
		}
	}
}