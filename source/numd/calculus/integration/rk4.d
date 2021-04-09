module numd.calculus.integration.rk4;

/++
	Fourth order Runge-Kutta time integration scheme
+/
class RK4(T, sizediff_t static_size = -1)
{
	static if(static_size == -1) {
		alias StateArray = T[];
	} else {
		alias StateArray = T[static_size];
	}
	
	private StateArray tmp;
	private StateArray k;

	this(size_t data_size)
	{
		if(static_size == -1) {
			tmp = new T[data_size];
			k = new T[data_size];
		}
	}

	@nogc void step(alias func, Args...)(ref StateArray x_next, ref StateArray x, T t, T dt, Args args)
	{
		func(k, x, t, dt, args);
		tmp[] = x[] + 0.5*dt*k[];
		x_next[] = k[];

		func(k, tmp, t + 0.5*dt, dt, args);
		tmp[] = x[] + 0.5*dt*k[];
		x_next[] += 2.0*k[];

		func(k, tmp, t, dt + 0.5*dt, args);
		tmp[] = x[] + dt*k[];
		x_next[] += 2.0*k[];

		func(k, tmp, t + dt, dt, args);
		x_next[] += k[];
		x_next[] *= (dt/6.0);
		x_next[] += x[];
		//x_next[] = x[] + (dt/6.0)*(k1[] + 2.0*k2[] + 2.0*k3[] + k4[]);
	}

	auto integrate(alias func, Args...)(StateArray x_0, T t_end, T dt, Args args) {
		import std.range : enumerate, iota;
		auto time = iota(0, t_end, dt);

		static if(static_size == -1) {
			auto x = new T[][](time.length, tmp.length);
		} else {
			auto x = new StateArray[time.length];
		}

		x[0][] = x_0[];

		foreach(idx, t; time.enumerate[0..$-1]) {
			step!func(x[idx + 1], x[idx], t, dt, args);
		}

		import std.typecons : tuple;
		return tuple!("x", "t")(x, time);
	}
}
