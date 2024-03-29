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

	this(size_t data_size) {
		if(static_size == -1) {
			tmp = new T[data_size];
			k = new T[data_size];
		}
	}

	@nogc void step(alias func, Args...)(ref StateArray x_next, ref StateArray x, T t, T dt, Args args) {
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
	}

	@nogc void integrate(alias func, Args...)(StateArray x, StateArray x_0, T t_end, T dt, auto ref Args args) {
		assert(x.ptr != x_0.ptr);

		T t_start, t_stop;
		if(dt > 0) {
			t_start = 0;
			t_stop = t_end;
		} else {
			t_start = t_end + dt;
			t_stop = dt;
		}

		auto time = iota(0, t_end, dt);

		// Time march, ping-ponging which buffer we use for the output.
		foreach(idx, t; time.enumerate[0..$-1]) {
			if(idx % 2 == 0) {
				step!func(x, x_0, t, dt, args);
			} else {
				step!func(x_0, x, t, dt, args);
				
				// If this is the last time step, we want the actual result in x not x_0
				if(idx == time.length - 1) {
					x[] = x_0[];
				}
			}
		}
	}

	auto integrate(alias func, Args...)(StateArray x_0, T t_end, T dt, auto ref Args args) {
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
