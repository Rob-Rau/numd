module numd.calculus.integration.forwardeuler;

/++
	First order forwared euler integrator
+/
class ForwardEuler(T, sizediff_t static_size = -1)
{
	static if(static_size == -1) {
		alias StateArray = T[];
	} else {
		alias StateArray = T[static_size];
	}
	
	private StateArray x_dot;

	this(size_t data_size) {
		if(static_size == -1) {
			x_dot = new T[data_size];
		}
	}

	@nogc void step(alias func, Args...)(ref StateArray x_next, ref StateArray x, T t, T dt, Args args) {
		func(x_dot, x, t, dt, args);
		x_next[] = x[] + dt*x_dot[];
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

		x[] = x_0[];

		// Time march, ping-ponging which buffer we use for the output.
		foreach(idx, t; time.enumerate[0..$-1]) {
			step!func(x, x, t, dt, args);
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
