@compute(CompileFor.hostAndDevice) module numd.math;

import numd.math.cuda;

import ldc.dcompute;

public enum double fpTol = 1.0e-16;

pragma(inline, true)
@trusted @nogc T abs(T)(auto ref T t) {
	
	if(__dcompute_reflect(ReflectTarget.Host)) {
		//static import std.math;
		static import core.stdc.math;
		return core.stdc.math.fabs(t);
	} else if(__dcompute_reflect(ReflectTarget.CUDA)) {
		import numd.math.cuda : fabs_d, fabs_f;
		static if(is(T == float)) {
			return fabs_f(t);
		} else static if(is(T == double)) {
			return fabs_d(t);
		}
	} else if(__dcompute_reflect(ReflectTarget.OpenCL)) {
		assert(0);
	}
	assert(0);
}

pragma(inline, true)
@trusted @nogc T sqrt(T)(auto ref T t) {
	pragma(inline, true);
	if(__dcompute_reflect(ReflectTarget.Host)) {
		static import std.math;
		return std.math.sqrt(t);
	} else if(__dcompute_reflect(ReflectTarget.CUDA)) {
		import numd.math.cuda : fabs_d, fabs_f;
		static if(is(T == float)) {
			return sqrt_f(t);
		} else static if(is(T == double)) {
			return sqrt_d(t);
		}
	} else if(__dcompute_reflect(ReflectTarget.OpenCL)) {
		assert(0);
	}
	assert(0);
}

pragma(inline, true)
@trusted @nogc T round(T)(auto ref T t) {
	if(__dcompute_reflect(ReflectTarget.Host)) {
		static import std.math;
		return std.math.round(t);
	} else if(__dcompute_reflect(ReflectTarget.CUDA)) {
		import numd.math.cuda : fabs_d, fabs_f;
		static if(is(T == float)) {
			return round_f(t);
		} else static if(is(T == double)) {
			return round_d(t);
		}
	} else if(__dcompute_reflect(ReflectTarget.OpenCL)) {
		assert(0);
	}
	assert(0);
}
