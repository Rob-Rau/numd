module numd.linearalgebra.matrixarray;

import numd.linearalgebra.matrix;
import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.meta;
import std.range;
import std.stdio;
import std.traits;
import std.typecons;

unittest
{
	immutable size_t arraySize = 10_000_000;
	//immutable size_t arraySize = 10000;
	immutable size_t vecSize = 5;
	auto n0 = new double[arraySize];
	n0[] = 0.5;
	auto n1 = new double[arraySize];
	n1[] = 0.5;
	auto dqdx = VectorArray!vecSize(arraySize, 0.1);
	auto dqdy = VectorArray!vecSize(arraySize, 1.0);
	//auto A00 = MatrixArray!(vecSize, vecSize)(arraySize, 2.0);
	//auto A01 = MatrixArray!(vecSize, vecSize)(arraySize, 3.0);
	//auto A10 = MatrixArray!(vecSize, vecSize)(arraySize, 4.0);
	//auto A11 = MatrixArray!(vecSize, vecSize)(arraySize, 5.0);
	auto flux = VectorArray!vecSize(arraySize, 5.0);

	@nogc diffuseFlux()
	{
		/+
		immutable size_t chunkSize = 65;
		for(size_t i = 0; i < flux.length; i += chunkSize)
		{
			size_t sliceEnd = i + chunkSize > flux.length ? flux.length : i+chunkSize;
			auto G0 = A00[i..sliceEnd]*dqdx[i..sliceEnd] + A01[i..sliceEnd]*dqdy[i..sliceEnd];
			//flux = G0;
			auto G1 = A10[i..sliceEnd]*dqdx[i..sliceEnd] + A11[i..sliceEnd]*dqdy[i..sliceEnd];
			//flux[i..sliceEnd] = n0[i..sliceEnd]*G0 + n1[i..sliceEnd]*G1;
			flux[i..sliceEnd] = n0[i..sliceEnd]*G0 + n1[i..sliceEnd]*G1 + n0[i..sliceEnd]*G1 + n1[i..sliceEnd]*G0 + n0[i..sliceEnd]*(G0 + G1);
		}+/

		immutable size_t chunkSize = 136;
		for(size_t i = 0; i < flux.length; i += chunkSize)
		{
			size_t sliceEnd = i + chunkSize > flux.length ? flux.length : i+chunkSize;
			auto A00 = dqdy[i..sliceEnd]*dqdx[i..sliceEnd].transpose;
			auto A01 = dqdx[i..sliceEnd]*dqdy[i..sliceEnd].transpose;
			auto A10 = dqdy[i..sliceEnd]*dqdx[i..sliceEnd].transpose;
			auto A11 = dqdx[i..sliceEnd]*dqdy[i..sliceEnd].transpose;

			auto G0 = A00*dqdx[i..sliceEnd] + A01*dqdy[i..sliceEnd];
			auto G1 = A10*dqdx[i..sliceEnd] + A11*dqdy[i..sliceEnd];

			auto G2 = A00*dqdy[i..sliceEnd] + A01*dqdx[i..sliceEnd];
			auto G3 = A10*dqdy[i..sliceEnd] + A11*dqdx[i..sliceEnd];

			flux[i..sliceEnd] = n0[i..sliceEnd]*G0 + n1[i..sliceEnd]*G1 + n0[i..sliceEnd]*G1 + n1[i..sliceEnd]*G0 + n0[i..sliceEnd]*(G2 + G3) - n1[i..sliceEnd]*(G2 - G3);
		}
		/+
		auto A00 = dqdy*dqdx.transpose;
		//flux = A00*dqdy;
		//pragma(msg, typeof(A00));
		auto A01 = dqdx*dqdy.transpose;
		auto A10 = dqdy*dqdx.transpose;
		auto A11 = dqdx*dqdy.transpose;

		auto G0 = A00*dqdx + A01*dqdy;
		//flux = G0;
		auto G1 = A10*dqdx + A11*dqdy;
		//flux = n0*G0 + n1*G1;
		flux = n0*G0 + n1*G1 + n0*G1 + n1*G0 + n0*(G0 + G1);
		+/
		//flux = n0*(A00*dqdx + A01*dqdy) + n1*(A10*dqdx + A11*dqdy) + n0*(A10*dqdx + A11*dqdy) + n1*(A00*dqdx + A01*dqdy);
	}

	auto dqdxo = new Vector!vecSize[arraySize];//, 0.1);
	auto dqdyo = new Vector!vecSize[arraySize];//, 1.0);
	//auto A00o = new Matrix!(vecSize, vecSize)[arraySize];//, 2.0);
	//auto A01o = new Matrix!(vecSize, vecSize)[arraySize];//, 3.0);
	//auto A10o = new Matrix!(vecSize, vecSize)[arraySize];//, 4.0);
	//auto A11o = new Matrix!(vecSize, vecSize)[arraySize];//, 5.0);
	auto fluxo = new Vector!vecSize[arraySize];//, 5.0);

	for(size_t i = 0; i < arraySize; i++)
	{
		dqdxo[i] = 0.1;
		dqdyo[i] = 1.0;
		//A00o[i] = 2;
		//A01o[i] = 3;
		//A10o[i] = 4;
		//A11o[i] = 5;
		fluxo[i] = 5;
	}

	@nogc diffuseFluxOld()
	{
		for(size_t i = 0; i < arraySize; i++)
		{
			//auto G0 = A00o[i]*dqdxo[i] + A01o[i]*dqdyo[i];
			//auto G1 = A10o[i]*dqdxo[i] + A11o[i]*dqdyo[i];
			//fluxo[i] = n0[i]*G0 + n1[i]*G1 + n0[i]*G1 + n1[i]*G0;
			//fluxo[i] = n0[i]*G0 + n1[i]*G1 + n0[i]*G1 + n1[i]*G0 + n0[i]*(G0 + G1);
			auto A00 = dqdyo[i]*dqdxo[i].transpose;
			auto A01 = dqdxo[i]*dqdyo[i].transpose;
			auto A10 = dqdyo[i]*dqdxo[i].transpose;
			auto A11 = dqdxo[i]*dqdyo[i].transpose;
			auto G0 = A00*dqdxo[i] + A01*dqdyo[i];
			auto G1 = A10*dqdxo[i] + A11*dqdyo[i];

			auto G2 = A00*dqdyo[i] + A01*dqdxo[i];
			auto G3 = A10*dqdyo[i] + A11*dqdxo[i];

			fluxo[i] = n0[i]*G0 + n1[i]*G1 + n0[i]*G1 + n1[i]*G0 + n0[i]*(G2 + G3) - n1[i]*(G2 - G3);
			//fluxo[i] = n0[i]*G0 + n1[i]*G1;
			//pragma(msg, typeof(A00));
			//fluxo[i] = A00*dqdyo[i];
		}
	}

	writeln("Running diffuse flux test");
	import std.datetime.stopwatch : benchmark;
	//uint compTimes = 5_000_000;
	uint compTimes = 1000;
	//auto r = benchmark!nogcComputeAdd(compTimes);
	//writeln("Addition computation time ", r[0]);

	writeln;
	writeln("Running new diffuse flux test");
	auto r = benchmark!diffuseFlux(compTimes);
	writeln("New diffuse flux time ", r[0]);

	writeln;
	writeln("Running old diffuse flux test");
	r = benchmark!diffuseFluxOld(compTimes);
	writeln("Old diffuse flux time ", r[0]);

	foreach(k; 0..arraySize)
	{
		auto qMat = flux.matrixAt(k);
		static foreach(i; 0..vecSize)
		{
			static foreach(j; 0..1)
			{
				assert(abs(qMat[i,j] - fluxo[k][i,j]) < 1.0e-15,
					"Differing results in matrix "~k.to!string~":"~
					"\nQ[k]["~i.to!string~", "~j.to!string~"] = "~qMat[i,j].to!string~
					"\nQold[k]["~i.to!string~", "~j.to!string~"] = "~fluxo[k][i,j].to!string~
					"\nQ[k] = "~qMat.to!string~
					"\nQold[k] = "~fluxo[k].to!string~
					"\ndiff = "~(qMat[i,j] - fluxo[k][i,j]).to!string
				);
			}
		}
	}
}
/+
unittest
{
	import std.stdio : writeln;
	
	writeln;
	writeln("=======================");
	writeln("Allocating small arrays");
	size_t arraySize = 100;
	//size_t arraySize = 10_000_000;
	immutable size_t vecSize = 5;
	auto ma1 = MatrixArray!(vecSize,vecSize)(arraySize, 1.0);
	auto ma2 = MatrixArray!(vecSize,vecSize)(arraySize, 2.0);
	auto ma4 = MatrixArray!(vecSize,vecSize)(arraySize, 3.0);
	auto ma5 = MatrixArray!(vecSize,vecSize)(arraySize, 4.0);
	
	/+
	@nogc void nogcComputeAdd()
	{
		auto ma3 = ma1 + ma2;
		ma4 = ma3;
		ma5 = ma1 + ma2 + ma4;
	}
	+/
	@nogc void nogcComputeMult()
	{
		pragma(inline, true);
		/+
		for(size_t i = 0; i < ma5.length; i += 100)
		{
			size_t sliceEnd = i + 100 > ma5.length ? ma5.length : i+100;
			ma5[i..sliceEnd] = ma4[i..sliceEnd] + ma1[i..sliceEnd]*ma2[i..sliceEnd];
		}+/
		ma5 = ma4 + ma1*ma2;
	}

	auto ma1o = new Matrix!(vecSize,vecSize)[arraySize];
	auto ma2o = new Matrix!(vecSize,vecSize)[arraySize];
	auto ma4o = new Matrix!(vecSize,vecSize)[arraySize];
	auto ma5o = new Matrix!(vecSize,vecSize)[arraySize];
	for(size_t i = 0; i < arraySize; i++)
	{
		ma1o[i] = 1.0;
		ma2o[i] = 2.0;
		ma4o[i] = 3.0;
		ma5o[i] = 4.0;
	}
	@nogc void nogcComputeMultOld()
	{
		pragma(inline, true);
		for(size_t i = 0; i < arraySize; i++)
		{
			ma5o[i] = ma4o[i] + ma1o[i]*ma2o[i];
		}
	}

	writeln("Running matrix array test");
	import std.datetime.stopwatch : benchmark;
	//uint compTimes = 5_000_000;
	uint compTimes = 100000;
	//auto r = benchmark!nogcComputeAdd(compTimes);
	//writeln("Addition computation time ", r[0]);

	writeln;
	writeln("Running new mat mult test");
	auto r = benchmark!nogcComputeMult(compTimes);
	writeln("New Multiplication computation time ", r[0]);

	writeln;
	writeln("Running old mat mult test");
	r = benchmark!nogcComputeMultOld(compTimes);
	writeln("Old Multiplication computation time ", r[0]);

	foreach(k; 0..arraySize)
	{
		auto qMat = ma5.matrixAt(k);
		static foreach(i; 0..vecSize)
		{
			static foreach(j; 0..vecSize)
			{
				assert(abs(qMat[i,j] - ma5o[k][i,j]) < 1.0e-16,
					"Differing results in matrix "~k.to!string~":"~
					"\nQ[k]["~i.to!string~", "~j.to!string~"] = "~qMat[i,j].to!string~
					"\nQold[k]["~i.to!string~", "~j.to!string~"] = "~ma5o[k][i,j].to!string~
					"\nQ[k] = "~qMat.to!string~
					"\nQold[k] = "~ma5o[k].to!string~
					"\ndiff = "~(qMat[i,j] - ma5o[k][i,j]).to!string
				);
			}
		}
	}
}+/
/+
unittest
{
	import std.stdio : writeln;

	writeln;
	writeln("=======================");
	writeln("Allocating big arrays");
	size_t arraySize = 100_000_000;
	auto ma1 = MatrixArray!(4,4)(arraySize);
	auto ma2 = MatrixArray!(4,4)(arraySize);
	auto ma4 = MatrixArray!(4,4)(arraySize);
	auto ma5 = MatrixArray!(4,4)(arraySize);
	
	@nogc void nogcComputeAdd()
	{
		auto ma3 = ma1 + ma2;
		ma4 = ma3;
		ma5 = ma1 + ma2 + ma4;
	}

	@nogc void nogcComputeMult()
	{
		ma4 = ma1*ma2;
	}

	writeln("Running matrix array test");
	import std.datetime.stopwatch : benchmark;
	uint compTimes = 2;
	auto r = benchmark!nogcComputeAdd(compTimes);
	writeln("Addition computation time ", r[0]);

	r = benchmark!nogcComputeMult(compTimes);
	writeln("Multiplication computation time ", r[0]);

}
+/
/+
unittest
{
	writeln;
	writeln("=======================");
	writeln("Allocating integrator arrays");
	immutable size_t vecSize = 4;
	//immutable size_t arraySize = 1_000_000;
	immutable size_t arraySize = 1_000;
	auto dts = new double[arraySize];

	double dt = 0.000002;
	dts[] = dt;

	auto Q = VectorArray!vecSize(arraySize);
	auto k1 = VectorArray!vecSize(arraySize);
	auto k2 = VectorArray!vecSize(arraySize);
	auto k3 = VectorArray!vecSize(arraySize);
	auto k4 = VectorArray!vecSize(arraySize);

	static foreach(i; 0..vecSize)
	{
		Q.assignAllAt!i(1.0 + 0.1*i.to!double);
		k1.assignAllAt!i(2.0 + 0.1*i.to!double);
		k2.assignAllAt!i(3.0 + 0.1*i.to!double);
		k3.assignAllAt!i(4.0 + 0.1*i.to!double);
		k4.assignAllAt!i(5.0 + 0.1*i.to!double);
	}

	void rk4new()
	{
		Q = Q + dts*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
	}

	void rk4newglobal()
	{
		Q = Q + dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
	}

	auto k1o = new Vector!vecSize[arraySize];
	auto k2o = new Vector!vecSize[arraySize];
	auto k3o = new Vector!vecSize[arraySize];
	auto k4o = new Vector!vecSize[arraySize];
	auto Qo = new Vector!vecSize[arraySize];

	for(size_t i = 0; i < arraySize; i++)
	{
		for(size_t j = 0; j < vecSize; j++)
		{
			Qo[i][j] = 1.0 + 0.1*j;
		}
		for(size_t j = 0; j < vecSize; j++)
		{
			k1o[i][j] = 2.0 + 0.1*j;
		}
		for(size_t j = 0; j < vecSize; j++)
		{
			k2o[i][j] = 3.0 + 0.1*j;
		}
		for(size_t j = 0; j < vecSize; j++)
		{
			k3o[i][j] = 4.0 + 0.1*j;
		}
		for(size_t j = 0; j < vecSize; j++)
		{
			k4o[i][j] = 5.0 + 0.1*j;
		}
	}

	void rk4old()
	{
		for(size_t i = 0; i < arraySize; i++)
		{
			Qo[i] += (dts[i]/6.0)*(k1o[i] + 2.0*k2o[i] + 2.0*k3o[i] + k4o[i]);
		}
	}

	void rk4oldglobal()
	{
		for(size_t i = 0; i < arraySize; i++)
		{
			Qo[i] += (dt/6.0)*(k1o[i] + 2.0*k2o[i] + 2.0*k3o[i] + k4o[i]);
		}
	}

	import std.datetime.stopwatch : benchmark;
	//uint compTimes = 500;
	uint compTimes = 5_000_000;

	writeln("Running rk4 new test");
	auto r = benchmark!rk4new(compTimes);
	writeln("rk4new computation time ", r[0]);

	writeln;
	writeln("Running rk4 old test");
	r = benchmark!rk4old(compTimes);
	writeln("rk4old computation time ", r[0]);

	foreach(k; 0..arraySize)
	{
		auto qMat = Q.matrixAt(k);
		static foreach(i; 0..vecSize)
		{
			static foreach(j; 0..1)
			{
				assert(abs(qMat[i,j] - Qo[k][i,j]) < 1.0e-16,
					"Differing results in matrix "~k.to!string~":"~
					"\n\tQ[k]["~i.to!string~", "~j.to!string~"] = "~qMat[i,j].to!string~
					"\n\tQold[k]["~i.to!string~", "~j.to!string~"] = "~Qo[k][i,j].to!string~
					"\n\tQ[k] = "~qMat.to!string~
					"\n\tQold[k] = "~Qo[k].to!string~
					"\n\tdiff = "~(qMat[i,j] - Qo[k][i,j]).to!string
				);
			}
		}
	}

	static foreach(i; 0..vecSize)
	{
		Q.assignAllAt!i(1.0 + 0.1*i.to!double);
	}

	for(size_t i = 0; i < arraySize; i++)
	{
		for(size_t j = 0; j < vecSize; j++)
		{
			Qo[i][j] = 1.0 + 0.1*j;
		}
	}

	writeln;
	writeln("Running rk4 new global test");
	r = benchmark!rk4newglobal(compTimes);
	writeln("rk4new global computation time ", r[0]);

	writeln;
	writeln("Running rk4 old global test");
	r = benchmark!rk4oldglobal(compTimes);
	writeln("rk4old global computation time ", r[0]);

	foreach(k; 0..arraySize)
	{
		auto qMat = Q.matrixAt(k);
		static foreach(i; 0..vecSize)
		{
			static foreach(j; 0..1)
			{
				assert(abs(qMat[i,j] - Qo[k][i,j]) < 1.0e-16,
					"Differing results in matrix "~k.to!string~":"~
					"\n\tQ[k]["~i.to!string~", "~j.to!string~"] = "~qMat[i,j].to!string~
					"\n\tQold[k]["~i.to!string~", "~j.to!string~"] = "~Qo[k][i,j].to!string~
					"\n\tQ[k] = "~qMat.to!string~
					"\n\tQold[k] = "~Qo[k].to!string~
					"\n\tdiff = "~(qMat[i,j] - Qo[k][i,j]).to!string
				);
			}
		}
	}
}+/

auto assumeNogc(alias Func, T...)(T xs) @nogc {
	import std.traits;
	static auto assumeNogcPtr(T)(T f) if (
		isFunctionPointer!T || isDelegate!T
	) {
  	    enum attrs = functionAttributes!T | FunctionAttribute.nogc;
    	    return cast(SetFunctionAttributes!(T, functionLinkage!T, attrs)) f;
	}
	return assumeNogcPtr(&Func!T)(xs);
}

private struct MatrixArrayResult(size_t r, size_t c, T = double, alias operation, Args...)
{
	alias ThisType = typeof(this);

	//pragma(msg, "Tuple!"~Args.stringof~" args;");
	mixin("Tuple!"~Args.stringof~" args;");
	
	template hasLength(L)
	{
		enum bool hasLength = hasMember!(L, "length");
	}

	immutable size_t length;
	@trusted @nogc this(Tuple!Args args)
	{
		this.args = args;
		alias lenArgs = Filter!(hasLength, AliasSeq!Args);
		alias lenIndex = staticIndexOf!(lenArgs[0], AliasSeq!Args);
		this.length = args[lenIndex].length;
	}

	@trusted @nogc auto opBinary(string op, size_t ic)(ref inout MatrixArray!(c, ic, T) inRhs)
		if(op == "*")
	{
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){

			static foreach(a; aliasSeqOf!(iota(0, c)))
			{
				mixin(`enum string code`~a.to!string~` = 
					"{
						T[`~destSize~`]	temp`~(level+1).to!string~`;
						`~operation!("=", level+1, i, to!string(a), args[0..this.args.length], range, destSize)~`
						temp`~(level+1).to!string~`[] *= finalResult.args[`~args[$-1]~`].q`~a.to!string~j~`[`~range~`];
						temp`~level.to!string~`[] `~(a == 0? cop : "+=")~` temp`~(level+1).to!string~`[];
					}";
				`);
			}
			mixin(`enum string code = "return \""~code0`~iota(1, c).map!(a => "~code"~a.to!string).join~`~"\";";`);
			mixin(code);
		}

		return MatrixArrayResult!(r, ic, T, comp, Args, MatrixArray!(c, ic, T))(tuple(this.args.expand, cast(MatrixArray!(c, ic, T))inRhs));
	}

	@trusted @nogc auto opBinary(string op, size_t ic)(inout MatrixArray!(c, ic, T) inRhs)
		if(op == "*")
	{
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){

			static foreach(a; aliasSeqOf!(iota(0, c)))
			{
				mixin(`enum string code`~a.to!string~` = 
					"{
						T[`~destSize~`]	temp`~(level+1).to!string~`;
						`~operation!("=", level+1, i, to!string(a), args[0..this.args.length], range, destSize)~`
						temp`~(level+1).to!string~`[] *= finalResult.args[`~args[$-1]~`].q`~a.to!string~j~`[`~range~`];
						temp`~level.to!string~`[] `~(a == 0? cop : "+=")~` temp`~(level+1).to!string~`[];
					}";
				`);
			}
			mixin(`enum string code = "return \""~code0`~iota(1, c).map!(a => "~code"~a.to!string).join~`~"\";";`);
			mixin(code);
		}

		return MatrixArrayResult!(r, ic, T, comp, Args, MatrixArray!(c, ic, T))(tuple(this.args.expand, cast(MatrixArray!(c, ic, T))inRhs));
	}
	
	@trusted @nogc auto opBinary(string op)(ref inout MatrixArray!(r, c, T) inRhs)
		if((op == "+") || (op == "-"))
	{
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			enum string nextLevel = (level+1).to!string;
			enum string code = `
				return "
					{
						T[`~destSize~`] temp`~nextLevel~`;
						`~operation!(cop, level+1, i, j, args[0..this.args.length], range, destSize)~`
						temp`~nextLevel~`[] `~op~`= finalResult.args[`~args[$-1]~`].q`~i~j~`[`~range~`];
						temp`~level.to!string~`[] `~cop~` temp`~nextLevel~`[];
					}
				";
			`;
			mixin(code);
		}

		return MatrixArrayResult!(r, c, T, comp, Args, MatrixArray!(r, c, T))(tuple(this.args.expand, cast(MatrixArray!(r, c, T))inRhs));
	}

	@trusted @nogc auto opBinary(string op, alias inOperation, inArgs...)(MatrixArrayResult!(r, c, T, inOperation, inArgs) inRhs)
		if((op == "+") || (op == "-"))
	{
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			enum string nextLevel = (level+1).to!string;
			enum string code = `
				return "
					{
						T[`~destSize~`] temp`~nextLevel~`;
						`~operation!("=", level+1, i, j, args[0..this.args.length], range, destSize)~`
						temp`~level.to!string~`[] `~cop~` temp`~nextLevel~`[];
						`~inOperation!("=", level+1, i, j, args[this.args.length..$], range, destSize)~`
						temp`~level.to!string~`[] `~op~`= temp`~nextLevel~`[];
					}

				";
			`;
			debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			mixin(code);
		}

		return MatrixArrayResult!(r, c, T, comp, Args, inArgs)(tuple(this.args.expand, inRhs.args.expand));
	}
	
	@trusted @nogc auto opBinaryRight(string op)(inout T[] inLhs)
		if(op == "*")
	{
		assert(this.length == inLhs.length, "Matrix arrays are not the same size");
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			enum string code = `
				return "
					`~operation!(cop, level, i, j, args[0..this.args.length], range, destSize)~`
					temp`~level.to!string~`[] `~op~`= finalResult.args[`~args[$-1]~`][`~range~`];
				";
			`;
			debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			//debug pragma(msg, "args = "~args.join(","));
			//debug pragma(msg, code~"\n");
			mixin(code);
		}

		return MatrixArrayResult!(r, c, T, comp, Args, T[])(tuple(this.args.expand, cast(T[])inLhs));
	}

	@trusted @nogc auto opBinaryRight(string op)(ref inout T inLhs)
		if(op == "*")
	{
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			enum string code = `
				return "
					`~operation!(cop, level, i, j, args[0..this.args.length], range, destSize)~`
					temp`~level.to!string~`[] `~op~`= finalResult.args[`~args[$-1]~`];
				";
			`;
			debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			//debug pragma(msg, "args = "~args.join(","));
			//debug pragma(msg, code~"\n");
			mixin(code);
		}

		return MatrixArrayResult!(r, c, T, comp, Args, T)(tuple(this.args.expand, cast(T)inLhs));
	}

	@trusted @nogc auto opBinary(string op)(inout T inRhs)
		if((op == "*") || (op == "/"))
	{
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			enum string code = `
				return "
					{
						T[`~destSize~`] temp`~(level+1).to!string~`;
						`~operation!(cop, level+1, i, j, args[0..this.args.length], range, destSize)~`
						temp`~(level+1).to!string~`[] `~op~`= finalResult.args[`~args[$-1]~`];
						temp`~level.to!string~`[] `~cop~` temp`~(level+1).to!string~`[];
					}
				";
			`;
			mixin(code);
		}

		return MatrixArrayResult!(r, c, T, comp, Args, T)(tuple(this.args.expand, cast(T)inRhs));
	}
}

alias VectorArray(size_t l, T = double) = MatrixArray!(l, 1, T);

struct MatrixArray(size_t r, size_t c, T = double)
{
	alias ThisType = typeof(this);
	immutable size_t length;

	static foreach(i; 0..r)
	{
		static foreach(j; 0..c)
		{
			mixin("private T[] q"~i.to!string~j.to!string~";");
		}
	}

	this(size_t arraySize)
	{
		length = arraySize;
		static foreach(i; 0..r)
		{
			static foreach(j; 0..c)
			{
				mixin("q"~i.to!string~j.to!string~" = new T[arraySize];");
			}
		}
	}

	this(size_t arraySize, T initValue)
	{
		length = arraySize;
		
		static foreach(i; 0..r)
		{
			static foreach(j; 0..c)
			{
				mixin("q"~i.to!string~j.to!string~" = new T[arraySize];");
				mixin("q"~i.to!string~j.to!string~"[] = initValue;");
			}
		}
	}

	@nogc
	{
		private this(size_t arraySize, bool dummy)
		{
			length = arraySize;
		}

		@trusted MatrixArray!(c, r, T) transpose()
		{
			auto ret = MatrixArray!(c, r, T)(this.length, false);
			static foreach(i; aliasSeqOf!(iota(0, c)))
			{
				static foreach(j; aliasSeqOf!(iota(0, r)))
				{
					mixin("ret.q"~i.to!string~j.to!string~" = this.q"~j.to!string~i.to!string~";");
				}
			}
			return ret;
		}

		@trusted ThisType opSlice(size_t ii, size_t ij)
		{
			auto ret = ThisType(ij - ii, false);
			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					mixin("ret.q"~i.to!string~j.to!string~" = this.q"~i.to!string~j.to!string~"[ii..ij];");
				}
			}
			return ret;
		}

		@trusted void assignAllAt(size_t i, size_t j = 0)(inout T val)
		{
			mixin("q"~i.to!string~j.to!string~"[] = val;");
		}

		@trusted void assignAll(inout T val)
		{
			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					mixin("q"~i.to!string~j.to!string~"[] = val;");
				}
			}
		}

		@trusted void assignAll(ref inout Matrix!(r, c, T) mat)
		{
			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					mixin("q"~i.to!string~j.to!string~"[] = mat["~i.to!string~","~j.to!string~"];");
				}
			}
		}

		@trusted auto matrixAt(size_t idx)
		{
			Matrix!(r, c, T) newMat;
			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					mixin("newMat["~i.to!string~","~j.to!string~"] = q"~i.to!string~j.to!string~"[idx];");
				}
			}
			return newMat;
		}

		@trusted ref ThisType opAssign(alias operation, Args...)(MatrixArrayResult!(r, c, T, operation, Args) finalResult)
		{
			import core.stdc.stdio : printf;
			enum size_t stride16 = 16;
			enum size_t stride8 = 8;
			enum size_t stride4 = 4;
			enum size_t stride2 = 2;
			size_t k;
			size_t iterations;
			size_t remaining;
			pragma(inline, true);
			//enum size_t i = 0;
			//enum size_t j = 0;
			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					k = 0;
					iterations = finalResult.length/stride16;
					remaining = finalResult.length%stride16;
					//debug assumeNogc!writeln("iterations: ", iterations, " remaining: ", remaining);
					/+
					while(iterations)
					{
						T[stride16] temp0;
						debug pragma(msg, "finalResult.args.length = "~finalResult.args.length.to!string);
						enum string code = operation!("=", 0, to!string(i), to!string(j), iota(0, finalResult.args.length).map!(a => a.to!string).array, "k..k+"~to!string(stride16), to!string(stride16));
						debug pragma(msg, code);
						mixin(code);
						mixin("this.q"~i.to!string~j.to!string~"[k..k+stride16] = temp0[];");
						iterations--;
						k+=stride16;
					}
					
					if(remaining > stride8)
					{+/
						iterations = finalResult.length/stride8;
						remaining = finalResult.length%stride8;
						while(iterations)
						{
							T[stride8] temp0;
							enum string code = operation!("=", 0, to!string(i), to!string(j), iota(0, finalResult.args.length).map!(a => a.to!string).array, "k..k+"~to!string(stride8), to!string(stride8));
							mixin(code);
							//debug pragma(msg, code);
							mixin("this.q"~i.to!string~j.to!string~"[k..k+stride8] = temp0[];");
							iterations--;
							k+=stride8;
						}
					/+}

					if(remaining > stride4)
					{
						iterations = remaining/stride4;
						remaining = remaining%stride4;
						while(iterations)
						{
							T[stride4] temp0;
							enum string code = operation!("=", 0, to!string(i), to!string(j), iota(0, finalResult.args.length).map!(a => a.to!string).array, "k..k+"~to!string(stride4), to!string(stride4));
							mixin(code);
							//debug pragma(msg, code);
							mixin("this.q"~i.to!string~j.to!string~"[k..k+stride4] = temp0[0..stride4];");
							iterations--;
							k+=stride4;
						}
					}
					
					if(remaining > stride2)
					{
						iterations = remaining/stride2;
						remaining = remaining%stride2;
						while(iterations)
						{
							T[stride2] temp0;
							enum string code = operation!("=", 0, to!string(i), to!string(j), iota(0, finalResult.args.length).map!(a => a.to!string).array, "k..k+"~to!string(stride2), to!string(stride2));
							mixin(code);
							//debug pragma(msg, code);
							mixin("this.q"~i.to!string~j.to!string~"[k..k+stride2] = temp0[0..stride2];");
							iterations--;
							k+=stride2;
						}
					}+/
					
					while(remaining)
					{
						T[1] temp0;
						enum string code = operation!("=", 0, to!string(i), to!string(j), iota(0, finalResult.args.length).map!(a => a.to!string).array, "k..k+"~to!string(1), "1");
						mixin(code);
						//debug pragma(msg, code);
						mixin("this.q"~i.to!string~j.to!string~"[k..k+1] = temp0[0];");
						remaining--;
						k++;
					}
				}
			}
			return this;
		}
		
		@trusted auto opBinary(string op)(ref inout ThisType inRhs)
			if(op == "+" || op == "-")
		{
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
				enum string code = `
					return "
						temp`~level.to!string~`[] `~cop~` finalResult.args[`~args[0]~`].q`~i~j~`[`~range~`] `~op~` finalResult.args[`~args[1]~`].q`~i~j~`[`~range~`];
					";
				`;
				mixin(code);
			}

			return MatrixArrayResult!(r, c, T, comp, ThisType, ThisType)(tuple(this, cast(ThisType)inRhs));
		}

		@trusted auto opBinary(string op, alias operation, Args...)(MatrixArrayResult!(r, c, T, operation, Args) inRhs)
			if(op == "+" || op == "-")
		{
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
				enum string code = `
					return "
						{
							T[`~destSize~`] temp`~(level + 1).to!string~`;
							`~operation!(cop, level+1, i, j, args[0..inRhs.args.length], range, destSize)~`
							temp`~(level+1).to!string~`[] `~op~`= finalResult.args[`~args[$-1]~`].q`~i~j~`[`~range~`];
							temp`~level.to!string~`[] `~cop~` temp`~(level+1).to!string~`[];
						}
					";
				`;
				mixin(code);
			}

			return MatrixArrayResult!(r, c, T, comp, Args, ThisType)(tuple(inRhs.args.expand, this));
		}

		@trusted auto opBinary(string op, size_t ic)(ref inout MatrixArray!(r, ic, T) inRhs)
			if(op == "*")
		{
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destRange)(){
				enum string code = `
					return "
						temp`~level.to!string~`[] `~cop~` finalResult.args[`~args[0]~`].q`~i~`0[`~range~`] `~op~` finalResult.args[`~args[1]~`].q0`~j~`[`~range~`];
						`~iota(1, c).map!(
							a => "\ntemp"~level.to!string~"[] += finalResult.args["~args[0]~"].q"~i~to!string(a)~"["~range~"] "~op~" finalResult.args["~args[1]~"].q"~to!string(a)~j~"["~range~"];"
						).join~`
					";
				`;
				debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
				//debug pragma(msg, code);
				mixin(code);
			}

			return MatrixArrayResult!(r, ic, T, comp, ThisType, MatrixArray!(r, ic, T))(tuple(this, cast(MatrixArray!(r, ic, T))inRhs));
		}

		@trusted auto opBinary(string op, size_t ic)(inout MatrixArray!(c, ic, T) inRhs)
			if(op == "*")
		{
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destRange)(){
				enum string code = `
					return "
						temp`~level.to!string~`[] `~cop~` finalResult.args[`~args[0]~`].q`~i~`0[`~range~`] `~op~` finalResult.args[`~args[1]~`].q0`~j~`[`~range~`];
						`~iota(1, c).map!(
							a => "\ntemp"~level.to!string~"[] += finalResult.args["~args[0]~"].q"~i~to!string(a)~"["~range~"] "~op~" finalResult.args["~args[1]~"].q"~to!string(a)~j~"["~range~"];"
						).join~`
					";
				`;
				//pragma(msg, code);
				mixin(code);
			}

			return MatrixArrayResult!(r, ic, T, comp, ThisType, MatrixArray!(c, ic, T))(tuple(this, cast(MatrixArray!(c, ic, T))inRhs));
		}

		@trusted auto opBinaryRight(string op)(inout T inLhs)
			if(op == "*" || op == "/")
		{
			@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destRange)(){
				enum string code = `
					return "
						temp`~level.to!string~`[] `~cop~` finalResult.args[`~args[0]~`] `~op~` finalResult.args[`~args[1]~`].q`~i~j~`[`~range~`];
					";
				`;

				mixin(code);
			}

			return MatrixArrayResult!(r, c, T, comp, T, ThisType)(tuple(cast(T)inLhs, this));
		}

		void opOpAssign(string op, alias operation, Args...)(MatrixArrayResult!(r, c, T, operation, Args) finalResult)
			if((op == "+") || (op == "-"))
		{
			mixin("this = this "~op~" finalResult;");
		}
	}
}
