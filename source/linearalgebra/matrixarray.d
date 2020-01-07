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
	//immutable size_t arraySize = 10_000_000;
	immutable size_t arraySize = 16;
	//immutable size_t arraySize = 10000;
	immutable size_t vecSize = 4;
	//auto n0 = new double[arraySize];
	auto n0 = ScalarArray!double(arraySize, 0.5);
	//n0[] = 0.5;
	//auto n1 = new double[arraySize];
	auto n1 = ScalarArray!double(arraySize, 0.5);
	//n1[] = 0.5;
	//auto Q = VectorArray!vecSize(arraySize, 1.0);
	auto dqdx = VectorArray!vecSize(arraySize, 0.1);
	auto dqdy = VectorArray!vecSize(arraySize, 1.0);
	//auto A00 = MatrixArray!(vecSize, vecSize)(arraySize, 2.0);
	//auto A01 = MatrixArray!(vecSize, vecSize)(arraySize, 3.0);
	//auto A10 = MatrixArray!(vecSize, vecSize)(arraySize, 4.0);
	//auto A11 = MatrixArray!(vecSize, vecSize)(arraySize, 5.0);
	auto flux = VectorArray!vecSize(arraySize, 5.0);

	/+
	zip(Q, gradMat, indicies).map!((a) {
		//auto dq = a[2].map!()
	});
	+/

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

		
		/+immutable size_t chunkSize = 90;
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

			auto n0_slice = n0.slice(i, sliceEnd);
			auto n1_slice = n1.slice(i, sliceEnd);
			flux[i..sliceEnd] = n0_slice*G0 + n1_slice*G1 + n0_slice*G1 + n1_slice*G0 + n0_slice*(G2 + G3) - n1_slice*(G2 - G3);
		}+/
		
		//pragma(msg, "1");
		auto A00 = dqdy*dqdx.transpose;
		//flux = A00*dqdy;
		//pragma(msg, typeof(A00));
		auto A01 = dqdx*dqdy.transpose;
		auto A10 = dqdy*dqdx.transpose;
		auto A11 = dqdx*dqdy.transpose;

		//pragma(msg, "2");
		auto G0 = A00*dqdx + A01*dqdy;
		//flux = G0;
		auto G1 = A10*dqdx + A11*dqdy;
		//flux = n0*G0 + n1*G1;
		//flux = n0*G0 + n1*G1 + n0*G1 + n1*G0 + n0*(G0 + G1);
		auto G2 = A00*dqdy + A01*dqdx;
		auto G3 = A10*dqdy + A11*dqdx;
		//pragma(msg, "3");
		flux = n0*G0 + n1*G1 + n0*G1 + n1*G0 + n0*(G2 + G3) - n1*(G2 - G3);
		//pragma(msg, "4");
		
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

			fluxo[i] = n0[i]*G0 + n1[i]*G1 + n0 [i]*G1 + n1[i]*G0 + n0[i]*(G2 + G3) - n1[i]*(G2 - G3);
			//fluxo[i] = n0[i]*G0 + n1[i]*G1;
			//pragma(msg, typeof(A00));
			//fluxo[i] = A00*dqdyo[i];
		}
	}

	writeln("Running diffuse flux test");
	import std.datetime.stopwatch : benchmark;
	//uint compTimes = 5_000_000;
	uint compTimes = 40_000_000;
	//auto r = benchmark!nogcComputeAdd(compTimes);
	//writeln("Addition computation time ", r[0]);

	writeln;
	writeln("Running old diffuse flux test");
	auto r = benchmark!diffuseFluxOld(compTimes);
	writeln("Old diffuse flux time ", r[0]);

	writeln;
	writeln("Running new diffuse flux test");
	r = benchmark!diffuseFlux(compTimes);
	writeln("New diffuse flux time ", r[0]);

	

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
immutable ulong dataCacheBytes = 32*1024;
immutable ulong dataCacheBytes40percent = cast(ulong)(dataCacheBytes*0.4);
unittest
{
	import std.stdio : writeln;

	writeln;
	writeln("=======================");
	writeln("Allocating big arrays");
	size_t arraySize = 100_000_00;
	auto ma1 = MatrixArray!(4,4)(arraySize, 2);
	auto ma2 = MatrixArray!(4,4)(arraySize, 3);
	auto ma4 = MatrixArray!(4,4)(arraySize, 4);

	immutable ulong dataSizeBytes = 3*4*4*8*arraySize;
	//immutable ulong dataSizeInCache = dataSizeBytes/dataCacheBytes40percent;
	immutable ulong chunkSize = dataCacheBytes40percent/(3*4*4*8);
	immutable ulong remainderChunk = arraySize%chunkSize;
	writeln("40% Cache size: ", dataCacheBytes40percent);
	writeln("Data size: ", dataSizeBytes);
	//writeln("Data in 40% cache: ", dataSizeInCache);
	writeln("Chunk size: ", chunkSize);
	writeln("Remainder chunk: ", remainderChunk);


	//auto ma5 = MatrixArray!(4,4)(arraySize);
	
	/+@nogc void nogcComputeAdd()
	{
		auto ma3 = ma1 + ma2;
		ma4 = ma3;
		ma5 = ma1 + ma2 + ma4;
	}+/

	@nogc void nogcComputeMult()
	{
		//ma4 = ma1*ma2;
		for(ulong i = 0; i < ma4.length; i+= chunkSize)
		{
			ma4[i..i+chunkSize] = ma1[i..i+chunkSize]*ma2[i..i+chunkSize];
		}

		if(remainderChunk > 0)
		{
			ma4[arraySize - remainderChunk..arraySize] = ma1[arraySize - remainderChunk..arraySize]*ma2[arraySize - remainderChunk..arraySize];
		}
	}

	import std.datetime.stopwatch : benchmark;
	immutable uint compTimes = 200;
	writeln("Running matrix array test");
	/+
	auto r = benchmark!nogcComputeAdd(compTimes);
	writeln("Addition computation time ", r[0]);
	+/

	auto r = benchmark!nogcComputeMult(compTimes);
	writeln("Multiplication computation time ", r[0]/compTimes, " per iteration");

}+/

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

// Builds a named parameter tuple based on the args and typeid of the args provided.
string buildArgTuple(Args...)() {
	string tupleDecl = "Tuple!(";
	foreach(arg; AliasSeq!Args) {
		tupleDecl ~= arg.stringof ~ `, "arg` ~ arg.id ~ `",`;
	}
	tupleDecl = tupleDecl[0..$-1]; // remove trailing comma
	tupleDecl ~= ") args;";
	return tupleDecl;
}

string buildArgTupleConstructor(Args...)() {
	string tupleDecl = "tuple!(";
	foreach(arg; AliasSeq!Args) {
		tupleDecl ~= `"arg` ~ arg.id ~ `",`;
	}
	tupleDecl = tupleDecl[0..$-1]; // remove trailing comma
	tupleDecl ~= ")";
	return tupleDecl;
}

string buildArgTupleConstructor2(string code, Args...)() {
	string tupleDecl = "tuple!(";
	foreach(arg; AliasSeq!Args) {
		tupleDecl ~= `"arg` ~ arg.id ~ `",`;
	}
	tupleDecl = tupleDecl[0..$-1]; // remove trailing comma
	tupleDecl ~= ")("~code~")";
	return tupleDecl;
}

struct ScalarArray(T = double, uint line = __LINE__, string file = __FILE__) {
	
	alias ThisType = typeof(this);

	enum string id = to!string(ThisType.mangleof.hashOf);
	//pragma(msg, "array id (str): "~ThisType.mangleof);
	//pragma(msg, "array id: "~id.to!string);

	T[] data;
	alias data this;
	
	@trusted this(size_t arraySize) {
		data = new T[arraySize];
	}

	@trusted this(size_t arraySize, T initialValue) {
		data = new T[arraySize];
		data[] = initialValue;
	}

	@nogc {
		@trusted private this(T[] _data) {
			data = _data;
		}

		/+@trusted auto opBinary(string op, uint line, string file)(ref inout ScalarArray!(T, line, file) inRhs)
			if(op == "+" || op == "-")
		{

			alias InType = ScalarArray!(T, line, file);
			assert(this.length == inRhs.length, "Scalar arrays are not the same size");
			
			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `return "T[`
					~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.data[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.data[`~range~`];\n";`;
				mixin(code);
			}
			return MatrixArrayResult!(r, c, T, computationFunc, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}+/

		/+@trusted auto opBinary(string op, uint line, string file)(ScalarArray!(T, line, file) inRhs)
			if(op == "+" || op == "-")
		{

			alias InType = ScalarArray!(T, line, file);
			assert(this.length == inRhs.length, "Scalar arrays are not the same size");
			
			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `return "T[`
					~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.data[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.data[`~range~`];\n";`;
				mixin(code);
			}
			return MatrixArrayResult!(r, c, T, computationFunc, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}+/
	}
}


private struct MatrixArrayResult(size_t r, size_t c, T = double, alias operation, LhsType, RhsType, Args...)
{
	alias ThisType = typeof(this);

	enum string id = to!string(ThisType.mangleof.hashOf);
	alias _lhsId = RhsType;
	alias _rhsId = LhsType;

	//pragma(msg, "array result id: "~id.to!string);

	mixin(buildArgTuple!Args);
	
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

	@trusted @nogc auto opBinary(string op, size_t ic, uint line, string file)(ref inout MatrixArray!(c, ic, T, line, file) inRhs)
		if(op == "*")
	{
		alias InType = MatrixArray!(c, ic, T, line, file);
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){

			static foreach(a; aliasSeqOf!(iota(0, c)))
			{
				mixin(`enum string code`~a.to!string~` = "`~operation!(ThisType.id, i, to!string(a), range, destSize)~`";`);
			}
			mixin(`enum string code_required = code0`~iota(1, c).map!(a => "~code"~a.to!string).join~`;`);
			enum string code = `return "`
				~ code_required
				~ `T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~`0[] * finalResult.args.arg`~InType.id~`.q0`~j~`[`~range~`] `~iota(1, c).map!(a => " + temp_"~ThisType.id~"_"~i~a.to!string~"[]"~op~"finalResult.args.arg"~InType.id~".q"~a.to!string~j~"["~range~"]").join~`;\n";`;
			mixin(code);
		}

		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input matrix array is already in our argument list
			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inRhs));
		}
	}

	@trusted @nogc auto opBinary(string op, size_t ic, uint line, string file)(inout MatrixArray!(c, ic, T, line, file) inRhs)
		if(op == "*")
	{
		alias InType = MatrixArray!(c, ic, T, line, file);
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){

			static foreach(a; aliasSeqOf!(iota(0, c)))
			{
				mixin(`enum string code`~a.to!string~` = "`~operation!(ThisType.id, i, to!string(a), range, destSize)~`";`);
			}
			mixin(`enum string code_required = code0`~iota(1, c).map!(a => "~code"~a.to!string).join~`;`);
			enum string code = code_required ~ `return "T[`
				~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~`0[] * finalResult.args.arg`~InType.id~`.q0`~j~`[`~range~`] `~iota(1, c).map!(a => " + temp_"~ThisType.id~"_"~i~a.to!string~"[]"~op~"finalResult.args.arg"~InType.id~".q"~a.to!string~j~"["~range~"]").join~`;\n";`;
			mixin(code);
		}

		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input matrix array is already in our argument list
			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inRhs));
		}
	}
	
	@trusted @nogc auto opBinary(string op, uint line, string file)(ref inout MatrixArray!(r, c, T, line, file) inRhs)
		if((op == "+") || (op == "-"))
	{
		alias InType = MatrixArray!(r, c, T, line, file);
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			enum string code = `return "`
					~operation!(ThisType.id, i, j, range, destSize)~
					`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~j~`[] `~op~` finalResult.args.arg`~InType.id~`.q`~i~j~`[`~range~`];\n";`;
			mixin(code);
		}

		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input matrix array is already in our argument list
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inRhs));
		}
	}

	@trusted @nogc auto opBinary(string op, uint line, string file)(inout MatrixArray!(r, c, T, line, file) inRhs)
		if((op == "+") || (op == "-"))
	{
		alias InType = MatrixArray!(r, c, T, line, file);
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			enum string code = `return "`
					~operation!(ThisType.id, i, j, range, destSize)~
					`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~j~`[] `~op~` finalResult.args.arg`~InType.id~`.q`~i~j~`[`~range~`];\n";`;
			mixin(code);
		}

		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input matrix array is already in our argument list
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inRhs));
		}
	}

	@trusted @nogc auto opBinary(string op, alias inOperation, InLhsType, InRhsType, inArgs...)(MatrixArrayResult!(r, c, T, inOperation, InLhsType, InRhsType, inArgs) inRhs)
		if((op == "+") || (op == "-"))
	{
		alias InType = MatrixArrayResult!(r, c, T, inOperation, LhsType, RhsType, inArgs);
		assert(this.length == inRhs.length, "Matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
			enum string code = `return "`
					~operation!(ThisType.id, i, j, range, destSize)
					~inOperation!(InType.id, i, j, range, destSize)~
					`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~j~`[]`~op~`temp_`~InType.id~`_`~i~j~`[];\n";`;
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			mixin(code);
		}


		//return MatrixArrayResult!(r, c, T, comp, ThisType, InType, NoDuplicates!(AliasSeq!(Args, inArgs)))(tuple(this.args.expand, inRhs.args.expand));
		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, inArgs));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input matrix array is already in our argument list
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			//mixin(buildArgTuple!(FilteredArgs));
			static foreach(i, filteredArg; FilteredArgs)
			{
				static if(__traits(compiles, mixin("this.args.arg"~filteredArg.id)))
				{
					static if(i == FilteredArgs.length - 1)
					{
						mixin("enum string code"~to!string(i)~` = "this.args.arg"~filteredArg.id;`);
					}
					else
					{
						mixin("enum string code"~to!string(i)~` = "this.args.arg"~filteredArg.id~",";`);
					}
					
				}
				else
				{
					static if(i == FilteredArgs.length - 1)
					{
						mixin("enum string code"~to!string(i)~` = "inRhs.args.arg"~filteredArg.id;`);
					}
					else
					{
						mixin("enum string code"~to!string(i)~` = "inRhs.args.arg"~filteredArg.id~",";`);
					}
					//mixin("args.arg"~filteredArg.id~" = inRhs.args.arg"~filteredArg.id~";");
					
				}
			}
			mixin("enum string code = code0 "~iota(1,FilteredArgs.length).map!(a => "~ code"~to!string(a)).join~";");
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, FilteredArgs)(mixin(buildArgTupleConstructor2!(code, FilteredArgs)));
		}
	}
	
	//@trusted @nogc auto opBinaryRight(string op)(inout T[] inLhs)
	@trusted @nogc auto opBinaryRight(string op, uint line, string file)(ScalarArray!(T, line, file) inLhs)
		if(op == "*")
	{
		alias InType = ScalarArray!(T, line, file);
		assert(this.length == inLhs.length, "Scalar and matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			enum string code = `return "`
					~operation!(ThisType.id, i, j, range, destSize)~
					`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~j~`[] `~op~` finalResult.args.arg`~InType.id~`[`~range~`];\n";`;
			mixin(code);
		}

		//return MatrixArrayResult!(r, c, T, comp, InType, ThisType, Args, InType)(tuple(this.args.expand, cast(T[])inLhs));
		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input scalar array is already in our argument list
			return MatrixArrayResult!(r, c, T, comp, InType, ThisType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, c, T, comp, InType, ThisType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inLhs));
		}
	}

	@trusted @nogc auto opBinaryRight(string op, uint line, string file)(ref ScalarArray!(T, line, file) inLhs)
		if(op == "*")
	{
		alias InType = ScalarArray!(T, line, file);
		assert(this.length == inLhs.length, "Scalar and matrix arrays are not the same size");

		@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
			enum string code = `return "`
					~operation!(ThisType.id, i, j, range, destSize)~
					`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = temp_`~ThisType.id~`_`~i~j~`[] `~op~` finalResult.args.arg`~InType.id~`[`~range~`];\n";`;
			mixin(code);
		}

		//return MatrixArrayResult!(r, c, T, comp, InType, ThisType, Args, InType)(tuple(this.args.expand, cast(T[])inLhs));
		alias FilteredArgs = NoDuplicates!(AliasSeq!(Args, InType));
		static if(FilteredArgs.length == AliasSeq!(Args).length)
		{
			// Input scalar array is already in our argument list
			return MatrixArrayResult!(r, c, T, comp, InType, ThisType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand));
		}
		else
		{
			return MatrixArrayResult!(r, c, T, comp, InType, ThisType, FilteredArgs)(mixin(buildArgTupleConstructor!(FilteredArgs))(this.args.expand, cast(InType)inLhs));
		}
	}

	@trusted @nogc auto opBinaryRight(string op)(ref inout T inLhs)
		if(op == "*")
	{
		@nogc static immutable(string) comp(string cop, uint level, string i, string j, string[] args, string range, string destSize)(){
			enum string code = `return "`~operation!(cop, level, i, j, args[0..this.args.length], range, destSize)~
				`temp`~level.to!string~`[] `~op~`= finalResult.args[`~args[$-1]~`];\n";`;
			//debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
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

alias VectorArray(size_t l, T = double, uint line = __LINE__, string file = __FILE__) = MatrixArray!(l, 1, T, line, file);

private struct EnumResult
{
	this(ulong _line, string _code) {
		line = _line;
		code = _code;
	}
	ulong line;
	string code;
}

static EnumResult[] enumerateCode(string[] code) {
	return zip(iota(0, code.length), code)
		.map!(a => EnumResult(a[0], a[1]))
		.array;
}

struct MatrixArray(size_t r, size_t c, T = double, uint line = __LINE__, string file = __FILE__)
{

	alias ThisType = typeof(this);

	//enum ulong id = ThisType.mangleof.hashOf;
	enum string id = to!string(ThisType.mangleof.hashOf);
	enum size_t rows = r;
	enum size_t columns = c;
	
	//pragma(msg, "array id (str): "~ThisType.mangleof);
	//pragma(msg, "array id: "~id.to!string);

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

		@trusted MatrixArray!(c, r, T, line, file) transpose()
		{
			auto ret = MatrixArray!(c, r, T, line, file)(this.length, false);
			static foreach(i; aliasSeqOf!(iota(0, c)))
			{
				static foreach(j; aliasSeqOf!(iota(0, r)))
				{
					mixin("ret.q"~i.to!string~j.to!string~" = this.q"~j.to!string~i.to!string~";");
				}
			}
			return ret;
		}

		@trusted auto scalarArray(size_t i, size_t j = 0)()
		{
			auto ret = MatrixArray!(1, 1, T, line, file~i.to!string~j.to!string)(this.length, false);
			//return MatrixArray!(1, 1, T, line, file~i.to!string~j.to!string)(mixin("this.q"~i.to!string~j.to!string));
			mixin("ret.q00 = this.q"~i.to!string~j.to!string~";");
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

		@trusted ref ThisType opAssign(alias operation, LhsType, RhsType, Args...)(MatrixArrayResult!(r, c, T, operation, LhsType, RhsType, Args) finalResult)
		{

			alias InType = MatrixArrayResult!(r, c, T, operation, LhsType, RhsType, Args);
			import core.stdc.stdio : printf;
			enum size_t stride16 = 16;
			enum size_t stride8 = 8;
			enum size_t stride4 = 4;
			enum size_t stride2 = 2;
			size_t k;
			size_t iterations;
			size_t remaining;
			pragma(inline, true);

			static foreach(i; aliasSeqOf!(iota(0, r)))
			{
				static foreach(j; aliasSeqOf!(iota(0, c)))
				{
					k = 0;
					iterations = finalResult.length/stride16;
					remaining = finalResult.length%stride16;

					/++
						There are 2 reasons why we stride like this:
							1) So we can have stack apporpriate sized intermediate calculations
							2) This forces the compiler to properly vectorize. If we do single element by single element,
							   ldc2 doesn't want to vectorize (although perhaps I need to explicitly pass in force-unroll-loops)
							   (checked against ldc2 1.18 and earlier)
					+/
					iterations = finalResult.length/stride8;
					remaining = finalResult.length%stride8;
					import std.string : strip;
					while(iterations)
					{
						enum string code = operation!(
								InType.id,
								to!string(i),
								to!string(j),
								"k..k+"~to!string(stride8),
								to!string(stride8)
							)
							.split("\n")
							.enumerateCode
							.array
							.multiSort!("a.code.hashOf < b.code.hashOf", "a.line < b.line")
							.uniq!((a, b) => a.code == b.code)
							.array
							.sort!((a, b) => a.line < b.line)
							.map!(a => a.code)
							.join("\n")
							;
					
						//pragma(msg, code);
						mixin(code);
						mixin("this.q"~i.to!string~j.to!string~"[k..k+stride8] = temp_"~InType.id.to!string~"_"~i.to!string~j.to!string~"[];");
						iterations--;
						k+=stride8;
					}

					if(remaining > stride4)
					{
						iterations = remaining/stride4;
						remaining = remaining%stride4;
						while(iterations)
						{
							enum string code =
								operation!(
									InType.id,
									to!string(i),
									to!string(j),
									"k..k+"~to!string(stride4),
									to!string(stride4)
								)
								.split("\n")
								.enumerateCode
								.array
								.multiSort!("a.code.hashOf < b.code.hashOf", "a.line < b.line")
								.uniq!((a, b) => a.code == b.code)
								.array
								.sort!((a, b) => a.line < b.line)
								.map!(a => a.code)
								.join("\n")
								;
							mixin(code);
							mixin("this.q"~i.to!string~j.to!string~"[k..k+stride4] = temp_"~InType.id.to!string~"_"~i.to!string~j.to!string~"[];");
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
							enum string code =
								operation!(
									InType.id,
									to!string(i),
									to!string(j),
									"k..k+"~to!string(stride2),
									to!string(stride2)
								)
								.split("\n")
								.enumerateCode
								.array
								.multiSort!("a.code.hashOf < b.code.hashOf", "a.line < b.line")
								.uniq!((a, b) => a.code == b.code)
								.array
								.sort!((a, b) => a.line < b.line)
								.map!(a => a.code)
								.join("\n")
								;
							mixin(code);
							mixin("this.q"~i.to!string~j.to!string~"[k..k+stride2] = temp_"~InType.id.to!string~"_"~i.to!string~j.to!string~"[];");
							iterations--;
							k+=stride2;
						}
					}
					
					while(iterations)
					{
						enum string code =
							operation!(
								InType.id,
								to!string(i),
								to!string(j),
								"k..k+1",
								"1"
							)
							.split("\n")
							.enumerateCode
							.array
							.multiSort!("a.code.hashOf < b.code.hashOf", "a.line < b.line")
							.uniq!((a, b) => a.code == b.code)
							.array
							.sort!((a, b) => a.line < b.line)
							.map!(a => a.code)
							.join("\n")
							;
						mixin(code);
						mixin("this.q"~i.to!string~j.to!string~"[k..k+1] = temp_"~InType.id.to!string~"_"~i.to!string~j.to!string~"[];");
						remaining--;
						k++;
					}
				}
			}
			return this;
		}
		
		@trusted auto opBinary(string op, uint line, string file)(ref inout MatrixArray!(r, c, T, line, file) inRhs)
			if(op == "+" || op == "-")
		{

			alias InType = MatrixArray!(r, c, T, line, file);
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			
			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `return "T[`
					~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.q`~i~j~`[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.q`~i~j~`[`~range~`];\n";`;
				mixin(code);
			}
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}

		@trusted auto opBinary(string op, uint line, string file)(inout MatrixArray!(r, c, T, line, file) inRhs)
			if(op == "+" || op == "-")
		{

			alias InType = MatrixArray!(r, c, T, line, file);
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");
			
			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `return "T[`
					~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.q`~i~j~`[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.q`~i~j~`[`~range~`];\n";`;
				mixin(code);
			}
			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}

		@trusted auto opBinary(string op, alias operation, LhsType, RhsType, Args...)(MatrixArrayResult!(r, c, T, operation, LhsType, RhsType, Args) inRhs)
			if(op == "+" || op == "-")
		{
			alias InType = MatrixArrayResult!(r, c, T, operation, LhsType, RhsType, Args);
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");

			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `return "`
						~operation!(InType.id, i, j, range, destSize)~
						`T[`~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.q`~i~j~`[`~range~`] `~op~` temp_`~InType.id~`_`~i~j~`[];\n";`;
				mixin(code);
			}


			return MatrixArrayResult!(r, c, T, comp, ThisType, InType, NoDuplicates!(AliasSeq!(Args, ThisType)))(tuple(inRhs.args.expand, this));
		}

		@trusted auto opBinary(string op, size_t ic, uint line, string file)(ref inout MatrixArray!(r, ic, T, line, file) inRhs)
			if(op == "*")
		{
			alias InType = MatrixArray!(r, ic, T, line, file);
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");

			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `
					return "T[`~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.q`~i~`0[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.q0`~j~`[`~range~`] `~iota(1, c).map!(
							a => " + finalResult.args.arg"~ThisType.id~".q"~i~to!string(a)~"["~range~"] "~op~" finalResult.args.arg"~InType.id~".q"~to!string(a)~j~"["~range~"]"
						).join~`;\n";`;

				mixin(code);
			}

			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}

		@trusted auto opBinary(string op, size_t ic, uint line, string file)(inout MatrixArray!(c, ic, T, line, file) inRhs)
			if(op == "*")
		{
			alias InType = MatrixArray!(c, ic, T, line, file);
			assert(this.length == inRhs.length, "Matrix arrays are not the same size");

			@nogc static immutable(string) comp(string resultId, string i, string j, string range, string destSize)(){
				enum string code = `
					return "T[`~destSize~`] temp_`~resultId~`_`~i~j~` = finalResult.args.arg`~ThisType.id~`.q`~i~`0[`~range~`] `~op~` finalResult.args.arg`~InType.id~`.q0`~j~`[`~range~`]`~iota(1, c).map!(
							a => " + finalResult.args.arg"~ThisType.id~".q"~i~to!string(a)~"["~range~"] "~op~" finalResult.args.arg"~InType.id~".q"~to!string(a)~j~"["~range~"]"
						).join~`;\n";`;

				debug pragma(msg, __LINE__.to!string~": level = "~level.to!string);
				//debug pragma(msg, code);
				mixin(code);
			}

			return MatrixArrayResult!(r, ic, T, comp, ThisType, InType, ThisType, InType)(tuple(this, cast(InType)inRhs));
		}

		/+@trusted auto opBinaryRight(string op)(inout T inLhs)
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
		}+/
	}
}
