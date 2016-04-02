/*
14:00:46          evenex | rrau this is awesome!! Have you seen Ilya Yaroshenko's ndslice? It would help make a standardized and expressive frontend for all of your code.
14:01:17     WebFreak001 | ketmar: thanks! Gonna save that somewhere
14:01:49          ketmar | WebFreak001: took that from my iv.strex, as i tired of copypasting it to each project. i wonder why std.string doesn't have it
14:02:02            rrau | evenex: I've seen it but haven't looked to in depth at it, but I was thinking the same thing
14:02:49            rrau | evenex: The main problem I'm faced with is that with the current syntax, using opBinary and whatnot, it looks as if it would be problematic to add a BLAS backend
14:03:04          evenex | rrau: can you describe the problem specifically?
14:09:23            rrau | evenex: so a BLAS matrix multiply looks roughly like: gemm(alpha, A, B, beta, C) which is C = alpha*A*B + beta*C
14:09:41          evenex | rrau C is a result point
14:09:43            rrau | evenex: and I'm not sure how to translate that to D syntax
14:09:43          evenex | pointer*
14:10:10            rrau | evenex: it is, but it can also be apart of the equation, if it's like an update
14:10:28          evenex | rrau : you mean update like for iterated functions?
14:11:20            rrau | evenex: yeah. The BLAS spec specifically says the total operation is C = alpha*A*B + beta*C. so If you don't want C to be apart of the result you set beta = 0
14:11:44            rrau | evenex: and it just does C = alpha*A*B
14:12:46          evenex | gotcha. Well, if you want to have A*B where A,B : matrices go to gemm, you can just say obBinary(B) { mat C; gemm(1, &A, &B, 0, &C); return C; }
14:13:35          evenex | then for multiplying by a scalar, opBinaryRight(S) { mat I = mat.identity; mat C; gemm(S, &A, &I, 0, &C); return C; }
14:13:35          evenex | and so on
14:13:44          evenex | so basically you would provide overloads for the different options
14:15:11          evenex | if you wanted to optimize a bit and let s*a*b+beta*c throw everything into gemm at once, you will have to make it lazy... i would just make matrices functorial
14:16:02            rrau | evenex: ok
14:16:22          evenex | so a*M = M.fmap!(x => x*a), let element accesses evaluate the lifted function per-element, and let a subsequent matrix multiplication pull a into the gemm call
14:17:21          evenex | that should keep the number of unnecessary blas calls down while keeping the semantics consistent from the client perspective    

*/

template TestHeader(string test)
{
	pragma(msg, "");
	pragma(msg, "__________________________________________________________________________");
	pragma(msg, test);
	pragma(msg, "");
}

unittest
{
	alias header = TestHeader!("Basic ops test");
	
	auto mat1 = Matrix!(2, 2)(0);
	auto mat2 = Matrix!(2, 2)(0);
	Matrix!(2, 2) mat3 = Matrix!(2, 2)(0);
	
	//mat3 = mat2.opBinary(mat1);
	mat3 = mat2*mat1;

	pragma(msg, "");
	//Matrix!(2, 2) mat4 = mat3.opBinary(mat2).opBinary(mat1);
	Matrix!(2, 2) mat4 = mat3*mat2*mat1;
	
	pragma(msg, "");
	//Matrix!(2, 2) mat5 = mat4.opBinaryRight(4.0).opBinary(mat3).opBinary(mat2).opBinary(mat1);
	Matrix!(2, 2) mat5 = 4.0*mat4*mat3*mat2*mat1;
	
	//Matrix!(2, 2) mat6 = mat4.opBinary(mat3).opBinary(mat2.opBinary(mat1));
	Matrix!(2, 2) mat6 = mat4*mat3 + mat2*mat1;
	
	auto vec1 = Vector!2(0);
	auto vec2 = Vector!2(0);
	Matrix!(2, 2) mat7 = mat1*mat2*vec1 + mat3*vec2; 
}

unittest
{
	alias header = TestHeader!("Function pass test");
	auto mat1 = Matrix!(2, 2)(0);
	Vector!2 vec1;
	Vector!2 vec2 = Vector!2(0);
	Vector!2 vec3 = Vector!2(0);
	import std.conv : to;
	import std.stdio : readln;
	import std.string : strip;
	double s = readln.strip.to!double;

	auto doathing(T)(T vec)
		if(isMatrix!T)
	{
		Matrix!(2, 2) mat;
		// return mat.opBinary(vec);
		return mat*vec;
	}
	//vec1 = doathing(vec3.opBinaryRight(s).opBinary!"+"(mat1.opBinary(vec2)).opBinary!"-"(vec1.opBinaryRight(s)));
	vec1 = doathing(s*vec3 + mat1*vec2 - s*vec1);
	
	pragma(msg, "");
	auto doanotherthing(T)(T vec)
		if(isMatrix!T)
	{
		Matrix!(2, 2) mat;
		Vector!(2) vec1;
		// return mat.opBinary(vec).opBinary!"+"(vec1);
		// return mat.opBinary(vec3.opBinaryRight(s).opBinary!"+"(mat1.opBinary(vec2)).opBinary!"-"(vec1.opBinaryRight(s))).opBinary!"+"(vec1);
		return mat*vec + vec1;
	}
	//vec1 = doathing(vec3.opBinaryRight(s).opBinary!"+"(mat1.opBinary(vec2)).opBinary!"-"(vec1.opBinaryRight(s)));
	vec1 = doanotherthing(s*vec3 + mat1*vec2 - s*vec1);
}

import std.algorithm.comparison : among;
import std.string : indexOf;

enum isMatrix(rhsT) = (("data".among(__traits(allMembers, rhsT)) || "dataStack".among(__traits(allMembers, rhsT))) && (rhsT.stringof.indexOf("Matrix") != -1));
enum isFunctor(rhsT) = isMatrix!rhsT && "dataStack".among(__traits(allMembers, rhsT));

private auto filterScalarOps(ops...)()
{
	import std.algorithm.iteration : filter;
	import std.range : array;
	return [ops].filter!(a => (a.indexOf('s') != -1)).array;
}

private auto filterMatrixOps(ops...)()
{
	import std.algorithm.iteration : filter;
	import std.range : array;
	return [ops].filter!(a => (a.indexOf('s') == -1)).array;
}
	
alias Vector(size_t r, T = double) = Matrix!(r, 1, T);
struct Matrix(size_t r, size_t c, T = double, operations...)
{
@nogc:

	import std.meta : AliasSeq, aliasSeqOf;
	alias ops = AliasSeq!operations;
	alias ThisType = Matrix!(r, c, T, operations);

	static if(ops.length > 0)
	{
		immutable string[] scalarOps = filterScalarOps!(ops);
		immutable string[] matrixOps = filterMatrixOps!(ops);
		T[][matrixOps.length+1] dataStack;
		T[scalarOps.length] scalars;
		//pragma(msg, scalarOps);
		//pragma(msg, matrixOps);
		//pragma(msg, [ops]);
	}
	else
	{
		// destination storage
		T[] data;
		
		// temporary scratch data only ever alloced when needed.
		// This is for long expressions that have intermediate
		// steps that need more space than the destination has.
		static T[] scratchData;
		
		static ~this()
		{
			if(scratchData !is null)
			{
				// TODO: dealloc
			}
		}
	}

	import std.traits : isIntegral, isFloatingPoint;
	this(rhsT)(rhsT matFunctor)
		if(!isIntegral!rhsT && !isFloatingPoint!rhsT && isFunctor!rhsT)
	{
		import std.conv : to;
		pragma(msg, "functor ctor "~ThisType.stringof);

		this = opAssign(matFunctor);
	}

	this(T init)
	{
		
	}
	
	this(this)
	{
		//pragma(msg, "postblit: "~ThisType.stringof);
	}

	//ref ThisType opAssign(size_t ir, size_t ic, rhsT, opers...)(Matrix!(ir, ic, rhsT, opers) rhs)
	//	if((rhs.ops.length > 0) && (ops.length == 0))
	ref ThisType opAssign(rhsT)(rhsT rhs)
		if(isFunctor!rhsT && !isFunctor!ThisType)
	{
		import std.conv : to;
		//pragma(msg, to!string(rhs.ops)~": got dem ops");
		pragma(msg, "opAssign "~ThisType.stringof);
		//pragma(msg, isFunctor!rhsT);
		//pragma(msg, rhsT.stringof);
		//pragma(msg, isFunctor!ThisType);
		//pragma(msg, ThisType.stringof);
		//pragma(msg, rhs.datas.length);
		
		//immutable int sidx0 = 0;
		//foreach(int i, op; aliasSeqOf!(rhs.matrixOps))
		foreach(int i, op; rhs.ops)
		{
			static if(op.indexOf('s') != -1)
			{
				
				import std.string : chompPrefix;
				pragma(msg, i.to!string~" scalar op: "~op.chompPrefix("s"));
			}
			else
			{
				pragma(msg, i.to!string~" "~op);
			}
		}
		return this;
	}
	
	/++
		A*B*C
		A.opBinary(B).opBinary(C);
	
		A*B*C*D
		tmp1 = A.opBinary(B)
		tmp2 = tmp1.opBinary(C)
		ans = tmp2.opBinary(D);
		A.opBinary(B).opBinary(C).opBinary(D)
		
		s*A*B*C*D
		A.opBinaryRight(s).opBinary(B).opBinary(C).opBinary(D)
		
		A*B + C*D
		A.opBinary(B).opBinary!"+"(C.opBinary(D))
	+/

	//Matrix!(r, ic, rhsType, ops, opers, op) opBinary(string op, size_t ic, rhsType, opers...)
	auto opBinary(string op, size_t ic, rhsType, opers...)(Matrix!(r, ic, rhsType, opers) rhs)
		if((op is "*") || (op == "+") || (op is "-"))
	{
		import std.format : format;
		mixin(format("Matrix!(r, ic, T, ops, %s) res;", (op == "+") || (op == "-") ? "op, opers" : "opers, op"));
		
		static if(ops.length > 0)
		{
			//pragma(msg, "transfering data ref");
			foreach(int i, ref subData; dataStack)
			{
				res.dataStack[i] = subData;
			}
			//res.datas[0..$-1] = datas[];
		}
		else static if((ops.length == 0) && (rhs.ops.length == 0))
		{
			//pragma(msg, "transfering this data ref");
			//if(rhs.ops.le)
			res.dataStack[0] = rhs.data;
			res.dataStack[1] = data;
		}
		else static if((ops.length == 0) && (rhs.ops.length > 0))
		{
			//static assert(false, "nope");
		}
		
		return res;
	}

	/++
		Scalar multiplication
	+/
	import std.traits : isFloatingPoint;
	auto opBinaryRight(string op, T)(T s)
		if(isFloatingPoint!T && (op == "*"))
	{
		Matrix!(r, c, T, ops, "s"~op) res;
		res.scalars[$-1] = s;
		//pragma(msg, "scalar op");
		//pragma(msg, ThisType.stringof);
		//res.datas[$-1] = [s];
		//Matrix!(T, op) res;
		//res.
		return res;
		
	}
}