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

s*A*B

A.opBinaryRight(s).opBinary(B)

*/

unittest
{
	import std.traits : fullyQualifiedName;
	import std.conv : to;
	auto mat1 = Matrix!(2, 2)(0);
	auto mat2 = Matrix!(2, 2)(0);
	Matrix!(2, 2) mat3 = Matrix!(2, 2)(0);
	mat3 = mat2*mat1;
	
	/+
	pragma(msg, fullyQualifiedName!mat1);
	pragma(msg, fullyQualifiedName!mat2);
	pragma(msg, fullyQualifiedName!mat3);
	+/
	
	Matrix!(2, 2) mat4 = mat3*mat2*mat1;
	pragma(msg, fullyQualifiedName!mat4);
	Matrix!(2, 2) mat5 = 4.0*mat4*mat3*mat2*mat1;
	

}
//struct MatFunctor(T, operations...)
struct Matrix(size_t r, size_t c, T = double, operations...)
{
	//T a;

	//Matrix* lhs;
	import std.meta : AliasSeq;
	import std.traits : hasMember;
	enum isFunctor(rhsT) = !hasMember!(T[], "mData");
	alias ops = AliasSeq!operations;
	alias ThisType = Matrix!(r, c, T, operations);
	
	static if(ops.length > 0)
	{
		T[][ops.length+1] datas;
	}
	else
	{
		T[] mData;
	}

	//this(ir, ic, rhsT, opers...)(ref Matrix!(ir, ic, rhsT, opers) funky)
	//	if(funky.ops.length > 0)
	import std.traits : isIntegral, isFloatingPoint;
	this(rhsT)(rhsT funky)
		if(isFunctor!rhsT)
	{
		import std.traits : fullyQualifiedName;
		import std.conv : to;
		pragma(msg, "funky ctor "~fullyQualifiedName!ThisType);
		foreach(int i, op; funky.ops)
		{
			pragma(msg, i.to!string~" "~op);
		}
	}

	this(T init)
	{
		
	}
	/+
	this(ThisType lhs)
	{
		//this.lhs = lhs;
	}
	+/
	this(this)
	{
		import std.traits : fullyQualifiedName;
		pragma(msg, "postblit: "~fullyQualifiedName!ThisType);
	}
	
	ref ThisType opAssign(size_t ir, size_t ic, rhsT, opers...)(Matrix!(ir, ic, rhsT, opers) rhs)
		if((rhs.ops.length > 0) && (ops.length == 0))
	{
		import std.traits : fullyQualifiedName;
		import std.conv : to;
		pragma(msg, to!string(rhs.ops)~": got dem ops");
		pragma(msg, fullyQualifiedName!ThisType~": got dem ops");
		
		foreach(int i, op; rhs.ops)
		{
			pragma(msg, i.to!string~" "~op);
		}
		return this;
	}


/+
	auto ref fmap(alias f)()
	{
		return MatFunctor
	}
+/
	//MatFunctor!(T, ops, op) opBinary(string op)(ref Matrix rhs)
	
	/+
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
	Matrix!(r, ic, rhsType, ops, opers, op) opBinary(string op, size_t ic, rhsType, opers...)(ref Matrix!(r, ic, rhsType, opers) rhs)
		if(op is "*")
	{
		//Matrix C;
		//gemm(a, lhs, rhs, 0, C);
		Matrix!(r, ic, T, ops, opers, op) res;
		
		static if(ops.length > 0)
		{
			pragma(msg, "transfering data ref");
			foreach(int i, ref subData; datas)
			{
				res.datas[i] = subData;
			}
			//res.datas[0..$-1] = datas[];
		}
		else static if((ops.length == 0) && (rhs.ops.length == 0))
		{
			pragma(msg, "transfering this data ref");
			//if(rhs.ops.le)
			res.datas[0] = rhs.mData;
			res.datas[1] = mData;
		}
		else
		{
			static assert(false, "nope");
		}
		
		return res;
	}

	/++
		Scalar multiplication
	+/
	import std.traits : isFloatingPoint;
	Matrix!(r, c, T, "s"~op) opBinaryRight(string op, T)(T s)
		if(isFloatingPoint!T)
	{
		import std.traits : fullyQualifiedName;
		Matrix!(r, c, T, "s"~op) res;
		pragma(msg, "scalar op");
		pragma(msg, fullyQualifiedName!ThisType);
		//Matrix!(T, op) res;
		//res.
		return res;
		
	}
}