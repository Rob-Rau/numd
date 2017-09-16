# Numd
A simple @nogc matrix math library for the
[D Programming language](https://dlang.org).

## Intro
Numd is currently just a matrix math library intended to support
[EbbCFD](https://github.com/Rob-Rau/EbbCFD). It can perform basic
matrix math operations in an intuitive API. It is not particularly
fast or efficient but it does get the job done. More work on
performance is intended (quite possibly by utilizing
[mir](https://github.com/libmir/mir))

## Usage
One of the primary goals of Numd's matrix math library was ease of
use. If one is familliar with blas and lapack, the matrix routines
are very cumbersom to use. On the flip side, a package like Matlab
has very intuitive matrix math support. I aimed for the later when
designing this library.

I also wanted to use D's type system to its fullest when performing
matrix operations. This means emmiting compile errors when operating
on incompatible matricies, not allowing certain functions to be called
on matricies that cannot support them (e.g. inverse on non-square
matricies, dot on non-vectors, and cross on non length 3 vectors).
The effect of this is that the matrix struct is a template that takes
in its layout at compile time. This does have limitations, but is
acceptable in many cases.

Further, I wanted to be able to use this library in a @nogc context.
All matrix operations are @nogc and the matrix object uses a statically
allocated array for its underlying data. For a heap allocated matrix
the aliases ```RMatrix``` and ```RVector``` provide easy heap matricies
that use **ALMOST** the same semantics. However, because the
```RefCounted``` wrapper has an ```opAssign```, when you assign one
matrix to another, it only assigns a reference and does not copy the
underlying data.

Some basic usage:
```D
auto m1 = Matrix!(3, 2)(1, 2,
                        3, 4,
                        5, 6);

auto m2 = Matrix!(2, 2)(7, 8,
                        9, 10);
auto m3 = m1 * m2;
writeln(m3);
/+
[  25  28 ]
[  57  64 ]
[  89 100 ]
+/
```

## Future
Originally Numd was to be a complete numerical methods library.
This original vision has not been lost but in order to reduce
dependancies that made it somewhat difficult to get running on
a few clusters, it has been stripped down to only include the
matrix math library for now.

The plan going forward is to move things from EbbCFD into Numd
when they reach a stable level with an interface that I am happy
with. The first things that will likely be brought over are the
time integrators that EbbCFD uses. They currently have a specific
interface to EbbCFD but this should be changing soon.
