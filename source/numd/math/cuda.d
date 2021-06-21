@compute(CompileFor.deviceOnly) module numd.math.cuda;

import ldc.dcompute;

pure: nothrow: @nogc:

pragma(LDC_intrinsic, "llvm.nvvm.fabs.f") float fabs_f(float);
pragma(LDC_intrinsic, "llvm.nvvm.fabs.d") double fabs_d(double);

pragma(LDC_intrinsic, "llvm.nvvm.sqrt.f") float sqrt_f(float);
pragma(LDC_intrinsic, "llvm.nvvm.sqrt.rn.d") double sqrt_d(double);

pragma(LDC_intrinsic, "llvm.nvvm.round.f") float round_f(float);
pragma(LDC_intrinsic, "llvm.nvvm.round.d") double round_d(double);