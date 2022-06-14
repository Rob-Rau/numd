project('numd', 'd')

sources = [
	'./source/numd/utility.d',
	'./source/numd/calculus/integration/forwardeuler.d',
	'./source/numd/calculus/integration/rk4.d',
	'./source/numd/linearalgebra/matrix.d'
]

numd_incdir = include_directories('source')

if get_option('native')
	add_global_arguments('-mcpu=native', language: 'd')
endif

if get_option('avx512')
	assert(meson.get_compiler('d').get_id() == 'llvm', 'The LLVM D Compiler (LDC) is required for avx512 support')
	add_global_arguments('-mattr=-prefer-256-bit,+prefer-mask-registers', language: 'd')
endif

numd_lib = static_library('numd', sources, include_directories: numd_incdir)

numd_dep = declare_dependency(link_with: numd_lib, include_directories: numd_incdir)
