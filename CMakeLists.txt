cmake_minimum_required(VERSION 3.0)
project(numd LANGUAGES D)

include_directories(source)
add_library(numd
        source/numd/linearalgebra/matrix.d
        source/numd/utility.d
        )
target_link_libraries(numd)

#set(view_dir ${CMAKE_CURRENT_SOURCE_DIR}/views)
#target_compile_options(matplotlib-d PRIVATE "-J ${view_dir}")
#set_target_properties(matplotlib-d PROPERTIES TEXT_INCLUDE_DIRECTORIES view_dir)
