cmake_minimum_required(VERSION 3.0)
project(numd LANGUAGES D)

include_directories(source)
add_library(numd
        source/numd/linearalgebra/matrix.d
        source/numd/utility.d
        )
target_link_libraries(numd)

set_target_properties(
    numd PROPERTIES
    #RUNTIME_OUTPUT_DIRECTORY ./
    #LIBRARY_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_NAME numd
    #LIBRARY_OUTPUT_NAME numd
    #RUNTIME_OUTPUT_NAME numd
    #COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    #LINK_FLAGS ""
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)

#set(view_dir ${CMAKE_CURRENT_SOURCE_DIR}/views)
#target_compile_options(matplotlib-d PRIVATE "-J ${view_dir}")
#set_target_properties(matplotlib-d PROPERTIES TEXT_INCLUDE_DIRECTORIES view_dir)
