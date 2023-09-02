cmake_minimum_required(VERSION 3.17.2...3.26)

project(chempotpy
  VERSION 1.0
  DESCRIPTION "chempotpy module"
  LANGUAGES C Fortran
  )

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message("BLAS found.")
    include_directories(${BLAS_INCLUDE_DIRS})
    separate_arguments(BLAS_LINKER_FLAGS_LIST UNIX_COMMAND "${BLAS_LINKER_FLAGS}")

    # Add -l to linker flags that do not start with -l
    set(BLAS_LINKER_FLAGS_LIST_FIXED "")
    foreach(FLAG ${BLAS_LINKER_FLAGS_LIST})
        if(NOT "${FLAG}" MATCHES "^-l")
            set(FLAG "-l${FLAG}")
        endif()
        list(APPEND BLAS_LINKER_FLAGS_LIST_FIXED "${FLAG}")
    endforeach()
    message("BLAS linker flag: ${BLAS_LINKER_FLAGS_LIST_FIXED}")
endif(BLAS_FOUND)

# f2py headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

file(GLOB_RECURSE sources "*.f" "*.f90")

# Compile each f2py module
foreach(fortran_src_file ${sources})
    get_filename_component(f2py_module_name ${fortran_src_file} NAME_WE)
    string(REGEX REPLACE "[.]f(90)?$" ".pyf" f2py_signature_file ${fortran_src_file})
    set(f2py_module_so "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.${Python_SOABI}.so")

    add_custom_command(
      OUTPUT "${f2py_module_so}"
      COMMAND ${PYTHON_EXECUTABLE}  -m "numpy.f2py"
                       "${fortran_src_file}"
                       -c "${f2py_signature_file}"
                       ${BLAS_LINKER_FLAGS_LIST_FIXED}
                       --f90flags="-std=legacy"
                       --f77flags="-std=legacy"
      DEPENDS "${fortran_src_file}"
    )

    add_custom_target(${f2py_module_name} ALL DEPENDS ${f2py_module_so})

    install(PROGRAMS ${f2py_module_so} DESTINATION ${CMAKE_PROJECT_NAME})

endforeach()

