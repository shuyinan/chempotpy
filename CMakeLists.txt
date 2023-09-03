cmake_minimum_required(VERSION 3.17.2...3.26)

project(chempotpy
  DESCRIPTION "chempotpy module"
  LANGUAGES C Fortran
  )

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

find_package(BLAS REQUIRED)
message("BLAS found. ${BLAS_LIBRARIES}")

# f2py headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

file(GLOB_RECURSE sources "*_*.f" "*_*.f90")

# Compile each f2py module
foreach(fortran_src_file ${sources})
    get_filename_component(f2py_module_name ${fortran_src_file} NAME_WE)
    get_filename_component(f2py_dir_name ${fortran_src_file} DIRECTORY)
    get_filename_component(system_name ${f2py_dir_name} NAME)
    string(REGEX REPLACE "[.]f(90)?$" ".pyf" f2py_signature_file ${fortran_src_file})
    set(f2py_module_so "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}/${f2py_module_name}.${Python_SOABI}.so")

    add_custom_command(
      OUTPUT "${f2py_module_so}"
      COMMAND ${PYTHON_EXECUTABLE}  -m "numpy.f2py"
                       "${fortran_src_file}"
                       -c "${f2py_signature_file}"
                       ${BLAS_LIBRARIES}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}"
      COMMAND_EXPAND_LISTS
      DEPENDS "${fortran_src_file}"
    )

    add_custom_target(${f2py_module_name} ALL DEPENDS ${f2py_module_so})

    install(PROGRAMS ${f2py_module_so} DESTINATION "${CMAKE_PROJECT_NAME}/${system_name}")

endforeach()

