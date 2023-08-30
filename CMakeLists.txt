cmake_minimum_required(VERSION 3.17.2...3.26)

project(addpkg
  VERSION 1.0
  DESCRIPTION "addpkg module"
  LANGUAGES C Fortran
  )

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

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
    set(f2py_module_c "${f2py_module_name}module.c")

    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
      COMMAND ${PYTHON_EXECUTABLE}  -m "numpy.f2py"
                       "${fortran_src_file}"
                       "${f2py_signature_file}"
                       -m "${f2py_module_name}"
                       --lower
      DEPENDS "${fortran_src_file}"
    )

    python_add_library(${f2py_module_name} MODULE
                "${f2py_module_c}"
                "${F2PY_INCLUDE_DIR}/fortranobject.c"
                "${fortran_src_file}" WITH_SOABI)

    target_include_directories(${f2py_module_name} PUBLIC
                               ${F2PY_INCLUDE_DIR}
                              )

    target_link_libraries(${f2py_module_name} PUBLIC Python::NumPy)


    install(TARGETS ${f2py_module_name} DESTINATION ${CMAKE_PROJECT_NAME})

endforeach()

