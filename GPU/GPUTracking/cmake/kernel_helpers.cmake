# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# file kernel_helpers.cmake
# author David Rohr

add_custom_target(O2_GPU_KERNELS)
define_property(TARGET PROPERTY O2_GPU_KERNELS)
define_property(TARGET PROPERTY O2_GPU_KERNEL_NAMES)
define_property(TARGET PROPERTY O2_GPU_KERNEL_INCLUDES)
define_property(TARGET PROPERTY O2_GPU_KERNEL_FILES)
set(O2_GPU_KERNEL_WRAPPER_FOLDER "${CMAKE_CURRENT_BINARY_DIR}/GPU/include_gpu_onthefly")
file(MAKE_DIRECTORY ${O2_GPU_KERNEL_WRAPPER_FOLDER})
set(O2_GPU_BASE_DIR "${CMAKE_CURRENT_LIST_DIR}/../")
function(o2_gpu_add_kernel kernel_name kernel_files kernel_bounds kernel_type)
  math(EXPR TMP_CHK "${ARGC} & 1")
  if(${TMP_CHK})
    message(FATAL_ERROR "Invalid number of arguments to kernel ${TMP_CHK}, must be odd to have pairs of argument type, argument name")
  endif()
  list(LENGTH ARGV n)
  set(OPT1 "")
  set(OPT2 "")
  set(OPT3 "")
  if(${n} GREATER 4)
    math(EXPR n "${n} - 1")
    foreach(i RANGE 4 ${n} 2)
      math(EXPR j "${i} + 1")
      if(${ARGV${i}} MATCHES "\\*$")
        string(APPEND OPT1 ",GPUPtr1(${ARGV${i}},${ARGV${j}})")
        string(APPEND OPT2 ",GPUPtr2(${ARGV${i}},${ARGV${j}})")
      else()
        string(APPEND OPT1 ",${ARGV${i}} ${ARGV${j}}")
        string(APPEND OPT2 ",${ARGV${j}}")
      endif()
      string(APPEND OPT3 ",${ARGV${i}}")
    endforeach()
  endif()
  if(kernel_bounds MATCHES "^LB")
    set(TMP_BOUNDS "_LB")
  elseif(kernel_bounds MATCHES "^NO")
    set(TMP_BOUNDS "")
  else()
    message(FATAL_ERROR "Invalid bounds")
  endif()
  set(TMP_PRE "")
  set(TMP_POST "")
  if(NOT kernel_bounds MATCHES "_OCL1")
    set(TMP_PRE "#ifdef GPUCA_KRNL_NOOCL1\n")
    set(TMP_POST "#endif\n")
  endif()
  set(TMP_KERNEL "GPUCA_KRNL${TMP_BOUNDS}((${kernel_name}), (${kernel_type}), (${OPT1}), (${OPT2}), (${OPT3}))\n")
  separate_arguments(kernel_files NATIVE_COMMAND ${kernel_files})
  list(GET kernel_files 0 TMP_KERNEL_CLASS_FILE)
  if (TMP_KERNEL_CLASS_FILE STREQUAL "=")
    string(REGEX REPLACE ",.*$" "" TMP_KERNEL_CLASS_FILE "${kernel_name}")
  endif()
  set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNELS "${TMP_PRE}${TMP_KERNEL}${TMP_POST}")
  set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNEL_NAMES "${kernel_name}")
  set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNEL_INCLUDES "${TMP_KERNEL_CLASS_FILE}")
  set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNEL_FILES "${TMP_KERNEL_CLASS_FILE}.cxx")
  # add_custom_command OUTPUT option does not support target-dependend generator expressions, thus this workaround

  set(O2_GPU_KERNEL_TEMPLATE_FILES "GPUConstantMem.h")
  list(LENGTH kernel_files n)
  if(n GREATER 1)
    math(EXPR n "${n} - 1")
    foreach(i RANGE 1 ${n})
      list(GET kernel_files ${i} TMP_KERNEL_FILE_LIST)
      get_target_property(TMP_FILE_LIST O2_GPU_KERNELS O2_GPU_KERNELS_FILE_LIST_${TMP_KERNEL_FILE_LIST})
      if(NOT TMP_FILE_LIST)
        message(FATAL_ERROR "Invalid file list ${TMP_KERNEL_FILE_LIST}")
      endif()
      list(APPEND O2_GPU_KERNEL_TEMPLATE_FILES ${TMP_FILE_LIST})
    endforeach()
  endif()
  list(APPEND O2_GPU_KERNEL_TEMPLATE_FILES "${TMP_KERNEL_CLASS_FILE}.cxx")
  list(REMOVE_DUPLICATES O2_GPU_KERNEL_TEMPLATE_FILES)
  list(TRANSFORM O2_GPU_KERNEL_TEMPLATE_FILES APPEND "\"")
  list(TRANSFORM O2_GPU_KERNEL_TEMPLATE_FILES PREPEND "#include \"")
  list(JOIN O2_GPU_KERNEL_TEMPLATE_FILES "\n" O2_GPU_KERNEL_TEMPLATE_FILES)

  string(REPLACE ", " "_" TMP_FILENAME "${kernel_name}")
  if(CUDA_ENABLED)
    set(TMP_FILENAMEA "${O2_GPU_KERNEL_WRAPPER_FOLDER}/krnl_${TMP_FILENAME}.cu")
    set(O2_GPU_KERNEL_TEMPLATE_REPLACE "${TMP_KERNEL}")
    configure_file(${O2_GPU_BASE_DIR}/Base/cuda/GPUReconstructionCUDAkernel.template.cu ${TMP_FILENAMEA})
  endif()

  if(HIP_ENABLED)
    set(TMP_FILENAMEA "${O2_GPU_KERNEL_WRAPPER_FOLDER}/krnl_${TMP_FILENAME}.hip")
    set(O2_GPU_KERNEL_TEMPLATE_REPLACE "${TMP_KERNEL}")
    configure_file(${O2_GPU_BASE_DIR}/Base/hip/GPUReconstructionHIPkernel.template.hip ${TMP_FILENAMEA})
  endif()
endfunction()

function(o2_gpu_kernel_file_list list)
  define_property(TARGET PROPERTY O2_GPU_KERNELS_FILE_LIST_${list})
  list(LENGTH ARGV n)
  if(2 GREATER ${n})
    message(FATAL_ERROR "File list must contain at least one file")
  endif()
  if(list MATCHES "[^A-Z0-9_]")
    message(FATAL_ERROR "Invalid character in file list name ${list}")
  endif()
  math(EXPR n "${n} - 1")
  foreach(i RANGE 1 ${n})
    if(NOT ARGV${i} MATCHES "[^A-Z0-9_]")
      get_target_property(TMP_SUB_FILE_LIST O2_GPU_KERNELS O2_GPU_KERNELS_FILE_LIST_${ARGV${i}})
      if(NOT TMP_SUB_FILE_LIST)
        message(FATAL_ERROR "Invalid file list ${ARGV${i}}")
      endif()
      list(APPEND TMP_FILE_LIST ${TMP_SUB_FILE_LIST})
    else()
      list(APPEND TMP_FILE_LIST ${ARGV${i}})
      set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNEL_FILES "${ARGV${i}}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES TMP_FILE_LIST)
  set_property(TARGET O2_GPU_KERNELS PROPERTY O2_GPU_KERNELS_FILE_LIST_${list} "${TMP_FILE_LIST}")
endfunction()
