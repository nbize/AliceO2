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
function(o2_gpu_add_kernel kernel_name kernel_bounds kernel_type)
  math(EXPR TMP_CHK "${ARGC} & 1")
  if(NOT ${TMP_CHK})
    message(FATAL_ERROR "Invalid number of arguments to kernel ${TMP_CHK}, must be odd to have pairs of argument type, argument name")
  endif()
  list(LENGTH ARGV n)
  set(OPT1 "")
  set(OPT2 "")
  if(${n} GREATER 3)
    math(EXPR n "${n} - 1")
    foreach(i RANGE 3 ${n} 2)
      math(EXPR j "${i} + 1")
      if(${ARGV${i}} MATCHES "\\*$")
        string(APPEND OPT1 ",GPUPtr1(${ARGV${i}},${ARGV${j}})")
        string(APPEND OPT2 ",GPUPtr2(${ARGV${i}},${ARGV${j}})")
      else()
        string(APPEND OPT1 ",${ARGV${i}} ${ARGV${j}}")
        string(APPEND OPT2 ",${ARGV${j}}")
      endif()
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
  set(TMP_KERNEL "${TMP_PRE}GPUCA_KRNL${TMP_BOUNDS}((${kernel_name}), (${kernel_type}), (${OPT1}), (${OPT2}))\n${TMP_POST}")
  set_property(TARGET O2_GPU_KERNELS APPEND PROPERTY O2_GPU_KERNELS "${TMP_KERNEL}")
endfunction()
