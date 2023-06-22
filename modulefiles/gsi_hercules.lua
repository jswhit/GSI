help([[
]])

prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/spack-stack/spack-stack-1.4.0-hercules/envs/unified-env-v2/install/modulefiles/Core")
prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/spack-stack/spack-stack-1.4.0-hercules/envs/unified-env-v2/install/modulefiles/intel-oneapi-mpi/2021.7.1/intel/2021.7.1")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.7.1"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.7.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_mkl_ver=os.getenv("stack_mkl_ver") or "2022.2.1"
load(pathJoin("intel-oneapi-mkl", stack_mkl_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.14"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("gsi_common")

---setenv("CC", "icc")
---setenv("CXX", "icpc")
---setenv("FC", "ifort")
---setenv("CMAKE_Platform", "hercules.intel")

pushenv("CFLAGS", "-xHOST")
pushenv("FFLAGS", "-xHOST")

pushenv("GSI_BINARY_SOURCE_DIR", "/work/noaa/global/glopara/fix/gsi/20230601")

whatis("Description: GSI environment on Orion with Intel Compilers")
