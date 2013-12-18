SET(CMAKE_CXX_COMPILER "/usr/lib64/openmpi/bin/mpicxx")
SET(CMAKE_CXX_COMPILER_ARG1 "")
SET(CMAKE_CXX_COMPILER_ID "GNU")
SET(CMAKE_CXX_COMPILER_VERSION "4.4.7")
SET(CMAKE_CXX_PLATFORM_ID "Linux")

SET(CMAKE_AR "/usr/bin/ar")
SET(CMAKE_RANLIB "/usr/bin/ranlib")
SET(CMAKE_LINKER "/usr/bin/ld")
SET(CMAKE_COMPILER_IS_GNUCXX 1)
SET(CMAKE_CXX_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_CXX_COMPILER_ID_RUN 1)
SET(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
SET(CMAKE_CXX_LINKER_PREFERENCE 30)
SET(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
SET(CMAKE_CXX_SIZEOF_DATA_PTR "8")
SET(CMAKE_CXX_COMPILER_ABI "ELF")
SET(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

IF(CMAKE_CXX_SIZEOF_DATA_PTR)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
ENDIF(CMAKE_CXX_SIZEOF_DATA_PTR)

IF(CMAKE_CXX_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
ENDIF(CMAKE_CXX_COMPILER_ABI)

IF(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  SET(CMAKE_LIBRARY_ARCHITECTURE "")
ENDIF()

SET(CMAKE_CXX_HAS_ISYSROOT "")


SET(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpi_cxx;mpi;dl;stdc++;m;pthread;c")
SET(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/usr/lib64/openmpi/lib;/usr/lib/gcc/x86_64-redhat-linux/4.4.7;/usr/lib64;/lib64;/usr/lib")
