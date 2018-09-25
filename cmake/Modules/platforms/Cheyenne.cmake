function (setCheyenne)
  message("Setting paths for Cheyenne")
# if( NOT DEFINED ENV{WRFPATH} )
#   set(WRFPATH "/glade/p/work/wrfhelp/PRE_COMPILED_CODE_CHEYENNE/WRFV3.9_intel_dmpar_large-file" PARENT_SCOPE )
# else()
#   set(WRFPATH $ENV{WRFPATH} PARENT_SCOPE )
# endif()
  set(BUILD_CORELIBS "ON" )
  set(BUILD_UTIL "OFF" CACHE INTERNAL "" )
  set(BUILD_BUFR "ON" CACHE INTERNAL "")
  set(BUILD_SFCIO "ON" CACHE INTERNAL "")
  set(BUILD_SIGIO "ON" CACHE INTERNAL "")
  set(BUILD_W3EMC "ON" CACHE INTERNAL "")
  set(BUILD_W3NCO "ON" CACHE INTERNAL "")
  set(BUILD_BACIO "ON" CACHE INTERNAL "")
  set(BUILD_CRTM "ON" CACHE INTERNAL "")
  set(BUILD_SP "ON" CACHE INTERNAL "")
  set(BUILD_NEMSIO "ON" CACHE INTERNAL "")
  set(ENV{MPI_HOME} $ENV{MPI_ROOT} )
endfunction()

