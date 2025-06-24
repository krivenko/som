# #############
# Debian Package

option(BUILD_DEBIAN_PACKAGE "Build a DEB package" OFF)
if(BUILD_DEBIAN_PACKAGE AND NOT IS_SUBPROJECT)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
  set(CPACK_GENERATOR "DEB")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "https://github.com/krivenko/${PROJECT_NAME}")
  execute_process(COMMAND dpkg --print-architecture
                  OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= ${REQUIRED_TRIQS_VERSION_MIN})")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  include(CPack)
endif()

# #############
# easyconfig file

configure_file(TRIQS-som-foss-2023a.eb.in
               TRIQS-som-${PROJECT_VERSION}-foss-2023a.eb)
