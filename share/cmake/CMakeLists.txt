configure_file(som-config.cmake.in som-config.cmake @ONLY)
configure_file(som-config-version.cmake.in som-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/som-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/som-config-version.cmake
  DESTINATION lib/cmake/som
)

install(EXPORT som-targets NAMESPACE som:: DESTINATION lib/cmake/som)
