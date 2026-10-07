# Catch2 and its test discovery come from cmake-tools (tadah_setup_tests()).
tadah_setup_tests()

file(GLOB tests_src CONFIGURE_DEPENDS "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_tadah.lammps ${tests_src})
target_link_libraries(test_tadah.lammps PRIVATE Catch2 tadah.lammps)

include(Catch)
catch_discover_tests(test_tadah.lammps)

# tadah.lammps_headers (ctest -L headers): every header of include/
# compiles on its own and twice in a row (cmake-tools'
# tadah_add_header_check). Before the fix: passes as well (a pin: no
# Tadah.LAMMPS header failed when the check was added).
tadah_add_header_check(tadah.lammps)
