
file(GLOB tests_src "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_tadah.lammps ${tests_src})
target_link_libraries(test_tadah.lammps PRIVATE Catch2 tadah.lammps tadah.core tadah.models)

find_package(LAPACK)
target_link_libraries (test_tadah.lammps PRIVATE ${LAPACK_LIBRARIES})

include(Catch.cmake)
catch_discover_tests(test_tadah.lammps)
