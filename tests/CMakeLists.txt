Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7
    )

FetchContent_MakeAvailable(Catch2)

file(GLOB tests_src "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_tadah.lammps ${tests_src})
target_link_libraries(test_tadah.lammps PRIVATE Catch2::Catch2 tadah.lammps tadah.core tadah.models)
    #if(BUILD_SHARED_LIBS)
    #    target_link_libraries(test_tadah.lammps PRIVATE Catch2::Catch2 tadah.lammps)
    #else()
    #    add_library(tadah.lammps STATIC IMPORTED )
    #    set_target_properties(tadah.lammps PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/libtadah.lammps.a )
    #    #target_link_libraries(test_tadah.lammps PRIVATE Catch2::Catch2 tadah.lammps)
    #    #target_link_libraries(test_tadah.lammps PRIVATE Catch2::Catch2 ${CMAKE_BINARY_DIR}/libtadah.lammps.a)
    #endif()

find_package(LAPACK)
target_link_libraries (test_tadah.lammps PRIVATE ${LAPACK_LIBRARIES})

include(Catch.cmake)
catch_discover_tests(test_tadah.lammps)
