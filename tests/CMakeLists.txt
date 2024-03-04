Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7
    )

FetchContent_MakeAvailable(Catch2)

file(GLOB tests_src "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_tadah_lammps ${tests_src})
target_link_libraries(test_tadah_lammps PRIVATE Catch2::Catch2 tadah_lammps)

include(Catch.cmake)
catch_discover_tests(test_tadah_lammps)
