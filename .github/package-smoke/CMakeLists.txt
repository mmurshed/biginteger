cmake_minimum_required(VERSION 3.16)
project(BigMathPackageSmoke LANGUAGES CXX)

enable_testing()
find_package(bigmath CONFIG REQUIRED)
add_executable(package_smoke main.cpp)
target_link_libraries(package_smoke PRIVATE bigmath::bigmath)
add_test(NAME package_smoke COMMAND package_smoke)
