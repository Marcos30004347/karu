cmake_minimum_required(VERSION 3.10)

project(CameraSample)
add_executable(CameraSample camera/Camera.cpp)
target_link_libraries(CameraSample karu)
add_test(NAME CameraSample COMMAND CameraSample)

