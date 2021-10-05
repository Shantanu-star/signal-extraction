
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(fit FitXicZerotoXiPiInvMass.C)
target_link_libraries(fit SignalExtraction ${ROOT_LIBRARIES} EG Minuit)

add_executable(test_exe test.cpp)
target_link_libraries(test_exe SignalExtraction ${ROOT_LIBRARIES} EG Minuit)
