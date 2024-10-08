cmake_minimum_required(VERSION 3.30)
project(lotka-volterra)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS multi_array)
find_package(GSL COMPONENTS gsl_errno gsl_odeiv)

if(Boost_FOUND AND GSL_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${GSL_INCLUDE_DIRS})
    add_executable(lotka-volterra run_lotka_voltera.cpp)
    target_link_libraries(lotka-volterra ${Boost_LIBRARIES})
    target_link_libraries(lotka-volterra ${GSL_LIBRARIES})
endif()
