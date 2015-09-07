set(OPENGLLIBS ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLEW_LIBRARIES} ${OPENGL_LIBRARY})

if(WIN32)
message(STATUS "Win32: adding glfw32dll and glew32")
set(OPENGLLIBS glfw3dll glew32 ${OPENGLLIBS})
else()
set(OPENGLLIBS glfw ${OPENGLLIBS})
endif()


set(01_srcs  01_raytrace.cpp)                               # 01_raytrace
add_executable(01_raytrace ${01_srcs})                      # 01_raytrace
target_link_libraries(01_raytrace common ${OPENGLLIBS})     # 01_raytrace
SOURCE_GROUP("" FILES ${01_srcs})                           # 01_raytrace








if(CMAKE_GENERATOR STREQUAL "Xcode")
    set_property(TARGET  01_raytrace      PROPERTY XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++11)
    set_property(TARGET  01_raytrace      PROPERTY XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)
endif(CMAKE_GENERATOR STREQUAL "Xcode")


