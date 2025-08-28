cmake_minimum_required(VERSION 3.15)
project(imgui_kline_viewer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(IMGUI_KLINE_BUILD_DEMO "Build with ImGui demo window" OFF)

# Paths (now relative to repo root)
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/imgui)
set(IMGUI_BACKENDS_DIR ${IMGUI_DIR}/backends)

# ImGui sources
set(IMGUI_SOURCES
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_BACKENDS_DIR}/imgui_impl_glfw.cpp
    ${IMGUI_BACKENDS_DIR}/imgui_impl_opengl3.cpp
)
if (IMGUI_KLINE_BUILD_DEMO)
    list(APPEND IMGUI_SOURCES ${IMGUI_DIR}/imgui_demo.cpp)
endif()

add_library(imgui STATIC ${IMGUI_SOURCES})
target_include_directories(imgui PUBLIC ${IMGUI_DIR} ${IMGUI_BACKENDS_DIR})
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GL3W)

# Find GLFW (pkg-config preferred on Linux)
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
    pkg_check_modules(GLFW3 QUIET glfw3)
endif()

set(GLFW_LIBS)
if (GLFW3_FOUND)
    message(STATUS "Using GLFW via pkg-config")
    include_directories(${GLFW3_INCLUDE_DIRS})
    link_directories(${GLFW3_LIBRARY_DIRS})
    set(GLFW_LIBS ${GLFW3_LIBRARIES} ${GLFW3_LDFLAGS_OTHER})
else()
    find_package(glfw3 QUIET)
    if (TARGET glfw)
        set(GLFW_LIBS glfw)
    elseif (TARGET glfw3)
        set(GLFW_LIBS glfw3)
    else()
        message(FATAL_ERROR "GLFW not found. Install libglfw3-dev and try again.")
    endif()
endif()

# App
file(GLOB_RECURSE APP_SOURCES CONFIGURE_DEPENDS kline_viewer/src/*.cpp)
add_executable(kline_viewer ${APP_SOURCES})

target_include_directories(kline_viewer PRIVATE kline_viewer/include ${IMGUI_DIR} ${IMGUI_BACKENDS_DIR})

# OpenGL
find_package(OpenGL REQUIRED)

# Link
target_link_libraries(kline_viewer PRIVATE imgui ${GLFW_LIBS} OpenGL::GL)

# SQLite3 (for DB data source)
find_package(SQLite3 QUIET)
if (SQLite3_FOUND)
    target_link_libraries(kline_viewer PRIVATE SQLite::SQLite3)
else()
    # Fallback: try linking with system sqlite3 (Linux)
    target_link_libraries(kline_viewer PRIVATE sqlite3)
endif()

# Linux specifics
if (UNIX AND NOT APPLE)
    target_link_libraries(kline_viewer PRIVATE dl pthread)
endif()

# Use builtin OpenGL loader from ImGui backend
target_compile_definitions(kline_viewer PRIVATE IMGUI_IMPL_OPENGL_LOADER_GL3W)

# Install
install(TARGETS kline_viewer RUNTIME DESTINATION bin)
