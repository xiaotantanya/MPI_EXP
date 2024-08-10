cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(MPI_EXP
        DESCRIPTION "MPI EXP"
        VERSION 0.0.1
        LANGUAGES C)

# 查找 MPI 包
find_package(MPI REQUIRED)
set(SRC_LIST)
# 生成动态库
add_library(MPI_EXP SHARED src/utils.c src/sort.c)
target_include_directories(MPI_EXP PUBLIC ${PROJECT_SOURCE_DIR}/include)

# add_subdirectory(src)
# message(STATUS ${SRC_LIST})
# set(MPI_EXP_LIB MPI_EXP)


target_link_libraries(MPI_EXP MPI::MPI_C)

# 添加可执行文件
add_executable(fast_sort_single fast_sort_single.c)
add_executable(odd_even_sort odd_even_sort_para.c)
add_executable(PSRS PSRS.c)

# 链接 MPI 库
target_link_libraries(fast_sort_single MPI_EXP)
target_link_libraries(odd_even_sort PUBLIC MPI_EXP)
target_link_libraries(PSRS PUBLIC MPI_EXP)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/host" "${CMAKE_CURRENT_BINARY_DIR}/host" COPYONLY)

add_subdirectory(test)