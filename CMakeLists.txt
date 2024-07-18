cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(serial_code)

# 指定要生成的可执行文件
add_executable(serial_code gen_data.c serial_code.c)

# 查找 MPI 包
find_package(MPI REQUIRED)

# 添加可执行文件
add_executable(odd_even_sort gen_data.c odd_even_sort_para.c)
# add_executable(hello_mpi hello_mpi.c)
# 链接 MPI 库
target_link_libraries(odd_even_sort MPI::MPI_C)
# target_link_libraries(hello_mpi MPI::MPI_C)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/host" "${CMAKE_CURRENT_BINARY_DIR}/host" COPYONLY)