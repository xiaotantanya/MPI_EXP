add_executable(test_send_recv test_send_recv_large.c)
target_link_libraries(test_send_recv MPI::MPI_C)

add_executable(send_recv send_recv.c)
target_link_libraries(send_recv MPI::MPI_C)