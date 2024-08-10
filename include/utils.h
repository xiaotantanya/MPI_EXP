#include <stdint.h>
#include <mpi.h>

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char **argv, int rank, int is_test);
void check_correct(float* sort_arr, char* check_file, uint64_t N);
void gen_correct_file(char* check_file, uint64_t N, int rank, int force, int is_test, float seed);

int compare(const void *a, const void *b);
int re_compare(const void *a, const void *b);

int MPI_Send_Large(void *buf, uint64_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int element_size);
int MPI_Recv_Large(void *buf, uint64_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status, int element_size);
int MPI_Sendrecv_Large(void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                        int dest, int sendtag, void *recvbuf, uint64_t recvcount,
                        MPI_Datatype recvtype, int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status, int element_size);
