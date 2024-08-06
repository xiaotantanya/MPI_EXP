mpiexec --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include enp13s0f0,eno1,enp13s0f1 -n 2 --hostfile host ./odd_even_sort 256M
mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include eno1 -H pc -n 1 ./odd_even_sort 256M : --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include enp13s0f0,enp13s0f1 -H ubuntu -n 1 ./odd_even_sort 256M

mpirun --mca pml ob1 --mca btl tcp,self \
    -H pc --mca btl_tcp_if_include eno1 -n 1 ./odd_even_sort 256M : \
    -H ubuntu --mca btl_tcp_if_include enp13s0f0,enp13s0f1 -n 1 ./odd_even_sort 256M


## 多节点运行指令
在203，204服务器上
mpiexec --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include 202.38.247.204/24 -n 8 --hostfile host ./odd_even_sort 256M


## 注意
需要看懂https://github.com/open-mpi/ompi/issues/4963（已看懂）
