mpiexec --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include enp13s0f0,eno1,enp13s0f1 -n 2 --hostfile host ./odd_even_sort 256M
mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include eno1 -H pc -n 1 ./odd_even_sort 256M : --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include enp13s0f0,enp13s0f1 -H ubuntu -n 1 ./odd_even_sort 256M

mpirun --mca pml ob1 --mca btl tcp,self \
    -H pc --mca btl_tcp_if_include eno1 -n 1 ./odd_even_sort 256M : \
    -H ubuntu --mca btl_tcp_if_include enp13s0f0,enp13s0f1 -n 1 ./odd_even_sort 256M

mpiexec --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include 202.38.247.204/24 -n 8 --hostfile host ./odd_even_sort 256M


## 注意
需要看懂https://github.com/open-mpi/ompi/issues/4963

## MPI_Send 和 MPI_Recv死锁问题
相关链接：

https://mpitutorial.com/tutorials/point-to-point-communication-application-random-walk/zh_cn/
https://blog.csdn.net/susan_wang1/article/details/50068439

尽管 MPI_Send 是一个阻塞调用，但是 MPI 规范 表明 MPI_Send 会一直阻塞，直到可以回收发送缓冲区为止。 这意味着当网络可以缓冲消息时，MPI_Send 将返回。 如果发送最终无法被网络缓冲，它们将一直阻塞直到发布匹配的接收。 在我们的例子中，有足够多的小发送和频繁匹配的接收而不必担心死锁，但是，永远不该假定有足够大的网络缓冲区。

结论：不一定会报错，缓冲区没满的时候就不会报错，传递数据太多就会导致死锁。

![死锁示意图](/image/image.png)