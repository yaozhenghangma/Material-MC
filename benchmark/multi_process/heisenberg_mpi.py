from mpi4py import MPI
import heisenberg as hs

def monte_carlo(T_max:int=100):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size ()
    energies = []
    T = rank
    for _ in range(0, int(100/size)):
        energy = hs.MonteCarlo(float(T))
        energies.append(energy)
        T += size

    energies = comm.gather(energies, root=0)
    if rank == 0:
        outfile = open("output"+str(size)+".txt", "w")
        outfile.write("T\tenergy\n")
        for j in range(0, int(100/size)):
            for i in range(0, size):
                outfile.write(str(j*size+i)+"\t"+str(energies[i][j])+"\n")
        outfile.close()

monte_carlo()