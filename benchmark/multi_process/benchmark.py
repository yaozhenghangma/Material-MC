import time
import os

time_start = time.time()
os.system("mpirun -n 1 python3 heisenberg_mpi.py")
time_end = time.time()
print("Time cost for 1 process: ", time_end-time_start, "s")

time_start = time.time()
os.system("mpirun -n 2 python3 heisenberg_mpi.py")
time_end = time.time()
print("Time cost for 2 process: ", time_end-time_start, "s")

time_start = time.time()
os.system("mpirun -n 4 python3 heisenberg_mpi.py")
time_end = time.time()
print("Time cost for 4 process: ", time_end-time_start, "s")