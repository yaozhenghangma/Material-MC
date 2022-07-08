import time
import os
import hybrid_cf
import hybrid_pf

time_start = time.time()
os.system("./pure_cpp")
time_end = time.time()
print("Time cost for pure c++: ", time_end-time_start, "s")

time_start = time.time()
hybrid_cf.monte_carlo_cpp_function()
time_end = time.time()
print("Time cost for hybrid with c++ Hamiltonian: ", time_end-time_start, "s")

time_start = time.time()
hybrid_pf.monte_carlo_python_function()
time_end = time.time()
print("Time cost for hybrid with python Hamiltonian: ", time_end-time_start, "s")