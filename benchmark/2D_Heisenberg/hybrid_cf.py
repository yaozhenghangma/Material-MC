import hybrid_cpp_function as hcf

def monte_carlo_cpp_function():
    T = 0
    outfile = open("output2.txt", "w")
    outfile.write("T\tenergy\n")
    for _ in range(0, 101):
        average_energy = hcf.MonteCarlo(T)
        outfile.write(T, "\t", average_energy, "\n")
        T += 1
    outfile.close()
