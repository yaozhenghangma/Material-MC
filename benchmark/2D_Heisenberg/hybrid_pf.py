import hybrid_python_function as hpf

def energy(lattice, x, y):
    n = 5
    J = -1.0
    energy = 0
    if x==0:
        energy += J*(
            lattice[x][y][0] * lattice[x+1][y][0] +
            lattice[x][y][1] * lattice[x+1][y][1] +
            lattice[x][y][2] * lattice[x+1][y][2])
        energy += J*(
            lattice[x][y][0] * lattice[n-1][y][0] +
            lattice[x][y][1] * lattice[n-1][y][1] +
            lattice[x][y][2] * lattice[n-1][y][2])
    elif x==n-1:
        energy += J*(
            lattice[x][y][0] * lattice[0][y][0] +
            lattice[x][y][1] * lattice[0][y][1] +
            lattice[x][y][2] * lattice[0][y][2])
        energy += J*(
            lattice[x][y][0] * lattice[x-1][y][0] +
            lattice[x][y][1] * lattice[x-1][y][1] +
            lattice[x][y][2] * lattice[x-1][y][2])
    else:
        energy += J*(
            lattice[x][y][0] * lattice[x+1][y][0] +
            lattice[x][y][1] * lattice[x+1][y][1] +
            lattice[x][y][2] * lattice[x+1][y][2])
        energy += J*(
            lattice[x][y][0] * lattice[x-1][y][0] +
            lattice[x][y][1] * lattice[x-1][y][1] +
            lattice[x][y][2] * lattice[x-1][y][2])

    if y==0:
        energy += J*(
            lattice[x][y][0] * lattice[x][y+1][0] +
            lattice[x][y][1] * lattice[x][y+1][1] +
            lattice[x][y][2] * lattice[x][y+1][2])
        energy += J*(
            lattice[x][y][0] * lattice[x][n-1][0] +
            lattice[x][y][1] * lattice[x][n-1][1] +
            lattice[x][y][2] * lattice[x][n-1][2])
    elif y==n-1:
        energy += J*(
            lattice[x][y][0] * lattice[x][0][0] +
            lattice[x][y][1] * lattice[x][0][1] +
            lattice[x][y][2] * lattice[x][0][2])
        energy += J*(
            lattice[x][y][0] * lattice[x][y-1][0] +
            lattice[x][y][1] * lattice[x][y-1][1] +
            lattice[x][y][2] * lattice[x][y-1][2])
    else:
        energy += J*(
            lattice[x][y][0] * lattice[x][y+1][0] +
            lattice[x][y][1] * lattice[x][y+1][1] +
            lattice[x][y][2] * lattice[x][y+1][2])
        energy += J*(
            lattice[x][y][0] * lattice[x][y-1][0] +
            lattice[x][y][1] * lattice[x][y-1][1] +
            lattice[x][y][2] * lattice[x][y-1][2])
    return energy

def monte_carlo_python_function():
    T = 0
    outfile = open("output3.txt", "w")
    outfile.write("T\tenergy\n")
    for _ in range(0, 101):
        average_energy = hpf.MonteCarlo(T, energy)
        outfile.write(str(T)+"\t"+str(average_energy)+"\n")
        T += 1
    outfile.close()
