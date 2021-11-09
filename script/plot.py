import numpy as np
import matplotlib.pyplot as plt

path = 'output.txt'
with open(path, encoding='utf-8') as f:
    data = np.loadtxt(path, delimiter='\t', skiprows=1)



font = {'family' : 'Times New Roman',
'weight' : 'bold',
'size' : 23
}

################################# moment
# More versatile wrapper
fig, host = plt.subplots(figsize=(8,5)) # (width, height) in inches
par1 = host.twinx()

xmin = np.min(data[1:, 0])
xmax = np.max(data[1:, 0])

moment_max = np.max(data[1:, 3])
Ki_max = np.max(data[1:, 4])
    
host.set_xlim(xmin, xmax)
host.set_ylim(0, moment_max*1.1)
par1.set_ylim(0, Ki_max*2.2)
    
host.set_xlabel(r'$T(K)$', font)
host.set_ylabel(r'$M(\mu_B/f.u.)$', font)
par1.set_ylabel(r'$\chi$', font)

color1 = plt.cm.viridis(0)
color2 = plt.cm.viridis(0.5)

par1.spines['right']
host.spines["left"]
host.tick_params(axis='both', direction='in', labelsize=20)
par1.tick_params(axis='both', direction='in', labelsize=20)

p1, = host.plot(data[1:, 0], data[1:, 3],    color=color1, label=r'$M$')
p2, = par1.plot(data[1:, 0], data[1:, 4],    color=color2, label=r'$\chi$')

lns = [p1, p2]
host.legend(handles=lns, loc='best', prop=font)

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())

# Adjust spacings w.r.t. figsize
fig.tight_layout()
# Alternatively: bbox_inches='tight' within the plt.savefig function 
#                (overwrites figsize)

# Best for professional typesetting, e.g. LaTeX
plt.savefig("moment.pdf")
# For raster graphics use the dpi argument. E.g. '[...].png", dpi=200)'




################################# moment x
# More versatile wrapper
fig, host = plt.subplots(figsize=(8,5)) # (width, height) in inches
par1 = host.twinx()

xmin = np.min(data[1:, 0])
xmax = np.max(data[1:, 0])

moment_max = np.max(data[1:, 5])
Ki_max = np.max(data[1:, 6])
    
host.set_xlim(xmin, xmax)
host.set_ylim(0, moment_max*1.1)
par1.set_ylim(0, Ki_max*2.2)
    
host.set_xlabel(r'$T(K)$', font)
host.set_ylabel(r'$M_x(\mu_B/f.u.)$', font)
par1.set_ylabel(r'$\chi_x$', font)

color1 = plt.cm.viridis(0)
color2 = plt.cm.viridis(0.5)

par1.spines['right']
host.spines["left"]
host.tick_params(axis='both', direction='in', labelsize=20)
par1.tick_params(axis='both', direction='in', labelsize=20)

p1, = host.plot(data[1:, 0], data[1:, 5],    color=color1, label=r'$M_x$')
p2, = par1.plot(data[1:, 0], data[1:, 6],    color=color2, label=r'$\chi_x$')

lns = [p1, p2]
host.legend(handles=lns, loc='best', prop=font)

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())

# Adjust spacings w.r.t. figsize
fig.tight_layout()
# Alternatively: bbox_inches='tight' within the plt.savefig function 
#                (overwrites figsize)

# Best for professional typesetting, e.g. LaTeX
plt.savefig("moment_x.pdf")
# For raster graphics use the dpi argument. E.g. '[...].png", dpi=200)'



################################# moment y
# More versatile wrapper
fig, host = plt.subplots(figsize=(8,5)) # (width, height) in inches
par1 = host.twinx()

xmin = np.min(data[1:, 0])
xmax = np.max(data[1:, 0])

moment_max = np.max(data[1:, 7])
Ki_max = np.max(data[1:, 8])
    
host.set_xlim(xmin, xmax)
host.set_ylim(0, moment_max*1.1)
par1.set_ylim(0, Ki_max*2.2)
    
host.set_xlabel(r'$T(K)$', font)
host.set_ylabel(r'$M_y(\mu_B/f.u.)$', font)
par1.set_ylabel(r'$\chi_y$', font)

color1 = plt.cm.viridis(0)
color2 = plt.cm.viridis(0.5)

par1.spines['right']
host.spines["left"]
host.tick_params(axis='both', direction='in', labelsize=20)
par1.tick_params(axis='both', direction='in', labelsize=20)

p1, = host.plot(data[1:, 0], data[1:, 7],    color=color1, label=r'$M_y$')
p2, = par1.plot(data[1:, 0], data[1:, 8],    color=color2, label=r'$\chi_y$')

lns = [p1, p2]
host.legend(handles=lns, loc='best', prop=font)

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())

# Adjust spacings w.r.t. figsize
fig.tight_layout()
# Alternatively: bbox_inches='tight' within the plt.savefig function 
#                (overwrites figsize)

# Best for professional typesetting, e.g. LaTeX
plt.savefig("moment_y.pdf")
# For raster graphics use the dpi argument. E.g. '[...].png", dpi=200)'



################################# moment z
# More versatile wrapper
fig, host = plt.subplots(figsize=(8,5)) # (width, height) in inches
par1 = host.twinx()

xmin = np.min(data[1:, 0])
xmax = np.max(data[1:, 0])

moment_max = np.max(data[1:, 9])
Ki_max = np.max(data[1:, 10])
    
host.set_xlim(xmin, xmax)
host.set_ylim(0, moment_max*1.1)
par1.set_ylim(0, Ki_max*2.2)
    
host.set_xlabel(r'$T(K)$', font)
host.set_ylabel(r'$M_z(\mu_B/f.u.)$', font)
par1.set_ylabel(r'$\chi_z$', font)

color1 = plt.cm.viridis(0)
color2 = plt.cm.viridis(0.5)

par1.spines['right']
host.spines["left"]
host.tick_params(axis='both', direction='in', labelsize=20)
par1.tick_params(axis='both', direction='in', labelsize=20)

p1, = host.plot(data[1:, 0], data[1:, 9],    color=color1, label=r'$M_z$')
p2, = par1.plot(data[1:, 0], data[1:, 10],    color=color2, label=r'$\chi_z$')

lns = [p1, p2]
host.legend(handles=lns, loc='best', prop=font)

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())

# Adjust spacings w.r.t. figsize
fig.tight_layout()
# Alternatively: bbox_inches='tight' within the plt.savefig function 
#                (overwrites figsize)

# Best for professional typesetting, e.g. LaTeX
plt.savefig("moment_z.pdf")
# For raster graphics use the dpi argument. E.g. '[...].png", dpi=200)'


########################### Cv
fig, host = plt.subplots(figsize=(8,5))
host.set_xlabel(r'$T(K)$', font)
host.set_ylabel(r'$Cv(meV/K)$', font)
host.tick_params(axis='both', direction='in', labelsize=20)
host.plot(data[1:, 0], data[1:, 2],    color=color1, label=r'$Cv$')
plt.savefig("Cv.pdf")