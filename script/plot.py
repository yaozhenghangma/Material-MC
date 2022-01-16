import numpy as np
import matplotlib.pyplot as plt

###################################################### Data preprocess
path = 'output.txt'
with open(path, encoding='utf-8') as f:
    data = np.loadtxt(path, delimiter='\t', skiprows=1)

size = data.shape[1]
if size==5:
    projection = False
else:
    projection = True

temperature = data[1:, 0]
energy = data[1:, 1]
Cv = data[1:, 2]
moment = data[1:, 3]
chi = data[1:, 4]
if projection:
    moment_projection = data[1:, 5]
    chi_projection = data[1:, 6]

temperature_start = np.min(temperature)
temperature_end = np.max(temperature)
energy_max = np.max(energy)
energy_min = np.min(energy)
Cv_max = np.max(Cv)
moment_max = np.max(moment)
chi_max = np.max(chi)
if projection:
    moment_projection_max = np.max(moment_projection)
    chi_projection_max = np.max(chi_projection)

###################################################### Font and size
font = {'family' : 'Times New Roman',
'weight' : 'regular',
'size' : 23
}

font_italic = {'family' : 'Times New Roman',
'style' : "italic",
'weight' : 'regular',
'size' : 23
}

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times"],
})

figure_size = (8,6)

###################################################### Moment and magnetic susceptibility
fig, ax = plt.subplots(figsize=figure_size) # (width, height) in inches
par1 = ax.twinx()

# tick setting
ax.spines["left"]
ax.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
ax.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
par1.spines["right"]
par1.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
par1.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
ax.set_xlim(temperature_start, temperature_end)
ax.set_ylim(0, moment_max)
par1.set_ylim(0, chi_max)

# label setting
ax.set_xlabel(r'$T$ $(\mathrm{K})$', font)
ax.set_ylabel(r'$M$ $(\mu_\mathrm{B})$', font)
par1.set_ylabel(r'$\chi$', font)

# frame setting
frameSize = 1.5
ax.spines['left'].set_linewidth(frameSize)
ax.spines['right'].set_linewidth(frameSize)
ax.spines['top'].set_linewidth(frameSize)
ax.spines['bottom'].set_linewidth(frameSize)
par1.spines['left'].set_linewidth(frameSize)
par1.spines['right'].set_linewidth(frameSize)
par1.spines['top'].set_linewidth(frameSize)
par1.spines['bottom'].set_linewidth(frameSize)

p1, = ax.plot(temperature, moment, color="blue", label=r'$M$')
p2, = par1.plot(temperature, chi, color="red", label=r'$\chi$')

lns = [p1, p2]
ax.legend(handles=lns, loc='best', frameon=False, handlelength=1.2, prop=font)

fig.tight_layout()
plt.savefig("moment.png", format='png', dpi=600)

###################################################### Projected moment and magnetic susceptibility
if projection:
    fig, ax = plt.subplots(figsize=figure_size) # (width, height) in inches
    par1 = ax.twinx()

    # tick setting
    ax.spines["left"]
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
    ax.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
    par1.spines["right"]
    par1.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
    par1.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
    ax.set_xlim(temperature_start, temperature_end)
    ax.set_ylim(0, moment_projection_max)
    par1.set_ylim(0, chi_projection_max)

    # label setting
    ax.set_xlabel(r'$T$ $(\mathrm{K})$', font)
    ax.set_ylabel(r'$M$ $(\mu_\mathrm{B})$', font)
    par1.set_ylabel(r'$\chi$', font)

    # frame setting
    frameSize = 1.5
    ax.spines['left'].set_linewidth(frameSize)
    ax.spines['right'].set_linewidth(frameSize)
    ax.spines['top'].set_linewidth(frameSize)
    ax.spines['bottom'].set_linewidth(frameSize)
    par1.spines['left'].set_linewidth(frameSize)
    par1.spines['right'].set_linewidth(frameSize)
    par1.spines['top'].set_linewidth(frameSize)
    par1.spines['bottom'].set_linewidth(frameSize)

    p1, = ax.plot(temperature, moment_projection, color="blue", label=r'$M$')
    p2, = par1.plot(temperature, chi_projection, color="red", label=r'$\chi$')

    lns = [p1, p2]
    ax.legend(handles=lns, loc='best', frameon=False, handlelength=1.2, prop=font)

    fig.tight_layout()
    plt.savefig("projected_moment.png", format='png', dpi=600)

###################################################### Energy and specific heat
fig, ax = plt.subplots(figsize=figure_size) # (width, height) in inches
par1 = ax.twinx()

# tick setting
ax.spines["left"]
ax.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
ax.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
par1.spines["right"]
par1.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
par1.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=20)
ax.set_xlim(temperature_start, temperature_end)
ax.set_ylim(energy_min, energy_max)
par1.set_ylim(0, Cv_max)

# label setting
ax.set_xlabel(r'$T$ $(\mathrm{K})$', font)
ax.set_ylabel(r'$E$ $(\mathrm{meV})$', font)
par1.set_ylabel(r'$C_\mathrm{V}$', font)

# frame setting
frameSize = 1.5
ax.spines['left'].set_linewidth(frameSize)
ax.spines['right'].set_linewidth(frameSize)
ax.spines['top'].set_linewidth(frameSize)
ax.spines['bottom'].set_linewidth(frameSize)
par1.spines['left'].set_linewidth(frameSize)
par1.spines['right'].set_linewidth(frameSize)
par1.spines['top'].set_linewidth(frameSize)
par1.spines['bottom'].set_linewidth(frameSize)

p1, = ax.plot(temperature, energy, color="blue", label=r'$E$')
p2, = par1.plot(temperature, Cv, color="red", label=r'$C_\mathrm{V}$')

lns = [p1, p2]
ax.legend(handles=lns, loc='best', frameon=False, handlelength=1.2, prop=font)

fig.tight_layout()
plt.savefig("energy.png", format='png', dpi=600)