Ba2NiIrO6 Tutorial
==================

This chapter shows a complete, practical workflow for running Material-MC with the provided ``Ba2NiIrO6`` example.

Prerequisites
-------------

Before running this example, make sure:

- Material-MC has been built (the executable is ``build/MMC``).
- MPI is available on your system.
- You are in the project root directory.

Example files
-------------

The example folder contains two required input files:

- ``example/Ba2NiIrO6/POSCAR``: crystal structure (lattice vectors + atomic positions)
- ``example/Ba2NiIrO6/input.toml``: Monte Carlo and Hamiltonian settings

In this example:

- Method: ``ptmc`` (parallel tempering Monte Carlo)
- Temperature range: ``0.0`` to ``400.0``
- Temperature points: ``200``
- Supercell size: ``[8, 8, 4]``
- Magnetic elements: ``Ni`` and ``Ir``
- Model: ``Heisenberg``

Run the simulation
------------------

From the project root, run:

.. code-block:: bash

   mpirun -np 8 ./build/MMC \
     -c example/Ba2NiIrO6/POSCAR \
     -i example/Ba2NiIrO6/input.toml \
     -o Ba2NiIrO6_output.txt \
     -s Ba2NiIrO6_spin

Command-line options:

- ``-c``: structure file path (POSCAR)
- ``-i``: parameter file path (TOML)
- ``-o``: output table file name
- ``-s``: spin-output filename prefix

Check outputs
-------------

After a successful run, check:

- ``log.txt``: parsed settings and simulation log
- ``Ba2NiIrO6_output.txt``: thermodynamic results vs. temperature
- ``Ba2NiIrO6_spin*``: spin structure output files (with your chosen prefix)

What to tune next
-----------------

A typical next step is editing ``example/Ba2NiIrO6/input.toml``:

- ``relaxing_steps`` / ``counting_steps`` for convergence and statistics
- ``temperature_points_number`` for temperature resolution
- ``cell_number`` for system size effects
- exchange parameters in ``[[Elements.Neighbors]]`` for model studies

This gives you a reproducible baseline workflow: keep ``POSCAR`` fixed, tune ``input.toml``, and compare output curves from ``*_output.txt``.
