How to Configure input.toml
===========================

This tutorial explains how to set up ``input.toml`` for Material-MC.

Overview
--------

``input.toml`` controls:

- Monte Carlo schedule
- supercell size and neighbor matching tolerance
- magnetic elements and exchange shells
- Hamiltonian/model selection
- spin initialization
- output options

A complete template
-------------------

.. code-block:: toml

   [MonteCarlo]
   start_temperature = 0.0
   end_temperature = 400.0
   temperature_points_number = 200
   relaxing_steps = 8000
   counting_steps = 6000
   flipping_number = 6000
   method = "ptmc"            # "classical" or "ptmc"
   exchange_step = 1          # used by parallel tempering

   [Lattice]
   cell_number = [ 8, 8, 4 ]
   tolerance = 0.01

   [[Elements]]
   name = "Ni"
   spin = 1.0
   anisotropic_factor = [ 1.0, 1.0, 1.0 ]
   [[Elements.Neighbors]]
   name = "Ir"
   exchange_parameter = -7.5
   distance = 4.02
   [[Elements.Neighbors]]
   name = "Ni"
   exchange_parameter = 0.0
   distance = 5.68

   [[Elements]]
   name = "Ir"
   spin = 1.5
   anisotropic_factor = [ 1.0, 1.0, 1.0 ]
   [[Elements.Neighbors]]
   name = "Ni"
   exchange_parameter = -7.5
   distance = 4.02
   [[Elements.Neighbors]]
   name = "Ir"
   exchange_parameter = 0.0
   distance = 5.68

   [Hamiltonian]
   model = "Heisenberg"      # "Heisenberg", "Ising", or "Kitaev-Heisenberg"
   custom = false
   magnetic_field = [ 0.0, 0.0, 0.0 ]
   anisotropy = [ 0.0, 0.0, 0.0 ]
   KitaevEpsilon = 1e-12       # KH-only direction epsilon (optional)
   KitaevTolerance = 1e-6      # KH-only direction grouping tolerance (optional)
   J = 0.0                    # KH-only global coupling (optional)
   K = 0.0                    # KH-only global coupling (optional)
   G = 0.0                    # KH-only global coupling (optional)
   Gp = 0.0                   # KH-only global coupling (optional)

   [Hamiltonian.BondTypeDirection]
   type1 = "x"               # KH-only required mapping
   type2 = "y"               # KH-only required mapping
   type3 = "z"               # KH-only required mapping

   [Initialization]
   angleA = [ 0.0, 0.0, 0.0 ]
   angleB = [ 0.0, 0.0, 0.0 ]
   angleC = [ 0.0, 0.0, 0.0 ]

   [Output]
   magnifying_factor = 2.0
   ground_state = false

Section-by-section guide
------------------------

MonteCarlo
~~~~~~~~~~

- ``start_temperature``, ``end_temperature``: temperature range.
- ``temperature_points_number``: number of temperature points.
- ``relaxing_steps``: thermalization steps at each temperature.
- ``counting_steps``: sampling steps for averages.
- ``flipping_number``: local updates per MC step.
- ``method``: ``classical`` or ``ptmc``.
- ``exchange_step``: replica-exchange interval for PTMC.

.. note::
   Keep ``temperature_points_number > 1`` to avoid invalid temperature-step calculation.

Lattice
~~~~~~~

- ``cell_number = [nx, ny, nz]``: supercell expansion.
- ``tolerance``: relative tolerance used when matching neighbor-shell distances.

Elements and Neighbors
~~~~~~~~~~~~~~~~~~~~~~

Each ``[[Elements]]`` block defines one magnetic species.

- ``name``: element symbol; must match POSCAR naming.
- ``spin``: spin magnitude for this species.
- ``anisotropic_factor``: per-species spin scaling along x/y/z.

Each ``[[Elements.Neighbors]]`` entry defines one interaction shell:

- ``name``: neighbor element type.
- ``distance``: target shell distance.
- ``exchange_parameter``: exchange constant for this shell.

.. important::
   The magnetic-element definitions should be consistent with POSCAR element order and naming.

Hamiltonian
~~~~~~~~~~~

- ``model``: ``Heisenberg``, ``Ising``, or ``Kitaev-Heisenberg``.
- ``custom``: use custom Hamiltonian from ``custom/`` when ``true``.
- ``magnetic_field = [Bx, By, Bz]``: external field vector.
- ``anisotropy = [Dx, Dy, Dz]``: anisotropy constants.

When ``model = "Kitaev-Heisenberg"`` (also accepts ``"KH"``), you can additionally provide one global coupling set:

- ``KitaevEpsilon``: KH bond-direction normalization/sign epsilon (optional, default ``1e-12``)
- ``KitaevTolerance``: KH near-parallel grouping tolerance (optional, default ``1e-6``; must satisfy ``0 < value < 1``)
- ``J``: global Heisenberg coupling (optional, default ``0`` when omitted)
- ``K``: global Kitaev coupling (optional, default ``0`` when omitted)
- ``G``: global Gamma coupling (optional, default ``0`` when omitted)
- ``Gp``: global Gamma-prime coupling (optional, default ``0`` when omitted)

For KH runs, you must also provide a bond-type-to-direction mapping table:

.. code-block:: toml

   [Hamiltonian.BondTypeDirection]
   type1 = "x"
   type2 = "y"
   type3 = "z"

Validation rules for ``[Hamiltonian.BondTypeDirection]``:

- ``type1/type2/type3`` are all required for KH model.
- each value must be one of ``x``, ``y``, ``z`` (case-insensitive in parser).
- the three mapped values must be unique.

The runtime classifier groups honeycomb bond geometry into three deterministic bond types
(from PBC-aware neighbor templates) and then maps those three types to ``x/y/z`` using
this table.

.. note::
   KH-only fields (``KitaevEpsilon``, ``KitaevTolerance``, ``J``, ``K``, ``G``, ``Gp``)
   are parsed only for the Kitaev-Heisenberg model. Each value must be numeric
   (integer or floating-point).

Minimal KH input example
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   [Hamiltonian]
   model = "Kitaev-Heisenberg"
   custom = false
   magnetic_field = [ 0.0, 0.0, 0.0 ]
   anisotropy = [ 0.0, 0.0, 0.0 ]
   KitaevEpsilon = 1e-12
   KitaevTolerance = 1e-6
   J = -1.0
   K = 2.0
   G = 0.5
   Gp = -0.25

   [Hamiltonian.BondTypeDirection]
   type1 = "x"
   type2 = "y"
   type3 = "z"

Reference reproducible sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a lightweight reproducible KH sanity run, use:

- ``example/kh_minimal/POSCAR``
- ``example/kh_minimal/input.toml``

This sample intentionally keeps the system and MC schedule minimal while still
covering:

- KH global couplings ``J``, ``K``, ``G``, ``Gp``
- KH direction thresholds ``KitaevEpsilon`` and ``KitaevTolerance``
- required ``[Hamiltonian.BondTypeDirection]`` mapping
- minimal honeycomb geometry classification into three KH bond types

Initialization
~~~~~~~~~~~~~~

- ``angleA``, ``angleB``, ``angleC``: Euler-angle rotations (in degrees) used while propagating spins along supercell directions.
- You can additionally provide ``[[Initialization.Elements]]`` and ``[[Initialization.Elements.Atoms]]`` with explicit ``spin_direction`` for atoms in the first cell.

Output
~~~~~~

- ``magnifying_factor``: scaling for spin structure output.
- ``ground_state``: when ``true``, exports a ground-state structure at ``T = 0`` workflow.

Practical tuning workflow
-------------------------

1. Start from the example values in ``example/Ba2NiIrO6/input.toml``.
2. Confirm neighbor shells first (distances and element pairing).
3. Increase ``relaxing_steps`` and ``counting_steps`` for better statistics.
4. Refine ``temperature_points_number`` near phase-transition regions.
5. Scale ``cell_number`` for finite-size checks.
