easistrain |release|
====================

A Python library for strain/stress analysis from X-ray diffraction data.

Getting started
---------------

Install the library

.. code:: bash

    python -m pip install git+https://github.com/EASI-STRESS/easistrain

Process energy-dispersive X-ray diffraction data wuth

.. code:: bash

    python -m easistrain.EDD.<task_name> <config_file>

`<task_name>` should be replaced by a valid task name (see below).

`<config_file>` should be replaced by the path to the configuration for the task. Examples of config files for each task can be found at https://github.com/EASI-STRESS/easistrain/tree/main/example_config or in the `example_config` folder of the source code.

The available tasks are

- Task 1 → `calibrationEDD` → calibration of the conversion from Channel → Energy
- Task 2 → `angleCalibEDD` → calibration of the diffraction angle
- Task 3 → `fitEDD` → fitting
- Task 4 → `coordTransformation` → transformation of the coordinates from the gonio reference to the sample reference
- Task 5 → `regroupPoints` → regroup all the points
- Task 6 → `preStraind0cstEDD` → calculates the strain in the measurement direction
- Task 7 → `strainStressd0cstEDD` → calculates the strain and stress tensors

Documentation
-------------

.. toctree::
    :maxdepth: 2

    id15_NXstress
    theory
    reference_frames
    notes
    api
