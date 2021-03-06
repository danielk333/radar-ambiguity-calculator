Example of utilization and expected outputs.
============================================

Solve the problem for a certain radar configuration
###################################################

As an example take one of the radar configurations, **JONES** in this case, with a frequency of 31 MHz.

The coordinates of the subarray in wave lengths are

==== ====
x    y
==== ====
0    2
0    -2.5
-2   0
2.5  0
0    0
==== ====

By running a python script with

::

    from ambiguity_calculator import ambiguities_calculate

    ambiguities_calculate(radar_name='Ydist', frequency=31)

a HDF5 called JONES.h5 containing the calculation results is generated in the folder ../processed_data/JONES.

JONES.h5 contains several items, data sets. Organized between two main HDF5 groups.

* root:

    * trivial_calculations: holds results from calculations which are straight forward and whose results can be used for to track the results.

        * *sensor_groups*
        * *subgroup_phase_center*
        * *linear_coefficients*
        * *base_numbers*
        * *k_length*

    * results_permutations: holds the results of the permutations and ambiguities.

        * *intersections_integers_complete*
        * *ambiguity_distances_INT_FORM_MAT*
        * *ambiguity_distances_INT_FORM_mean*
        * *ambiguity_distances_WAVE_FORM_MAT*
        * *ambiguity_distances_WAVE_FORM*
        * *intersection_line*
        * *survivors*

Use the results to see the ambiguities for a DOA.
##################################################

Continuing with JONES radar configuration, if one runs another script with

::

    from plots_generator import generate_plots

    generate_plots(radar_name='JONES', frequency=31, elevation=50, azimuth=270)

the same HDF5 file will the imported and plots with the results will be generated.

.. figure:: figures/figure1.png
    :scale: 80%
    :align: center

.. figure:: figures/figure2.png
    :scale: 80%
    :align: center

.. figure:: figures/figure3.png
    :scale: 80%
    :align: center

.. figure:: figures/figure4.png
    :scale: 80%
    :align: center

.. figure:: figures/figure6.png
    :scale: 80%
    :align: center

.. figure:: figures/figure7.png
    :scale: 80%
    :align: center


