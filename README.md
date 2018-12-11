An ambiguity problem arises then determining the position and motion of objects with a radar system. The ambiguity problem
is translated in the fact that different Direction of Arrival (DOA) can lead to the same response. In  the paper
by Daniel Kastinen [Determining all ambiguities in direction of arrival measured by radar systems](http://www.ursi.org/content/RSB/RSB_365_2018_06.pdf), 
a mathematical framework and practical method to find all ambiguities in any multichannel system is described. 
The new formulation allows for an efficient implementation using the numerical Moore-Penrose inverse to find all
ambiguities and approximate ambiguities.

Here two main functions composed by other several functions are developed to implement the algorithm proposed and
retrieve solutions. The *ambiguity_calculator.py* file contains a function that by entering the name of a listed radar
configuration, described in *radarconf.py* and operating frequency is able to calculate all ambiguities of a radar
system and saves the results into an HDF5 file for later processing. *plots_generator.py* contains a function that
generates the solution and plots for a certain radar configuration, with a certain operating frequency and for a given
signal wave (described as a wave vector). All of the use user defined functions in *functions.py*.

By running a python script with


    from ambiguity_calculator import ambiguities_calculate
    
    ambiguities_calculate(radar_name='Ydist', frequency=31)
    
a HDF5 called JONES.h5 containing the calculation results is generated in the folder ../processed_data/JONES.

JONES.h5 contains several items, data sets. Organized between two main HDF5 groups.

Continuing with JONES radar configuration, if one runs another script with


    from plots_generator import generate_plots
    
    generate_plots(radar_name='JONES', frequency=31, elevation=50, azimuth=270)

the same HDF5 file will the imported and plots with the results will be generated.

Written and tested in Python 3.7.