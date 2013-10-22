.. discsim documentation master file, created by
   sphinx-quickstart on Wed Oct 16 16:16:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================
Documentation for discsim
=========================

.. toctree::
   :maxdepth: 2

.. contents::

This is the documentation for the :mod:`discsim` module, an efficient coalescent 
simulator for the disc-based extinction/recolonisation continuum model. This 
module is intended as a companion to the :mod:`ercs` module, which provides 
support for simulating the history of a sample under more general models. 
The :mod:`ercs` module, however, is much less efficient for large neighbourhood
sizes, and the :mod:`discsim` module documented here uses a different 
underlying simulation algorithm. The interface presented by the two modules 
is very similar (although not identical), so those unfamiliar with the 
interface are recommended to first read the documentation for the 
:mod:`ercs` module. In particular, the representation of the history of 
a sample uses the same format, and the :class:`ercs.MRCACalculator` class 
can be used to find most recent common ancestors for both simulations.

Besides the specialisation of the simulation to the disc model, there are
some other differences between the simulations. Firstly, :mod:`discsim`
provides access to the state of the sample during simulations, allowing 
us to simulate the distribution of the location of ancestors. We also 
support simulating the history of a sample in one or two dimensions, 
which is achieved by simply providing a 1D or 2D sample of locations.
Finally, :mod:`discsim` allows us to simulate the history of either 
the genetic ancestors and their history (as before) or the pedigree 
ancestors of the sample.


-----------------------------------
:mod:`discsim` -- Module reference
-----------------------------------

.. module:: discsim


***************************
:class:`discsim.Simulator`
***************************
.. autoclass:: discsim.Simulator

    .. attribute:: sample
       
        The location of individuals at the beginning of the simulation. This 
        must be either a list of 2-tuples or numbers describing locations 
        within the space defined by the torus. If a list of numbers is 
        provided, the simulation is performed in a 1D environment; if a 
        list of tuples is provided, the simulation occurs in 2D. Dimensions 
        cannot be mixed.  The zero'th element of the list **must be** ``None``.

        **Default value:** ``None``.

    .. attribute:: event_classes
        
        The event classes to simulate. This must be a list of 
        :class:`ercs.DiscEventClass` instances. There must be at least one event 
        class specified. The underlying algorithm will take the first event class 
        and only simulate events in which there is a high 
        probability of a lineage jumping. Subsequent event classes are not treated
        in any special way, and are simulated directly without conditioning.
        
        **Default value:** ``None``.

    .. attribute:: torus_diameter

        The diameter of the torus we are simulating on. This defines the 
        size of the 1D or 2D space that lineages can move around in. 
        
        **Default value:** Specified at instantiation time.

    .. attribute:: num_parents 

       The number of parents in each event. For a single locus simulation 
       there must be at least one parent and for multi-locus simulations 
       at least two. 
        
       **Default value:** 1 if the simulation is single locus genetic simulation;
       otherwise 2.

    .. attribute:: recombination_probability
    	
       The probability of recombination between adjacent loci at an event.

       **Default value:** 0.5 (free recombination).
    
    .. attribute:: num_loci 
    	
       The number of loci we simulate the history of in a genetic simulation. 

       **Default value:** 1. 

    .. attribute:: pixel_size

       The length of one edge of a pixel. This is an important performance tuning 
       factor. For a single locus simulation, this should be around 2.25 for best 
       performance; for multilocus (or pedigree) simulations, this should be less 
       than 2.25. For pixel_size less than 1, memory requirements increase sharply
       and may not justify any increase in speed of simulation.

       For 1D simulations, this **must** be equal to 2.
   
       There is a strong requirement that ``torus_diameter`` / ``pixel_size``
       must be an integer **exactly**. Thus, there can be issues with floating 
       point division not providing exact values. The best solution to this 
       is to use pixel sizes that are exactly representable in binary (such as 
       2.125, 1.5, etc).

       **Default value:** 2

    .. attribute:: random_seed

       The random seed for the current simulation run.

       **Default value:** A random value chosen by the standard Python random number 
       generator.

    .. attribute:: max_population_size

       The maximum number of extant individuals in the simulation. 
       If the number of individuals we are tracking exceeds this limit, 
       the simulation aborts and raises an :exc:`_discsim.LibraryError`.  
       
       **Default value:** 1000 
    
    .. attribute:: max_occupancy

       The maximum number of individuals that may occupy a simulation 
       pixel. This is defined as the number of individuals that are within
       distance r of the pixel itself (**not** just those within the pixel).
       If this value is exceeded the simulation aborts and raises 
       an :exc:`_discsim.LibraryError`.  
       
       **Default value:** N times the catchment area of a pixel, where N 
       is the neighbourhood size. This should be sufficient for most purposes,
       but may need to be increased for very high sampling densities. 

    .. automethod:: discsim.Simulator.run

    .. automethod:: discsim.Simulator.get_population

    .. automethod:: discsim.Simulator.get_history

    .. automethod:: discsim.Simulator.reset

    .. automethod:: discsim.Simulator.get_time

    .. automethod:: discsim.Simulator.get_num_reproduction_events

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

