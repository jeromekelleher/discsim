.. discsim documentation master file, created by
   sphinx-quickstart on Wed Oct 16 16:16:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to discsim's documentation!
===================================

Contents:

.. toctree::
   :maxdepth: 2

-----------------------------------
:mod:`discsim` -- Module reference
-----------------------------------

.. module:: discsim

The ``discsim`` module provides a convenient interface to coalescent simulations 
for the extinction/recolonisation model, and some utilities to process 
the resulting genealogies.


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
        :class:`ercs.DiscEventClass`
        instances. There must be at least one event class specified.
        
        **Default value:** ``None``.

    .. attribute:: torus_diameter

        The diameter of the torus we are simulating on. This defines the 
        size of the space that lineages can move around in. 
        
        **Default value:** Specified at instantiation time.

    .. attribute:: num_parents 

       The number of parents in each event. For a single locus simulation 
       there must be at least one parent and for multi-locus simulations 
       at least two. 
        
       **Default value:** 1 if the simulation is single locus, otherwise 2.

    .. attribute:: recombination_probabilities
    	
       The list of inter-locus recombination probabilities; the length of 
       this list also determines the number of loci for each individual.
       At an event, the probability of locus ``j`` and ``j + 1`` descending 
       from different parents is ``recombination_probabilities[j]``.
       The number of loci in the simulation is therefore 
       ``len(recombination_probabilities) + 1``.
       
       **Default value:** The empty list [] (so, we have a single locus 
       simulation by default).

    .. attribute:: max_time 

       The maximum amount of time (in simulation units) that we simulate. If 
       this is set to `0.0` the simulation continues until all loci have 
       coalesced.
       
       **Default value:** 0.0

    .. attribute:: max_lineages

       The maximum number of extant lineages in the simulation. 
       If the number of lineages we are tracking exceeds this limit, 
       the simulation aborts and raises an :exc:`_discsim.LibraryError`.  
       
       **Default value:** 1000 

    .. automethod:: discsim.Simulator.run



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

