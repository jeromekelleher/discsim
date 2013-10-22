===================================================
Efficient coalescent simulation in continuous space
===================================================

Simulates the coalescent for populations evolving in a spatial 
continuum under the extinction/recolonisation model. This 
package is a specialisation of the `ercs <https://pypi.python.org/pypi/ercs>`_
package, and provides a much more efficient method of 
simulating the spatial coalescent for the disc model. 
A very similar (but not identical) interface to ``ercs``
is provided.


The simulations 
support:
        
- A sample of ``n`` individuals with ``m`` loci at arbitrary locations on a 
  torus of diameter ``L``.
- A fixed recombination rate between pairs of adjacent loci. 
- An arbitrary number of classes of event occuring at fixed
  rates. 
- Simulations in one and two dimensions.
- Access to the locations of ancestors at any time in the past.
- Simulations of the locations of all pedigree ancestors, as well 
  as the genetic ancestors.

The ``discsim`` module supports Python 2 and 3.

-------------
Documentation
-------------

Here's a quick example for the impatient::

        import ercs
        import discsim
        sim = discsim.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.run()
        pi, tau = sim.get_history()


Full documentation for ``discsim`` is available at 
`<http://pythonhosted.org/discsim>`_.

------------
Installation
------------

*******************************
Quick install for Debian/Ubuntu
*******************************

If you are running Debian or Ubuntu, this should get you up and running quickly::

        $ sudo apt-get install python-dev libgsl0-dev 
        $ sudo pip install ercs discsim

For Python 3, use ``python3-dev`` and ``pip3``.

********************
General instructions
********************

The ``discsim`` module 
depends on the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_,
which must be installed before it can be built.
Fortunately, this is straightforward on most platforms. For example, 
on Debian or Ubuntu use::

        $ sudo apt-get install libgsl0-dev

or on Fedora::

        $ sudo yum install gsl-devel

GSL is available on most packaging systems; if it is not available on your
platform, it can be installed from source.

The ``discsim`` module also depends on the 
`ercs <https://pypi.python.org/pypi/ercs>`_
Python module, which 
must also be installed, using the same methods as outlined below.

Once GSL has been installed we can build the ``discsim`` module using the 
standard Python `methods <http://docs.python.org/install/index.html>`_. For 
example, using pip we have ::
        
        $ sudo pip install discsim

Or, we can manually download the package, unpack it and then run::
        
        $ python setup.py build
        $ sudo python setup.py install

Most of the time this will compile and install the module without difficulty.

It is also possible to download the latest development version of 
``discsim`` from `github <https://github.com/jeromekelleher/discsim>`_. 

******************
Potential problems
******************

On platforms that GSL is not available as part of the native packaging 
system (or GSL was installed locally because of non-root access)
there can be issues with finding the correct headers and libraries
when compiling ``ercs`` and ``discsim``. For example, on FreeBSD we get something 
like this::

        $ python setup.py build
        ... [Messages cut for brevity] ...
        _discsimmodule.c:515: error: 'sim_t' has no member named 'time'
        _discsimmodule.c: In function 'Simulator_get_num_reproduction_events':
        _discsimmodule.c:529: error: 'sim_t' has no member named 'num_reproduction_events'
        _discsimmodule.c: In function 'Simulator_get_history':
        _discsimmodule.c:743: error: 'sim_t' has no member named 'pi'
        _discsimmodule.c:748: error: 'sim_t' has no member named 'tau'
        _discsimmodule.c: In function 'Simulator_run':
        _discsimmodule.c:789: error: 'sim_t' has no member named 'time'
        error: command 'cc' failed with exit status 1
        
This can be remedied by using the ``gsl-config`` program to set the 
the ``LDFLAGS`` and ``CFLAGS`` environment variables to 
their correct values::
        
         $ CFLAGS=`gsl-config --cflags` LDFLAGS=`gsl-config --libs` python setup.py build

*****
Tests
*****

``discsim`` provides some test cases to ensure that the installation has gone smoothly.
It is a good idea to run these immediately after installation::

        $ python tests.py


****************
Tested platforms
****************

Discsim has been successfully built and tested on the following platforms:

================        ========        ======          ========
Operating system        Platform        Python          Compiler
================        ========        ======          ========
Debian wheezy           x86_64          2.7.3           gcc 4.7.2	
Debian wheezy           x86_64          3.2.3           gcc 4.7.2	
Debian wheezy           x86             2.7.3           gcc 4.7.2	
Debian squeeze          ppc64           2.6.6           gcc 4.4.5	
Debian squeeze          ppc64           3.1.3           gcc 4.4.5	
Debian squeeze          x86_64          2.6.6           gcc 4.4.5	
Debian squeeze          x86_64          3.1.3           gcc 4.4.5	
FreeBSD 9.2             x86_64          2.7.5           gcc 4.2.1
FreeBSD 9.2             x86_64          3.3.2           gcc 4.2.1
Fedora 19               x86_64          2.7.5           gcc 4.8.1
Fedora 19               x86_64          3.3.2           gcc 4.8.1
SunOS 5.10              SPARC           3.3.2           gcc 4.8.0 
================        ========        ======          ========


