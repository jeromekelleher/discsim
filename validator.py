"""
Module used to validate the results of the simulations using various 
means. These are not quite tests, since we don't have exact values 
to check against, and everything is necessarily approximate.
"""
from __future__ import print_function
from __future__ import division 


import sys
import time
import math
import numpy as np
import random
import multiprocessing
from matplotlib import ticker
from matplotlib import pyplot

import ercs
import discsim
import _discsim

class ErcsSingleLocusIdentitySimulator(ercs.Simulator):
    """
    Class that calculates identity in state for genes separated by a range
    of distances.
    """
    def setup(self, num_points, max_distance, mutation_rate, accuracy_goal):
        """
        Sets up the simulation so that we calculate identity at the specified
        number of points, the maximum distance between points is
        max_distance and mutation happens at the specified rate. Also
        set the max_time attribute to reflect the specified accuracy_goal.
        """
        self.mutation_rate = mutation_rate
        self.distances = np.linspace(0, max_distance, num_points)
        self.sample = [None, (0, 0)] + [(0, x) for x in self.distances]
        self.max_time = math.log(accuracy_goal) / (-2 * mutation_rate)

    def get_identity(self, seed):
        """
        Returns the probability of identity at all distance classes
        in this replicate.
        """
        pi, tau = self.run(seed)
        mc = ercs.MRCACalculator(pi[0])
        n = len(self.distances)
        F = [0.0 for j in range(n)]
        for j in range(n):
            mrca = mc.get_mrca(1, j + 2)
            if mrca != 0:
                F[j] = math.exp(-2 * self.mutation_rate * tau[0][mrca])
        return F

class SingleLocusIdentitySimulator(discsim.Simulator):
    """
    Class that calculates identity in state for genes separated by a range
    of distances.
    """
    def __init__(self, torus_diameter, distances, mutation_rate, accuracy_goal):
        super(SingleLocusIdentitySimulator, self).__init__(torus_diameter)
        self.__accuracy_goal = accuracy_goal
        self.__mutation_rate = mutation_rate
        self.__distances = distances
        self.__max_time = math.log(accuracy_goal) / (-2 * mutation_rate)
        self.sample = [None, (0, 0)] + [(0, x) for x in self.__distances]

    def get_identity(self, seed):
        """
        Returns the probability of identity at all distance classes
        in this replicate.
        """
        self.random_seed = seed
        self.run(self.__max_time)
        pi, tau = self.get_history()
        # reset the simulation so we can get another replicate.
        self.reset()
        mc = ercs.MRCACalculator(pi[0])
        n = len(self.__distances)
        F = [0.0 for j in range(n)]
        for j in range(n):
            mrca = mc.get_mrca(1, j + 2)
            if mrca != 0:
                F[j] = math.exp(-2 * self.__mutation_rate * tau[0][mrca])
        return F

def subprocess_identity_worker(t):
    sim, seed = t
    return sim.get_identity(seed)

def run_identity_replicates(sim, num_replicates, worker_pool):
    args = [(sim, random.randint(1, 2**31)) for j in range(num_replicates)]
    replicates = worker_pool.map(subprocess_identity_worker, args)
    mean_identity = np.mean(np.array(replicates), axis=0)
    return mean_identity


def simple_identity_check(r=1, u=0.125, rate=1, num_parents=1, 
        num_replicates=10000, mutation_rate=1e-6):
    """
    Checks identity using very simple model parameters.
    """
    events = [ercs.DiscEventClass(r=r, u=u, rate=rate)]
    ll_events = [e.get_low_level_representation() for e in events]
    torus_diameter = 100
    s = _discsim.IdentitySolver(ll_events,
            torus_diameter=torus_diameter,
            num_quadrature_points=512,
            integration_abserr=1e-6,
            integration_relerr=0,
            integration_workspace_size=1000,
            max_x=50, mutation_rate=mutation_rate,
            num_parents=num_parents)
    s.solve()
    # Set up the simulations
    num_points = 10
    distances = np.linspace(0, 10, num_points)
    sim = SingleLocusIdentitySimulator(torus_diameter, distances, 
            mutation_rate, 1e-6)
    sim.num_parents = num_parents
    sim.event_classes = events
    workers = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    F_sim = run_identity_replicates(sim, num_replicates, workers)
    F_num = [s.interpolate(x) for x in distances]
    for x, fs, fn in zip(distances, F_sim, F_num):
        print("{0:.1f}\t{1:.6f}\t{2:.6f}".format(x, fs, fn)) 
    pyplot.plot(distances, F_sim, label="Simulation")
    pyplot.plot(distances, F_num, label="Numerical")
    pyplot.legend()
    pyplot.show()

    

def mixed_events_identity_check(num_replicates):
    torus_diameter = 100
    num_points = 50
    max_x = 20
    mutation_rate = 1e-6
    accuracy_goal = 1e-3
    small_events = ercs.DiscEventClass(rate=1.0, r=1, u=0.5)
    large_events = ercs.DiscEventClass(rate=0.1, r=10, u=0.05)
    sim = ErcsSingleLocusIdentitySimulator(torus_diameter)
    sim.setup(num_points, max_x, mutation_rate, accuracy_goal)
    workers = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    l = [small_events, large_events]
    sim.event_classes = l
    before = time.time()
    ercs_F = run_identity_replicates(sim, num_replicates, workers)
    duration = time.time() - before
    print("ercs done...", duration)
    distances = np.linspace(0, max_x, num_points)
    sim = SingleLocusIdentitySimulator(torus_diameter, distances, 
            mutation_rate, 1e-6)
    sim.event_classes = l
    before = time.time()
    discsim_F = run_identity_replicates(sim, num_replicates, workers)
    duration = time.time() - before
    print("discsim done...", duration)
    pyplot.plot(distances, ercs_F, label="ercs")
    pyplot.plot(distances, discsim_F, label="discsim")
    pyplot.legend()
    pyplot.show()


def get_mean_squared_displacement(z, pop):
    """
    Returns the mean squared displacement of the specified population from 
    the specified point.
    """
    d2 = 0.0 
    for p, a in pop:
        d2 += (p[0] - z[0])**2
        d2 += (p[1] - z[1])**2
    n = len(pop)
    return d2 / (n * 2) 


def single_locus_diffusion(u, r, rate):
    """
    Measure the mean squared displacement of lineages for a single 
    locus simulation.
    """
    z = (100, 100)
    sample_size = 10000
    s = 2.25
    L = 100 * s
    sim = discsim.Simulator(L)
    sim.pixel_size = s 
    sim.sample = [None] + [z for j in range(sample_size)]
    sim.event_classes = [ercs.DiscEventClass(r=r, u=u, rate=rate)]
    sim.max_occupancy = 2 * sample_size
    sim.max_population_size = 2 * sample_size
    sim.print_state()
    T = []
    X = []
    D = []
    S = []
    for j in range(100):
        t = j * 100 * L**2
        sim.run(t)
        pop = sim.get_population()
        msd = get_mean_squared_displacement(z, pop) 
        t = sim.get_time() / L**2
        T.append(t)
        X.append(msd)
        S.append(t * (r**4) * rate * u * math.pi / 2)
        print(T[-1], X[-1], S[-1])
    pyplot.plot(T, X, T, S)
    pyplot.show()


def subprocess_wave_worker(args):
    sim, times, seed = args
    sim.random_seed = seed
    L = int(sim.torus_diameter)
    n = np.zeros((len(times), L))
    for j, t in enumerate(times):
        sim.run(t)
        pop = sim.get_population()
        for tup in pop:
            if sim.simulate_pedigree:
                k = int(tup)
            else:
                k = int(tup[0])
            n[j, k] += 1
            n[j, k + 1] += 1
    sim.reset()
    return n

def run_wave_replicates(sim, times, num_replicates, worker_pool=None):
    args = [(sim, times, random.randint(1, 2**31)) for j in range(num_replicates)]
    if worker_pool is None:
        replicates = [subprocess_wave_worker(a) for a in args]
    else:
        replicates = worker_pool.map(subprocess_wave_worker, args)
    mean_n = []
    for j in range(len(times)):
        n = []
        for r in replicates:
            n.append(r[j])
        mean_n.append(np.mean(n, axis=0))

    return mean_n

def wave_1d(u, num_loci=0):
    """
    Simulates the wave of pedigree ancestors in 1D.
    """
    N = int(2 / u)
    L = 100
    s = discsim.Simulator(L, num_loci==0)
    if num_loci != 0:
        s.num_loci = num_loci
    s.max_population_size = 10000
    s.event_classes = [ercs.DiscEventClass(r=1, u=u)]
    s.sample = [None, L/2, L/2]
    workers = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    #workers = None
    t = [j * 500 * L for j in range(5)]
    x = [j for j in range(L)]
    for n in run_wave_replicates(s, t, 100, workers):
        pyplot.plot(x, n)
    pyplot.axhline(0.797 * N)
    pyplot.show()
    
def subprocess_structure_worker(args):
    sim, max_time, seed = args
    sim.random_seed = seed
    before = time.time()
    sim.run(max_time)
    sim.print_state()
    sim.reset()
    after = time.time() - before
    return after

def run_structure_replicates(sim, n_reps, max_time, worker_pool=None):
    args = [(sim, max_time, random.randint(1, 2**31)) for j in range(n_reps)]
    if worker_pool is None:
        replicates = [subprocess_structure_worker(a) for a in args]
    else:
        replicates = worker_pool.map(subprocess_structure_worker, args)
    mean_time = np.mean(np.array(replicates), axis=0)
    return mean_time
    
def fine_structure():
    rep_num = 1 # number of simulations over which to average runtimes
    sample_size = 4 # number of sampled lineages
    num_loci = 10 # number of loci per lineage
    arg_effective_population_size = 2 * 10**4 # effective human population size
    # We fix recombination probability by moment matching to human chr 22:
    # Chr 22 has ~ 5e7 loci with recombination probability 1e-8.
    # Thus we expect 0.5 recombinations per reproduction event.
    # rho = 1 / (2 num_loci) gives the same expectation for our locus count.
    rho = 1 / (2 * num_loci)
    r = 1 # event raduius normalised to 1
    gen_no = 10**4 # number of generations to simulate
    
    # Variance of lineage displacement per unit time = u pi r^4 / 2.
    # Time per generation = 1 / (u pi r^2).
    # Combining these yields:
    sigma2 = r**2 / 2 # lineage displacement variance per generation
    
    # Equating mixing time of Brownian motion on torus of side length L
    # and lineage displacement variance over gen_no generations:
    # L^2 log L = sigma2 gen_no
    # yields torus side length with reasonable wrap around probability.
    # Ignore log L correction for simplicity and round down to an even integer.
    L = math.sqrt(gen_no * sigma2)
    L = L - L % 2
    
    # 2 * Ne * pi * sigma2 = Wrights neighbouhood size = 2 / u
    # Thus the migration rate is a = 2 / u log L (c.f. Theorem 5.13 of
    # Durrett, Probability models for DNA sequence evolution).
    # By Theorem 5.12 of Durrett, Ne = (1 + a) L^2 log L / (2 pi sigma2).
    # Substituting in for a and solving for u yields
    u = 2 * L**2 / (2 * arg_effective_population_size * math.pi * sigma2 
        - L**2 * math.log(L))
    gen_time = 1 / (math.pi * r**2 * u) # time per generation
    # total time for gen_no generations scaled for total event rate 1.
    max_time = gen_no * gen_time * L**2
    s = discsim.Simulator(L, False, True)
    s.num_loci = num_loci
    s.recombination_probability = rho
    #ARG recombination rate moment matched to rho
    s.arg_recombination_rate = 1 / (num_loci * arg_effective_population_size) 
    s.arg_effective_population_size = arg_effective_population_size
    s.max_population_size = sample_size * num_loci + 2
    s.sample = [None] + [(random.random()*L, random.random()*L) for j in range(sample_size)]
    s.event_classes = [ercs.DiscEventClass(r=1, u=u)]
    #print(run_structure_replicates(s, rep_num, max_time, None))
    #workers = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    #print(run_structure_replicates(s, rep_num, max_time, workers))
    s.run(max_time)
    s.print_state()
    pop = s.get_population()
    for entry in pop:
        print(entry[0][0], entry[0][1], entry[1])
    return 1

def main():
    #simple_identity_check(rate=0.5)
    #simple_identity_check(r=0.93, u=0.133, rate=0.5, num_parents=2, 
    #        num_replicates=10**6, mutation_rate=1e-7)
    #mixed_events_identity_check(100000)
    #plot_mixed_events_identity()
    #single_locus_diffusion(u=0.0000125, r=1, rate=1.0)
    #wave_1d(u=0.005)
    #wave_1d(u=0.005, num_loci=100000)
    fine_structure()

if __name__ == "__main__":
    main()
