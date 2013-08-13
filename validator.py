"""
Module used to validate the results of the simulations using various 
means. These are not quiet tests, since we don't have exact values 
to check against, and everything is necessarily approximate.
"""
from __future__ import print_function
from __future__ import division 


import sys
import ercs
import math
import multiprocessing
import numpy as np
import random
import _discsim

class SingleLocusIdentitySimulator(object):
    """
    Class that calculates identity in state for genes separated by a range
    of distances.
    """
    def __init__(self, distances, mutation_rate, accuracy_goal, events, **kwargs):
        self.__params = kwargs 
        self.__events = events
        self.__mutation_rate = mutation_rate
        self.__distances = distances
        self.__sample = [(0, 0)] + [(0, x) for x in self.__distances]
        self.__max_time = math.log(accuracy_goal) / (-2 * mutation_rate)

    def run(self, seed):
        """
        Runs the simulation and returns the simulated history.
        """
        self.__params["max_time"] = self.__max_time
        self.__params["random_seed"] = seed 
        s = _discsim.Simulator(self.__sample, self.__events, **self.__params)
        s.run()
        return s.get_history()

    def get_identity(self, seed):
        """
        Returns the probability of identity at all distance classes
        in this replicate.
        """
        pi, tau = self.run(seed)
        mc = ercs.MRCACalculator(pi[0])
        n = len(self.__distances)
        F = [0.0 for j in range(n)]
        for j in range(n):
            mrca = mc.get_mrca(1, j + 2)
            if mrca != 0:
                F[j] = math.exp(-2 * self.__mutation_rate * tau[0][mrca])
        return F

def subprocess_worker(t):
    sim, seed = t
    return sim.get_identity(seed)

def run_replicates(sim, num_replicates, worker_pool):
    args = [(sim, random.randint(1, 2**31)) for j in range(num_replicates)]
    replicates = worker_pool.map(subprocess_worker, args)
    mean_identity = np.mean(np.array(replicates), axis=0)
    return mean_identity


def simple_identity_check(r=1, u=0.125, num_parents=1, num_replicates=10000, 
        mutation_rate=1e-6):
    """
    Checks identity using very simple model parameters.
    """
    events = [{"r":r, "u":u, "rate":1}]
    torus_diameter = 100
    s = _discsim.IdentitySolver(events, 
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
    sim = SingleLocusIdentitySimulator(distances, mutation_rate, 1e-6,
            events, torus_diameter=torus_diameter,
            num_parents=num_parents)
    workers = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    identity = run_replicates(sim, num_replicates, workers)
    for x, f in zip(distances, identity):
        print("{0:.1f}\t{1:.6f}\t{2:.6f}".format(x, f, s.interpolate(x)))

def main():
    #simple_identity_check()
    simple_identity_check(0.73, 0.0133, 3, 10**4, 1e-9)

if __name__ == "__main__":
    main()
