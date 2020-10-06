import mdtraj as md
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.optimize
import unyt as u

class BondCalculator:
    def __init__(self, traj, T):
        self.traj = traj
        graph = traj.top.to_bondgraph()
        bonds = self.identify_bonds(graph)
        angles = self.identify_angles(graph)

        bond_params = dict()
        angle_params = dict()

        for bond_type, pairs in bonds.items():
            bond_lengths, bond_prob = self.calc_lengths(pairs, range=[0, 1.0])
            params = self.calc_parameters(bond_lengths, bond_prob)
            k = 2 * u.kb * (T*u.K) / (params[0] * u.nm)**2 * u.Na
            l0 = params[1] * u.nm
            bond_params[bond_type]= {"k": k, "x0": l0}

        for angle_type, triplets in angles.items():
            bond_angles, angle_prob = self.calc_angles(triplets, range=[0, 2*np.pi])
            params = self.calc_parameters(bond_angles, angle_prob)
            k = 2 * u.kb * (T*u.K) / (params[0] * u.rad)**2 * u.Na
            t0 = params[1] * u.rad
            angle_params[angle_type]= {"k": k, "x0": t0}

        self.bond_params = bond_params
        self.angle_params = angle_params

    def identify_bonds(self, graph):
        all_bonds = [edge for edge in graph.edges]
        bonds = defaultdict(list)
        for bond in all_bonds:
            index = tuple(sorted([bond[0].name, bond[1].name]))
            pair = tuple([particle.index for particle in bond])
            bonds[index].append(pair)

        return bonds

    def identify_angles(self, graph):
        angle_subgraph = nx.Graph()
        angle_subgraph.add_edge(0, 1)
        angle_subgraph.add_edge(1, 2)

        matcher = nx.algorithms.isomorphism.GraphMatcher(graph, angle_subgraph)

        all_angles = []

        for m in matcher.subgraph_isomorphisms_iter():
            all_angles.append(tuple(k for k in m.keys()))

        angles = defaultdict(list)
        for angle in all_angles:
            index = tuple(particle.name for particle in angle)
            if angle[0].name < angle[2].name:
                index = tuple(reversed(index))

            triplet = tuple(particle.index for particle in angle)
            angles[index].append(triplet)

        return angles

    def calc_lengths(self, pairs, range=None):
        quantity = md.compute_distances(self.traj, pairs)
        hist, edges = np.histogram(quantity, density=True, range=range, bins=200)
        bins = (edges[1:]+edges[:-1]) * 0.5
        return bins, hist

    def calc_angles(self, triplets, range=None):
        quantity = md.compute_angles(self.traj, triplets)
        hist, edges = np.histogram(quantity, density=True, range=range, bins=200)
        bins = (edges[1:]+edges[:-1]) * 0.5
        hist /= np.sin(bins)
        hist /= np.sum(hist)*(bins[1]-bins[0])
        return bins, hist

    def cost_function(self, args, x, y):
        w, x0 = args
        return np.sum((self.gaussian(w, x0, x) - y)**2)

    def gaussian(self, w, x0, x):
        return ((w * np.sqrt(np.pi / 2))**-1)*(np.exp(-2 * (x - x0)**2 / (w**2)))

    def calc_parameters(self, x, y):
        res = scipy.optimize.minimize(lambda args: self.cost_function(args, x, y), [np.ptp(x)/10, x[np.argmax(y)]])
        return res.x
