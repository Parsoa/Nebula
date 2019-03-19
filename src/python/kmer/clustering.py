from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import time
import argparse
import operator
import traceback
import subprocess

from kmer import (
    bed,
    config,
    gapped,
    counter,
    reduction,
    simulator,
    counttable,
    map_reduce,
    production,
    statistics,
    visualizer,
    programming,
)

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

import numpy as np
np.set_printoptions(threshold = np.inf)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class RemnantCounterJob(production.MixCounterJob, map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'RemnantCounterJob'
    _category = 'programming'
    _previous_job = production.MixCounterJob
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = RemnantCounterJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.half_mers = {}
        self.inner_kmers = {}
        self.gapped_kmers = {}
        self.load_inner_kmers()
        self.load_gapped_kmers()
        self.round_robin()

    def load_inner_kmers(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'inner_kmers.json'), 'r') as json_file:
            self.inner_kmers = json.load(json_file)
        print('Counting', green(len(self.inner_kmers)), 'inner kmers')
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            _kmers = {}
            for kmer in self.inner_kmers:
                _kmers[kmer] = {}
                _kmers[kmer]['loci'] = {}
                _kmers[kmer]['tracks'] = self.inner_kmers[kmer]['tracks']
                for locus in self.inner_kmers[kmer]['loci']:
                    _kmers[kmer]['loci'][locus] = {}
                    _kmers[kmer]['loci'][locus]['masks'] = self.inner_kmers[kmer]['loci'][locus]['masks']
            json.dump(_kmers, json_file, indent = 4)

    def load_gapped_kmers(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'gapped_kmers.json'), 'r') as json_file:
            self.gapped_kmers = json.load(json_file)
            for kmer in self.gapped_kmers:
                left = kmer[:c.hsize]
                right = kmer[-c.hsize:]
                if not left in self.half_mers:
                    self.half_mers[left] = {}
                self.half_mers[left][right] = kmer 
                left = reverse_complement(left)
                right = reverse_complement(right)
                if not right in self.half_mers:
                    self.half_mers[right] = {}
                self.half_mers[right][left] = kmer
        print('Counting', green(len(self.gapped_kmers)), 'gapped kmers')
        with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
            json.dump(self.half_mers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'pre_gapped_kmers.json'), 'w') as json_file:
            json.dump(self.gapped_kmers, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ClusterIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'ClusterIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = RemnantCounterJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = ClusterIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def transform(self, track, track_name):
        c = config.Configuration()
        print(green(track_name))
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            _lp_kmers = {}
            inner_kmers = kmers['inner_kmers']
            _kmers = {kmer: inner_kmers[kmer] for kmer in list(filter(lambda kmer: len(inner_kmers[kmer]['tracks']) == 1, inner_kmers))}
            unique_kmers = list(filter(lambda kmer: len(_kmers[kmer]['loci']) == 1, _kmers))
            non_unique_kmers = list(filter(lambda kmer: len(_kmers[kmer]['loci']) != 1, _kmers))
            if len(unique_kmers) >= 5:
                inner_kmers = {kmer: _kmers[kmer] for kmer in unique_kmers[:min(50, len(unique_kmers))]}
            else:
                inner_kmers = {kmer: _kmers[kmer] for kmer in unique_kmers}
                for kmer in non_unique_kmers:
                    inner_kmers[kmer] = _kmers[kmer]
                    if len(inner_kmers) > 50:
                        break
            for kmer in inner_kmers:
                _lp_kmers[kmer] = True
                self.lp_kmers[kmer] = {} 
                self.lp_kmers[kmer]['type'] = 'inner'
                self.lp_kmers[kmer]['count'] = inner_kmers[kmer]['count']
                self.lp_kmers[kmer]['doubt'] = inner_kmers[kmer]['doubt']
                self.lp_kmers[kmer]['total'] = inner_kmers[kmer]['total']
                self.lp_kmers[kmer]['weight'] = 1.0
                self.lp_kmers[kmer]['coverage'] = c.coverage#inner_kmers[kmer]['coverage']
                self.lp_kmers[kmer]['reference'] = len(inner_kmers[kmer]['loci'])
                self.lp_kmers[kmer]['reduction'] = inner_kmers[kmer]['reference']
                self.lp_kmers[kmer]['tracks'] = inner_kmers[kmer]['tracks']
                self.lp_kmers[kmer]['loci'] = list(map(lambda l: l, inner_kmers[kmer]['loci']))
            gapped_kmers = kmers['gapped_kmers']
            for kmer in gapped_kmers:
                _lp_kmers[kmer] = True
                self.lp_kmers[kmer] = {
                    'gap': gapped_kmers[kmer]['gap'],
                    'side': 'outer',
                    'type': 'gapped', 
                    'count': gapped_kmers[kmer]['count'],
                    'tracks': {},
                    'weight': 1.0,
                    'reference': gapped_kmers[kmer]['tracks'][track_name],
                }
                self.lp_kmers[kmer]['tracks'][track_name] = 1
        path = os.path.join(self.get_current_job_directory(), 'lp_kmers_' + track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in _lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['coverage'] = c.coverage# if kmer['type'] == 'gapped' else kmer['coverage']
            kmer['residue'] = 0 if kmer['type'] == 'gapped' else kmer['reference'] - r
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
            l = len(list(filter(lambda i: self.lp_kmers[i]['type'] == 'inner', self.tracks[kmer['tracks'].keys()[0]]['kmers'])))
            l = l if l else 1
            kmer['weight'] = 2 * l if kmer['type'] == 'gapped' else 1.0

    def generate_linear_program(self):
        c = config.Configuration()
        globals()['cplex'] = __import__('cplex')
        problem = cplex.Cplex()
        problem.objective.set_sense(problem.objective.sense.minimize)
        # the coverage of each event
        names = [''] * len(self.tracks)
        for track in self.tracks:
            tokens = track.split('_')
            names[self.tracks[track]['index']] = 'c' + tokens[1]
        problem.variables.add(names = names,
            ub = [1.0] * len(self.tracks),
        )
        # the real-valued error parameter for kmer
        problem.variables.add(names = ['e' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            lb = [-100000000 for kmer in self.lp_kmers]
        )
        # absolute value of the kmer error parameter
        problem.variables.add(names = ['l' + str(index) for index, kmer in enumerate(self.lp_kmers)],
            obj = [kmer['weight'] for index, kmer in enumerate(self.lp_kmers)]
        )
        #self.add_snp_linear_constraints(problem)
        for index, kmer in enumerate(self.lp_kmers):
            self.add_error_absolute_value_constraints(problem, index)
            if kmer['type'] == 'inner':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                ind.append(len(self.tracks) + index) # Objective
                #ind.append(len(self.tracks) + 2 * len(self.lp_kmers) + index)
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                val.append(1.0) #Objective
                #val.append(kmer['count'] - kmer['coverage'] * kmer['residue'])
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
                    senses = ['E']
                )
            if kmer['type'] == 'gapped' and kmer['side'] == 'outer':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                #ind.append(len(self.tracks) + 2 * len(self.lp_kmers) + index)
                val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                #val.append(kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))))
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['E']
                )
        return problem

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ClusterIntegerProgrammingAnalysisJob(programming.IntegerProgrammingJob):

    _name = 'ClusterIntegerProgrammingAnalysisJob'
    _category = 'programming'
    _previous_job = RemnantCounterJob 

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = ClusterIntegerProgrammingAnalysisJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        self.tracks = self.load_previous_job_results()
        print(self.tracks)
        paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/simulation/HG00513.hg19.chr17.DEL.bed/Seed' + str(i) + 'SnpError0.001/20x/ClusterIntegerProgrammingJob' for i in range(2001, 3000)]
        for track in self.tracks:
            self.tracks[track] = {}
            for r in ['00', '10', '11']:
                for p in ['00', '10', '11']:
                    self.tracks[track][r + '_as_' + p] = 0
            print(self.tracks[track])
        for path in paths:
            for track in self.tracks:
                self.tracks[track][path] = {}
            self.verify_genotypes(path)
        for track in self.tracks:
            x = []
            y = []
            print(track)
            x.append('00')
            y.append(1.0)
            x.append('10')
            y.append(0.5)
            x.append('11')
            y.append(0.0)
            x.append('00')
            y.append(1.0 - 0.01)
            x.append('10')
            y.append(0.5 + 0.01)
            x.append('11')
            y.append(0.0 + 0.01)
            for path in paths:
                #if self.tracks[track][path]['lp_genotype'] == self.tracks[track][path]['actual_genotype']:
                x.append(self.tracks[track][path]['actual_genotype'])
                y.append(self.tracks[track][path]['lp_value'])
                lp_genotype = self.tracks[track][path]['lp_genotype']
                actual_genotype = self.tracks[track][path]['actual_genotype']
                self.tracks[track][actual_genotype + '_as_' + lp_genotype] += 1
            #for p in paths:
            #    for q in paths:
            #        if p != q:
            #            if self.tracks[track][p]['actual_genotype'] > self.tracks[track][q]['actual_genotype'] and self.tracks[track][p]['lp_value'] < self.tracks[track][q]['lp_value']
            print(len(x))
            print(len(y))
            #visualizer.violin(x, y, 'cluster_' + track, self.get_current_job_directory(), 'genotype', 'lp value')
            #visualizer.histogram(y, 'cluster_' + track, self.get_current_job_directory(), 'genotype', 'lp value', step = 0.05)
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4, sort_keys = True)
        exit()

    def verify_genotypes(self, path):
        print(path)
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        #output = subprocess.call('verify_genotypes', shell = True, stderr = subprocess.STDOUT, cwd = path)
        for r in ['00', '10', '11']:
            for p in ['00', '10', '11']:
                for track in bed.load_tracks_from_file(os.path.join(path, r + '_as_' + p + '.bed'), [('lp_genotype', None, str), ('lp_value', None, float)]):
                    self.tracks[str(track)][path]['lp_value'] = track.lp_value
                    self.tracks[str(track)][path]['lp_genotype'] = p
                    self.tracks[str(track)][path]['actual_genotype'] = r

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class KMeansClusteringJob(map_reduce.Job):

    _name = 'KMeansClusteringJob'
    _category = 'programming'
    _previous_job = RemnantCounterJob 

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = KMeansClusteringJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        self.tracks = self.load_previous_job_results()
        print(self.tracks)
        self.paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/simulation/HG00513.hg19.chr17.DEL.bed/Seed' + str(i) + 'SnpError0.001/20x/ClusterIntegerProgrammingJob' for i in range(2001, 3000)]
        for track in self.tracks:
            self.tracks[track] = {}
            for r in ['00', '10', '11']:
                for p in ['00', '10', '11']:
                    self.tracks[track][r + '_as_' + p] = 0
            print(self.tracks[track])
        for path in self.paths:
            self.verify_genotypes(path)
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        c = config.Configuration()
        features = []
        _kmers = {}
        for index, path in enumerate(self.paths):
            features.append([track[path]['lp_value']])
        features = np.array(features)
        genotypes, centroids = self.kmeans(track_name, features)
        for index, path in enumerate(self.paths):
            track[path]['kmeans_genotype'] = genotypes[index]
            #track[path]['features'] = features[index]
            actual_genotype = track[path]['actual_genotype']
            track[actual_genotype + '_as_' + genotypes[index]] += 1
        for r in ['00', '10', '11']:
            for p in ['00', '10', '11']:
                print(yellow(r + '_as_' + p), ':', track[r + '_as_' + p])
        self.calculate_lp_false_rate(track_name)
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(track, json_file, indent = 4, sort_keys = True)
        return path

    #def transform(self, track, track_name):
    #    if track_name != 'chr17_10886864_10895981':
    #        return None
    #    print(green(track_name))
    #    c = config.Configuration()
    #    features = []
    #    _kmers = {}
    #    for index, path in enumerate(self.paths):
    #        with open(os.path.join(path, track), 'r') as json_file:
    #            if index == 0:
    #                _kmers = json.load(json_file)
    #                kmers = _kmers
    #            else:
    #                kmers = json.load(json_file)
    #            features.append([min(kmers[kmer]['count'] / float(c.coverage), 1.0) for kmer in _kmers if _kmers[kmer]['type'] == 'inner'])
    #            #print(len(_kmers), len(features[0]))
    #            if len(features[0]) == 0:
    #                return None
    #    features = np.array(features)
    #    genotypes, centroids = self.kmeans(features)
    #    path = os.path.join(self.get_current_job_directory(), track_name + '.json')
    #    with open(path, 'w') as json_file:
    #        json.dump({'genotypes': genotypes, 'centroids': centroids.tolist(), 'features': features.tolist()}, json_file, indent = 4)
    #    return path

    def get_num_inversions(self, track):
        n = 0
        for p in self.paths:
            for q in self.paths:
                if p != q:
                    if self.tracks[track][p]['lp_value'] > self.tracks[q]['lp_value']:
                        if self.tracks[track][p]['actual_genptype'] > self.tracks[track][q]['actual_genotype']:
                            n += 1
        self.tracks[track]['inversions'] = n
        return n

    def calculate_lp_false_rate(self, track):
        n = 0
        l = 0
        k = 0
        for path in self.paths:
            n += 1
            if self.tracks[track][path]['actual_genotype'] != self.tracks[track][path]['lp_genotype']:
                l += 1
            if self.tracks[track][path]['actual_genotype'] != self.tracks[track][path]['kmeans_genotype']:
                k += 1
        return k / float(n), l / float(n)

    def kmeans(self, track, features):
        print(blue('clustering', track))
        n = len(features[0])
        f = [feature[0] for feature in features]
        K = 3
        if len(sorted(set(f))) < 3:
            K = len(sorted(set(f)))
        while True:
            indices = np.random.randint(0, len(features), size = K)
            print(indices)
            if len(sorted(set(list(indices)))) == K: #duplicate indices
                tmp = [features[j][0] for j in indices]
                if len(sorted(set(tmp))) == K: #duplicate values
                    break
        centroids = np.array([features[j] for j in indices])
        print(centroids)
        old_centroids = np.zeros(centroids.shape)
        candidates = []
        errors = []
        error = 0
        m = 0
        r = 0
        while True:
            m += 1
            if m == 1000 or m == -1 or (old_centroids == centroids).all():
                if m == -1:
                    print(red('empty cluster, skipping round'))
                else:
                    print(green('round', r), cyan('iteration', m))
                    m = 0
                    errors.append(error)
                    candidates.append(centroids)
                    while True:
                        indices = np.random.randint(0, len(features), size = K)
                        if len(sorted(set(list(indices)))) == K:
                            tmp = [features[j][0] for j in indices]
                            if len(sorted(set(tmp))) == K: #duplicate values
                                break
                    centroids = np.array([features[j] for j in indices])
                    old_centroids = np.zeros(centroids.shape)
                r += 1
                if r == 10:
                    break
            error = 0
            clusters = []
            for i in range(len(features)):
                distances = self.distance(features[i], centroids)
                cluster = np.argmin(distances)
                error += distances[cluster]
                clusters.append(cluster)
            old_centroids = copy.deepcopy(centroids)
            for k in range(K):
                points = np.array([features[i] for i in range(len(features)) if clusters[i] == k])
                if len(points) == 0:
                    m = -1
                centroids[k] = np.mean(points, axis = 0)
        #print(centroids)
        choice = np.argmin(errors)
        centroids = candidates[choice]
        clusters = []
        for i in range(len(features)):
            distances = self.distance(features[i], centroids)
            cluster = np.argmin(distances)
            clusters.append(cluster)
        genotypes = []
        norms = np.array([np.linalg.norm(centroid) for centroid in centroids])
        print(norms)
        centroids = []
        for i in range(K):
            centroid = np.argmin(norms)
            norms[centroid] += 1
            centroids.append(centroid)
        print(centroids)
        for i in range(len(features)):
            if clusters[i] == centroids[0]:
                genotypes.append('11')
            elif clusters[i] == centroids[1]:
                genotypes.append('10')
            else:
                genotypes.append('00')
        x = []
        y = []
        a = np.random.randint(0, n, size = 12) / float(n)
        for i in range(12):
            x.append(i % K)
            y.append(a[i])
        for i in range(len(features)):
            x.append(clusters[i])
            y.append(features[i][0])
        try:
            visualizer.violin(x, y, 'kmeans_' + track, self.get_current_job_directory(), 'cluster', 'value')
        except:
            pass
        return genotypes, centroids

    def distance(self, a, b):
        return np.linalg.norm(a - b, axis = 1)

    def reduce(self):
        improvement = []
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'r') as json_file:
                self.tracks[track] = json.load(json_file)
                k, l = self.calculate_lp_false_rate(track)
                improvement.append(l / k)
        visualizer.histogram(improvement, 'false_prediction_lp_vs_kmeans', self.get_current_job_directory(), 'number of events', 'improvement ratio')


    def verify_genotypes(self, path):
        print(path)
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        #output = subprocess.call('verify_genotypes', shell = True, stderr = subprocess.STDOUT, cwd = path)
        for r in ['00', '10', '11']:
            for p in ['00', '10', '11']:
                for track in bed.load_tracks_from_file(os.path.join(path, r + '_as_' + p + '.bed'), [('lp_genotype', None, str), ('lp_value', None, float)]):
                    self.tracks[str(track)][path] = {}
                    self.tracks[str(track)][path]['lp_value'] = track.lp_value
                    self.tracks[str(track)][path]['lp_genotype'] = p
                    self.tracks[str(track)][path]['actual_genotype'] = r

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
