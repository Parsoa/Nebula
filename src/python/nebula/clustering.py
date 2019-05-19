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

from nebula import (
    bed,
    cgc,
    config,
    counter,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
    preprocessor
)

from nebula.kmers import *
from nebula.commons import *
from nebula.chromosomes import *
print = pretty_print

from pulp import *
import numpy as np
np.set_printoptions(threshold = np.inf)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class RemnantCounterJob(cgc.CgcCounterJob, map_reduce.BaseGenotypingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'RemnantCounterJob'
    _category = 'preprocessing'
    _previous_job = cgc.CgcCounterJob
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = RemnantCounterJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inner_kmers(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'inner_kmers.json'), 'r') as json_file:
            self.inner_kmers = json.load(json_file)
        self.load_junction_kmers()
        print('Counting', green(len(self.inner_kmers)), 'inner kmers')
        with open(os.path.join(self.get_current_job_directory(), 'pre_inner_kmers.json'), 'w') as json_file:
            _kmers = {}
            for kmers in [self.inner_kmers, self.junction_kmers]:
                for kmer in kmers:
                    _kmers[kmer] = {}
                    _kmers[kmer]['loci'] = {}
                    _kmers[kmer]['tracks'] = kmers[kmer]['tracks']
                    for locus in kmers[kmer]['loci']:
                        _kmers[kmer]['loci'][locus] = {}
                        _kmers[kmer]['loci'][locus]['masks'] = kmers[kmer]['loci'][locus]['masks']
            json.dump(_kmers, json_file, indent = 4)

    def load_junction_kmers(self):
        c = config.Configuration()
        with open(os.path.join(self.get_previous_job_directory(), 'junction_kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                self.junction_kmers[kmer] = kmers[kmer]

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
    _category = 'preprocessing'
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
        print(green(track_name))
        c = config.Configuration()
        t = bed.track_from_name(track_name)
        with open(os.path.join(self.get_previous_job_directory(), track), 'r') as json_file:
            kmers = json.load(json_file)
            lp_kmers = {}
            #inner_kmers = self.select_inner_kmers(kmers['inner_kmers'])
            #TODO: these can be removed
            for kmer in kmers['inner_kmers']:
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers['inner_kmers'][kmer]
                self.lp_kmers[kmer]['reference'] = len(kmers['inner_kmers'][kmer]['loci'])
                self.lp_kmers[kmer]['reduction'] = kmers['inner_kmers'][kmer]['reference']
            for kmer in kmers['junction_kmers']:
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers['junction_kmers'][kmer]
                self.lp_kmers[kmer]['reference'] = len(kmers['junction_kmers'][kmer]['loci'])
                self.lp_kmers[kmer]['reduction'] = kmers['junction_kmers'][kmer]['reference']
            for kmer in kmers['gapped_kmers']:
                lp_kmers[kmer] = True
                self.lp_kmers[kmer] = kmers['gapped_kmers'][kmer]
                self.lp_kmers[kmer]['reference'] = 1 
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump({kmer: self.lp_kmers[kmer] for kmer in lp_kmers}, json_file, indent = 4, sort_keys = True)
        return path

    # Should make sure that non-unique ones are considered the same in all events they appear in
    def select_inner_kmers(self, kmers):
        inner_kmers = list(filter(lambda kmer: kmers[kmer]['type'] == 'inner', kmers))
        unique_kmers = list(filter(lambda kmer: len(kmers[kmer]['loci']) == 1, inner_kmers))
        non_unique_kmers = sorted(list(filter(lambda kmer: len(kmers[kmer]['loci']) != 1, inner_kmers)), key = lambda x: len(kmers[x]['loci']))
        if len(unique_kmers) >= 5:
            inner_kmers = {kmer: kmers[kmer] for kmer in unique_kmers[:min(50, len(unique_kmers))]}
        else:
            inner_kmers = {kmer: kmers[kmer] for kmer in unique_kmers}
            for kmer in non_unique_kmers:
                inner_kmers[kmer] = kmers[kmer]
                if len(inner_kmers) > 50:
                    break
        return inner_kmers

    def calculate_residual_coverage(self):
        c = config.Configuration()
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['coverage'] = c.coverage
            kmer['residue'] = kmer['reference'] - r
            kmer['weight'] = 1.0
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
            #l = len(self.lp_kmers)
            #l = len(list(filter(lambda i: self.lp_kmers[i]['type'] != 'gapped', self.tracks[kmer['tracks'].keys()[0]]['kmers'])))
            #l = l if l else 1
            #kmer['weight'] = l if kmer['type'] == 'gapped' else 1.0
        self.calculate_error_weights()

    def calculate_error_weights(self):
        # assuming no kmers are shared
        for track in self.tracks:
            self.add_weights_to_track(track)

    def add_weights_to_track(self, track):
        l = self.tracks[track]['kmers']
        inner = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'inner']
        gapped = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'gapped']
        junction = [self.lp_kmers[i] for i in l if self.lp_kmers[i]['type'] == 'junction']
        for kmers in [inner, gapped, junction]:
            for kmer in kmers:
                kmer['weight'] = float(len(inner) + len(gapped) + len(junction)) / (3.0 * len(kmers))

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
        for index, kmer in enumerate(self.lp_kmers):
            self.add_error_absolute_value_constraints(problem, index)
            if kmer['type'] == 'inner':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks'])) # Coverage
                ind.append(len(self.tracks) + index) # Objective
                val = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks'])) #Coverage corrected for errors
                val.append(1.0) #Objective
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - kmer['coverage'] * kmer['residue']],
                    senses = ['E']
                )
            if kmer['type'] == 'gapped' and kmer['side'] == 'outer' or kmer['type'] == 'junction':
                ind = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                ind.append(len(self.tracks) + index)
                val = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                val.append(1.0)
                problem.linear_constraints.add(
                    lin_expr = [cplex.SparsePair(
                        ind = ind,
                        val = val,
                    )],
                    rhs = [kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))],
                    senses = ['E']
                )
        return problem

    def generate_mps_linear_program(self):
        c = config.Configuration()
        problem = LpProblem("Nebula", LpMinimize)
        i = 0
        names = [''] * len(self.tracks)
        variables = [None] * (len(self.tracks) + 2 * len(self.lp_kmers))
        for track in self.tracks:
            tokens = track.split('_')
            variables[self.tracks[track]['index']] = LpVariable('c' + tokens[1], 0, 1)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintLE, 'c' + tokens[1] + '_ub', 1.0)
            problem += LpConstraint(LpAffineExpression([(variables[self.tracks[track]['index']], 1.0)]), LpConstraintGE, 'c' + tokens[1] + '_lb', 0.0)
            i += 1
        # error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('e' + str(index))
            i += 1
        # absolute value of the error variables
        for index, kmer in enumerate(self.lp_kmers):
            variables[i] = LpVariable('l' + str(index))
            i += 1
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], kmer['weight']) for index, kmer in enumerate(self.lp_kmers)])
        problem += expr
        for i, kmer in enumerate(self.lp_kmers):
            self.add_mps_error_absolute_value_constraints(problem, variables, i)
            if kmer['type'] == 'inner':
                indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                indices.append(len(self.tracks) + i)
                coeffs = list(map(lambda track: kmer['coverage'] * kmer['tracks'][track] * (1.0 - 0.03), kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['count'] - kmer['coverage'] * kmer['residue']
                expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
            else:
                indices = list(map(lambda track: self.tracks[track]['index'], kmer['tracks']))
                indices.append(len(self.tracks) + i)
                coeffs = list(map(lambda track: -1 * kmer['coverage'] * kmer['tracks'][track], kmer['tracks']))
                coeffs.append(1.0)
                rhs = kmer['count'] - sum(list(map(lambda track: kmer['coverage'] * kmer['tracks'][track], kmer['tracks'])))
                expr = LpAffineExpression([(variables[v], coeffs[index]) for index, v in enumerate(indices)])
                problem += LpConstraint(expr, LpConstraintEQ, 'k' + str(i), rhs)
        return problem, variables

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class ClusterIntegerProgrammingAnalysisJob(programming.IntegerProgrammingJob):

    _name = 'ClusterIntegerProgrammingAnalysisJob'
    _category = 'preprocessing'
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
        paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/simulation/HG00513.hg19.chr17.DEL.bed/Seed' + str(i) + 'SnpError0.001/20x/ClusterIntegerProgrammingJob' for i in range(2001, 2251)]
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
            visualizer.violin(x, y, 'cluster_' + track, self.get_current_job_directory(), 'genotype', 'lp value')
            visualizer.histogram(y, 'cluster_' + track, self.get_current_job_directory(), 'genotype', 'lp value', step = 0.05)
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
    _category = 'preprocessing'
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
        self.paths = ['/share/hormozdiarilab/Codes/NebulousSerendipity/output/simulation/HG00513.hg19.chr17.DEL.bed/Seed' + str(i) + 'SnpError0.001/20x/ClusterIntegerProgrammingJob' for i in range(2001, 2251)]
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
        self.get_num_inversions(track_name)
        self.get_num_mistakes(track_name)
        path = os.path.join(self.get_current_job_directory(), track_name + '.json')
        with open(path, 'w') as json_file:
            json.dump(track, json_file, indent = 4, sort_keys = True)
        return path

    def get_num_inversions(self, track):
        n = 0
        for p in self.paths:
            for q in self.paths:
                if p != q:
                    if self.tracks[track][p]['lp_value'] > self.tracks[track][q]['lp_value']:
                        if self.tracks[track][p]['actual_genotype'] > self.tracks[track][q]['actual_genotype']:
                            n += 1
        self.tracks[track]['inversions'] = n
        return n

    def get_num_mistakes(self, track):
        n = 0
        for p in self.paths:
            if self.tracks[track][p]['lp_genotype'] != self.tracks[track][p]['actual_genotype']:
                n += 1
        self.tracks[track]['mistakes'] = n
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
        for track in self.tracks:
            self.tracks[track] = json.load(open(os.path.join(self.get_current_job_directory(), track + '.json'), 'r'))
        performance = []
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'r') as json_file:
                self.tracks[track] = json.load(json_file)
                i, m = self.tracks[track]['inversions'], self.tracks[track]['mistakes']
                if m != 0:
                    performance.append(float(i) / float(m))
        visualizer.histogram(performance, 'false_prediction_lp_vs_kmeans', self.get_current_job_directory(), 'number of events', 'improvement ratio')

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
