from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import argparse
import operator
import traceback
import subprocess

from nebula import (
    bed,
    config,
    counter,
    map_reduce,
    statistics,
    visualizer,
)

from nebula.kmers import *
from nebula.logger import *
from nebula.chromosomes import *
print = pretty_print

from pulp import *

# ============================================================================================================================ #
# ============================================================================================================================ #
# Models the problem as an integer program and uses CPLEX to solve it
# This won't need any parallelization
# ============================================================================================================================ #
# ============================================================================================================================ #

class IntegerProgrammingJob(map_reduce.Job):

    _name = 'IntegerProgrammingJob'
    _category = 'preprocessing'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = IntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.tracks = self.load_previous_job_results()
        self.round_robin(self.tracks)
        self.resume_from_reduce = False
        self.lp_kmers = {}

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(self.lp_kmers, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def reduce(self):
        self.index_kmers()
        self.index_tracks()
        #self.cluster_kmers()
        #self.clustered_solve()
        self.solve()

    def index_kmers(self):
        print('Indexing kmers..')
        c = config.Configuration()
        self.tracks = {}
        self.lp_kmers = []
        index = {}
        for batch in self.load_output():
            kmers = batch
            for kmer in kmers:
                if not kmer in index:
                    index[kmer] = len(self.lp_kmers)
                    self.lp_kmers.append(copy.deepcopy(kmers[kmer]))
                    self.lp_kmers[-1]['kmer'] = kmer
                    for track in kmers[kmer]['tracks']:
                        if not track in self.tracks:
                            self.tracks[track] = bed.track_from_id(track)
                            self.tracks[track]['confidence'] = 'N/A'
                            #if kmers[kmer]['type'] == 'inner':
                            #    if self.tracks[track]['confidence'] == 'HIGH':
                            #        self.tracks[track]['confidence'] = 'MIX'
                            #    else:
                            #        self.tracks[track]['confidence'] = 'LOW'
                            #else:
                            #    if self.tracks[track]['confidence'] == 'LOW':
                            #        self.tracks[track]['confidence'] = 'MIX'
                            #    else:
                            #        self.tracks[track]['confidence'] = 'HIGH'
                        if 'svtype' in self.lp_kmers[-1] and self.lp_kmers[-1]['svtype'] != self.tracks[track].svtype:
                            system_print_warning('Kmer with multiple event types..')
                            pass
                        self.lp_kmers[-1]['svtype'] = self.tracks[track].svtype
        print('Selected', len(self.lp_kmers), 'kmers supporting', len(self.tracks), 'tracks for the LP..')
        return self.lp_kmers

    def index_tracks(self):
        c = config.Configuration()
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track]['index'] = n
            self.tracks[track]['kmers'] = []
            self.tracks[track]['cluster'] = -1
            n += 1
        m = 0
        ultrahigh_tracks = {}
        for index, kmer in enumerate(self.lp_kmers):
            self.calculate_residual_coverage(kmer)
            for track in kmer['tracks']:
                self.tracks[track]['kmers'].append(index)
                if kmer['count'] > kmer['coverage'] + 3 * c.std:
                    m += 1
                    if not track in ultrahigh_tracks:
                        ultrahigh_tracks[track] = 0
                    ultrahigh_tracks[track] += 1
        print(len(ultrahigh_tracks), 'events with ultra high-frequency kmers.')
        for track in self.tracks:
            self.tracks[track]['ultrahigh'] = 0
        for track in ultrahigh_tracks:
            self.tracks[track]['ultrahigh'] = ultrahigh_tracks[track]
        return self.tracks

    def cluster_kmers(self):
        print('Clustering kmers..')
        clusters = []
        cluster_index = {}
        for track in self.tracks:
            cluster_index[track] = -1
        n = 0
        for track in self.tracks:
            cluster = cluster_index[track]
            for i in self.tracks[track]['kmers']: # iterate kmers
                kmer = self.lp_kmers[i]
                for track_name in kmer['tracks']: # iterate tracks for kmer
                    if track_name != track: # kmer has track other than current one
                        if cluster_index[track_name] == -1: # other track not yet clustered
                            if cluster == -1: # current track not yet clustered, create cluster
                                clusters.append([])
                                cluster = len(clusters) - 1
                                cluster_index[track] = cluster
                                cluster_index[track_name] = cluster
                                clusters[cluster].append(track)
                                clusters[cluster].append(track_name)
                                cluster_index[track_name] = cluster
                                clusters[cluster].append(track_name)
                        else: # append current track to cluster
                            if cluster == -1:
                                cluster = cluster_index[track_name]
                                clusters[cluster].append(track)
                            else: # merge cluster with this cluster
                                c = cluster_index[track_name]
                                if c != cluster:
                                    for _track in clusters[c]:
                                        cluster_index[_track] = cluster
                                        clusters[cluster].append(_track)
                                    clusters[c] = None
            if cluster == -1:
                clusters.append([])
                cluster = len(clusters) - 1
                cluster_index[track] = cluster
                clusters[cluster].append(track)
            n += 1
            #if n % 100 == 0:
            #    print('Processed', n, 'out of', len(self.tracks), 'tracks..')
        l = len(list(filter(lambda x: x != None, clusters)))
        assert all([cluster_index[t] != -1 for t in self.tracks])
        print('Clustered', len(self.lp_kmers), 'kmers and', len(self.tracks), 'tracks into', l, 'clusters..')
        self.clusters = clusters

    def clustered_solve(self):
        _tracks = self.tracks
        _lp_kmers = self.lp_kmers
        self.tracks = {}
        self.lp_kmers = []
        for i, cluster in enumerate(self.clusters):
            print('Solving for cluster', i, 'with', len(cluster), 'tracks..')
            self.tracks = {track: _tracks[track] for track in cluster}
            self.lp_kmers = [kmer for kmer in _lp_kmers if any([track in kmer['tracks'] for track in cluster])]
            self.index_tracks()
            problem, variables = self.generate_mps_linear_program()
            problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin.lp'))
            problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin.mps'))
            command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
            output = subprocess.call(command, shell = True)
            self.import_lp_values()
            self.round_lp()
            for track in self.tracks:
                _tracks[track]['lp_value'] = self.tracks[track]['lp_value']
                _tracks[track]['lp_genotype'] = self.tracks[track]['lp_genotype']
        self.tracks = _tracks
        self.lp_kmers = _lp_kmers
        self.verify_genotypes()
        self.correct_genotypes()
        self.export_genotypes('merge')
        self.export_kmers()

    def solve(self):
        c = config.Configuration()
        problem, variables = self.generate_mps_linear_program()
        problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin.lp'))
        problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin.mps'))
        command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
        output = subprocess.call(command, shell = True)
        self.import_lp_values()
        self.round_lp()
        self.verify_genotypes()
        self.export_genotypes('merge')
        self.correct_genotypes()
        self.export_kmers()

    def add_mps_error_absolute_value_constraints(self, problem, variables, index):
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], 1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_1_' + str(index), 0) 
        expr = LpAffineExpression([(variables[len(self.tracks) + len(self.lp_kmers) + index], 1.0), (variables[len(self.tracks) + index], -1.0)])
        problem += LpConstraint(expr, LpConstraintGE, 'abs_2_' + str(index), 0) 

    # COIN doesn't supply values for certain variables
    def import_lp_values(self, path = 'solution.mps'):
        c = config.Configuration()
        var_index = {}
        self.solution = [0.0] * (len(self.tracks) + 2 * len(self.lp_kmers))
        regex = re.compile('[^a-zA-Z0-9]')
        for track in self.tracks:
            name = regex.sub('_', track)
            var_index[name] = self.tracks[track]['index']
        with open(os.path.join(self.get_current_job_directory(), path), 'r') as f:
            status = f.readline()
            objective = f.readline()
            line = f.readline()
            while(line):
                tokens = line.split()
                name = tokens[1]
                index = int(tokens[0])
                value = float(tokens[2])
                # tracks are indexed lexicographically, different from iteration index
                if name[0] == 'c':
                    self.solution[var_index[name[1:]]] = value
                else:
                    self.solution[index] = value
                line = f.readline()

    def round_lp(self):
        c = config.Configuration()
        system_print_normal('Assigning genotypes for', len(self.tracks), 'tracks..')
        for track in self.tracks:
            index = self.tracks[track]['index']
            self.tracks[track]['lp_value'] = self.solution[index]
            self.tracks[track]['lp_genotype'] = self.round_genotype(self.solution[index], self.tracks[track]['svtype'])[1]

    def export_genotypes(self, name):
        c = config.Configuration()
        name = name + '.bed' if not c.cgc else 'genotypes_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.bed'
        path = os.path.join(c.workdir if c.cgc else self.get_current_job_directory(), name)
        with open(path, 'w') as bed_file:
            bed_file.write('\t'.join(['#CHROM', 'BEGIN', 'END', 'ID', 'SVTYPE', 'LP', 'CONFIDENCE', 'GENOTYPE', 'COMPATIBILITY', 'RENOTYPE', 'ULTRAHIGH', '0/0', '1/0', '1/1']) + '\n')
            for t in bed.sort_tracks(self.tracks):
                #_t = c.tracks[str(t)]
                bed_file.write('\t'.join([str(x) for x in [t.chrom, t.begin, t.end, t.id, t.svtype, t['lp_value'], t['confidence'], t['lp_genotype'], t['compatibility'], t['renotype'], t['ultrahigh'], t['score_0/0'], t['score_1/0'], t['score_1/1']]]) + '\n')

    def export_kmers(self):
        c = config.Configuration()
        tracks = {}
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                if not track in tracks:
                    tracks[track] = {}
                tracks[track][kmer['kmer']] = kmer
        for track in tracks:
            path = os.path.join(self.get_current_job_directory(), track + '.json')
            with open(path, 'w') as json_file:
                json.dump(tracks[track], json_file, indent = 4)
        if c.cgc:
            name = 'lp_kmer_' + (c.fastq.split('/')[-1] if c.fastq else c.bam.split('/')[-1]) + '.json'
            with open(os.path.join(c.workdir, name), 'w') as json_file:
                json.dump({kmer['kmer']: kmer for kmer in self.lp_kmers}, json_file, indent = 4, sort_keys = True)

    def calc_kmer_genotype_likelihood(self, kmer, genotype, std):
        c = 0 if genotype == '0/0' else 0.5 if genotype == '1/0' else 1.0
        s = 0.25 if genotype == '0/0' else 0.50 if genotype == '1/0' else 1.00
        if kmer['type'] == 'junction':
            d = statistics.NormalDistribution(kmer['gc_coverage'] * c, std * s)
        else:
            if kmer['trend'] == 'upward':
                d = statistics.NormalDistribution(kmer['gc_coverage'] * c, std * s)
            else:
                d = statistics.NormalDistribution(kmer['gc_coverage'] * (1 - c), std * (0.25 / s))
        return d.log_pmf(kmer['count'])

    def verify_genotypes(self):
        c = config.Configuration()
        genotypes = ['0/0', '1/0', '1/1']
        for track in self.tracks:
            #print('Validating ', track)
            kmer_likelihood = {g: 0 for g in genotypes}
            n = 0
            for k in self.tracks[track]['kmers']:
                kmer = self.lp_kmers[k]
                likelihoods = {g: self.calc_kmer_genotype_likelihood(kmer, g, c.std) for g in genotypes}
                g = max(likelihoods, key = likelihoods.get)
                kmer_likelihood[g] += 1
                n += 1
            kmer_likelihood = {g: float(kmer_likelihood[g]) / n for g in genotypes}
            g = max(kmer_likelihood, key = kmer_likelihood.get)
            lp_genotype = self.tracks[track]['lp_genotype']
            score = kmer_likelihood[lp_genotype]
            self.tracks[track]['compatibility'] = score
            for g in kmer_likelihood:
                self.tracks[track]['score_' + g] = kmer_likelihood[g]
            if score < 0.7:
                kmer_likelihood.pop(lp_genotype, None)
                g = max(kmer_likelihood, key = kmer_likelihood.get)
                if ('1' in lp_genotype and not '1' in g) or ('1' in g and not '1' in lp_genotype):
                    self.tracks[track]['renotype'] = './.'
                else:
                    self.tracks[track]['renotype'] = lp_genotype
            else:
                self.tracks[track]['renotype'] = lp_genotype

    def correct_genotypes(self):
        c = config.Configuration()
        for track in self.tracks:
            if self.tracks[track]['ultrahigh'] > 10:
                self.tracks[track]['lp_genotype'] = './.'
            if float(self.tracks[track]['ultrahigh']) / len(self.tracks[track]['kmers']) > 0.5:
                self.tracks[track]['lp_genotype'] = './.'
            if self.tracks[track]['lp_genotype'] == '0/0' and float(self.tracks[track]['lp_value']) > 0.05 and self.tracks[track]['renotype'] != '0/0':
                self.tracks[track]['lp_genotype'] = '0/1'
        self.export_genotypes('corrected')

