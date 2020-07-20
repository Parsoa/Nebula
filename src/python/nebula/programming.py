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
        c = config.Configuration()
        self.index_kmers()
        self.index_tracks()
        #self.cluster_kmers()
        self.calculate_residual_coverage()
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
                            if kmers[kmer]['type'] == 'inner':
                                self.tracks[track]['confidence'] = 'LOW'
                            else:
                                self.tracks[track]['confidence'] = 'HIGH'
                        if 'svtype' in self.lp_kmers[-1] and self.lp_kmers[-1]['svtype'] != self.tracks[track].svtype:
                            system_print_warning('Kmer with multiple event types..')
                            pass
                        self.lp_kmers[-1]['svtype'] = self.tracks[track].svtype
        print('Selected', len(self.lp_kmers), 'kmers supporting', len(self.tracks), 'tracks for the LP..')
        return self.lp_kmers

    def index_tracks(self):
        n = 0
        tmp = sorted([t for t in self.tracks])
        for track in tmp:
            self.tracks[track]['index'] = n
            self.tracks[track]['kmers'] = []
            self.tracks[track]['cluster'] = -1
            n += 1
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['kmers'].append(index)
        return self.tracks

    #def cluster_kmers(self):
    #    print('Clustering kmers..')
    #    clusters = []
    #    cluster_index = {}
    #    for track in self.tracks:
    #        cluster_index[track] = -1
    #    n = 0
    #    for track in self.tracks:
    #        #print(track, len(self.tracks[track]['kmers']))
    #        cluster = cluster_index[track]
    #        for i in self.tracks[track]['kmers']: # iterate kmers
    #            kmer = self.lp_kmers[i]
    #            #print(kmer)
    #            for track_name in kmer['tracks']: # iterate tracks for kmer
    #                if track_name != track: # kmer has track other than current one
    #                    if cluster_index[track_name] == -1: # other track not yet clustered
    #                        if cluster == -1: # current track not yet clustered, create cluster
    #                            clusters.append([])
    #                            cluster = len(clusters) - 1
    #                            cluster_index[track] = cluster
    #                            cluster_index[track_name] = cluster
    #                            clusters[cluster].append(track)
    #                            clusters[cluster].append(track_name)
    #                            cluster_index[track_name] = cluster
    #                            clusters[cluster].append(track_name)
    #                    else: # append current track to cluster
    #                        if cluster == -1:
    #                            cluster = cluster_index[track_name]
    #                            clusters[cluster].append(track)
    #                        else:
    #                            # merge cluster with this cluster
    #                            c = cluster_index[track_name]
    #                            if c != cluster:
    #                                #print('Reclustering', len(clusters[c]), 'tracks..')
    #                                for _track in clusters[c]:
    #                                    cluster_index[_track] = cluster
    #                                    clusters[cluster].append(_track)
    #                                clusters[c] = None
    #        if cluster == -1:
    #            clusters.append([])
    #            cluster = len(clusters) - 1
    #            cluster_index[track] = cluster
    #            clusters[cluster].append(track)
    #        n += 1
    #        #if n % 100 == 0:
    #        #    print('Processed', n, 'out of', len(self.tracks), 'tracks..')
    #    l = len(list(filter(lambda x: x != None, clusters)))
    #    assert all([cluster_index[t] != -1 for t in self.tracks])
    #    print('Clustered', len(self.lp_kmers), 'kmers and', len(self.tracks), 'tracks into', l, 'clusters..')
    #    self.clusters = clusters

    #def clustered_solve(self):
    #    for cluster in self.clusters:
    #        if len(cluster) == 1:
    #            track = cluster[0] 
    #            n = len(self.tracks[track]['kmers'])
    #            cov = 0
    #            count = 0
    #            kmer_type = 'junction'
    #            #lp_file = open(os.path.join(self.get_current_job_directory(), track + '.lp'), 'w')
    #            for k in self.tracks[track]['kmers']:
    #                kmer = self.lp_kmers[k]
    #                cov += kmer['coverage'] * 0.97
    #                count += kmer['lp_count'] - kmer['coverage'] * kmer['residue']
    #                #lp_file.write(str(kmer['coverage'] * 0.97) + ' '+ track + ' + e = ' + str(kmer['lp_count'] - kmer['coverage'] * kmer['residue']) + '\n')
    #                kmer_type = kmer['type']
    #            if kmer['type'] == 'inner' and self.tracks[track].svtype == 'DEL':
    #                lp = 1.0 - float(count) / float(cov)
    #            if kmer['type'] == 'junction' or self.tracks[track].svtype == 'INS':
    #                lp = float(count) / float(cov)
    #            self.tracks[track]['lp_value'] = lp
    #            self.tracks[track]['lp_genotype'] = self.round_genotype(lp, self.tracks[track]['svtype'])[1]
    #        else:
    #            print('Skipping cluster..')
    #    self.export_solution('algerba')

    def solve(self):
        c = config.Configuration()
        problem, variables = self.generate_mps_linear_program()
        problem.writeLP(os.path.join(self.get_current_job_directory(), 'program_coin.lp'))
        problem.writeMPS(os.path.join(self.get_current_job_directory(), 'program_coin.mps'))
        command = '/share/hormozdiarilab/Codes/NebulousSerendipity/coin/build/bin/clp ' + os.path.join(self.get_current_job_directory(), 'program_coin.mps') + ' -dualsimplex -solution ' + os.path.join(self.get_current_job_directory(), 'solution.mps')
        #output = subprocess.call(command, shell = True)
        #self.import_lp_values()
        #self.round_lp()
        self.verify_genotypes()
        self.correct_genotypes()
        self.export_genotypes('merge')
        #self.export_kmers()

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
        #with open(os.path.join(self.get_current_job_directory(), 'solution_coin.json'), 'w') as json_file:
        #    json.dump({'variables': self.solution}, json_file, indent = 4, sort_keys = True)

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
            bed_file.write('#CHROM\tBEGIN\tEND\tID\tSVTYPE\tLP\tCONFIDENCE\tGENOTYPE\tCOMPATIBILITY\tRENOTYPE\tULTRAHIGH\t0/0\t0/1\t1/1\n')
            for t in bed.sort_tracks(self.tracks):
                #bed_file.write('\t'.join([str(x) for x in [t.chrom, t.begin, t.end, t.id, t.svtype, t['lp_value'], t['confidence'], t['score'], t['lp_genotype']]]) + '\n')
                _t = c.tracks[str(t)]
                bed_file.write('\t'.join([str(x) for x in [t.chrom, t.begin, t.end, t.id, t.svtype, _t['lp'], _t['confidence'], t['genotype'], t['compatibility'], t['renotype'], t['ultrahigh'], t['score_0/0'], t['score_1/0'], t['score_1/1']]]) + '\n')

    # BUG 0000
    # This is problematic
    # Assume an event was genotyped with inner kmers and one of these kmers also belongs to an event that
    # was genotyped with junction kmers alone. Although it won't appear in the second event's kmer set, the
    # first event will add it here. kmer will be mistakenly categorized as junction by Export later as Export
    # assumes all kmers of an event have the same type..
    # Fix: 1. Only add kmer to other event if it matches type
    #      2. Filter kmer entirely
    def export_kmers(self):
        c = config.Configuration()
        tracks = {}
        for index, kmer in enumerate(self.lp_kmers):
            j = all(self.tracks[track]['confidence'] == 'LOW' for track in kmer['tracks'])
            i = all(self.tracks[track]['confidence'] == 'HIGH' for track in kmer['tracks'])
            if not i and not j:
                # ignore kmer
                # or keep it but force the LP to use it
                continue
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

    def calc_kmer_genotype_likelihood(self, kmer, track, genotype, std):
        c = 0 if genotype == '0/0' else 0.5 if genotype == '1/0' else 1.0
        if kmer['type'] == 'junction':
            d = statistics.NormalDistribution(kmer['coverage'] * c, std)
        else:
            if track.svtype == 'INS':
                d = statistics.NormalDistribution(kmer['coverage'] * c, std)
            else:
                d = statistics.NormalDistribution(kmer['coverage'] * (1 - c), std)
        return d.log_pmf(kmer['lp_count'])

    #def verify_genotypes(self):
    #    c = config.Configuration()
    #    genotypes = ['0/0', '1/0', '1/1']
    #    for track in self.tracks:
    #        print('Validating ', track)
    #        likelihoods = {g: 1 for g in genotypes}
    #        for k in self.tracks[track]['kmers']:
    #            kmer = self.lp_kmers[k]
    #            #print(kmer['kmer'], kmer['count'])
    #            if kmer['count'] > kmer['coverage'] + 3 * c.std:
    #                continue
    #            for g in likelihoods:
    #                l = self.calc_kmer_genotype_likelihood(kmer, self.tracks[track], g, c.std)
    #                likelihoods[g] += l
    #                #print(g, l)
    #        g = max(likelihoods, key = likelihoods.get)
    #        #self.tracks[track]['score'] = likelihoods[g] 
    #        #self.tracks[track]['renotype'] = g 
    #        print(likelihoods)
    #        genotype = c.tracks[track].genotype
    #        #genotype = self.tracks[track].lp_genotype
    #        ## log likelihoods are always negative, but increasing as the log is increasing itslef, so higher log likelihood means higher likelihood
    #        ## if test is positive, then g is more likely
    #        ## if test is negative, then g is less likely
    #        ## if this selects zero, means that there was no positive genotype
    #        ratios = [int(likelihoods[g] - likelihoods[genotype]) for g in genotypes if g != genotype]
    #        #print(genotype, ratios)
    #        m = int(max(ratios))
    #        i = ratios.index(m)
    #        self.tracks[track]['score'] = float(max(ratios)) / likelihoods[genotype]
    #        self.tracks[track]['renotype'] = genotypes[i]
    #    #exit()

    def verify_genotypes(self):
        c = config.Configuration()
        genotypes = ['0/0', '1/0', '1/1']
        for track in self.tracks:
            t = c.tracks[track]
            print('Validating ', track)
            kmer_likelihood = {g: 0 for g in genotypes}
            n = 0
            for k in self.tracks[track]['kmers']:
                kmer = self.lp_kmers[k]
                likelihoods = {g: self.calc_kmer_genotype_likelihood(kmer, self.tracks[track], g, c.std) for g in genotypes}
                g = max(likelihoods, key = likelihoods.get)
                kmer_likelihood[g] += 1
                n += 1
            kmer_likelihood = {g: float(kmer_likelihood[g]) / n for g in genotypes}
            g = max(kmer_likelihood, key = kmer_likelihood.get)
            score = kmer_likelihood[t.genotype]
            self.tracks[track]['compatibility'] = score
            for g in kmer_likelihood:
                self.tracks[track]['score_' + g] = kmer_likelihood[g]
            if score < 0.7:
                kmer_likelihood.pop(t.genotype, None)
                g = max(kmer_likelihood, key = kmer_likelihood.get)
                if '1' in t.genotype and not '1' in g:
                    self.tracks[track]['renotype'] = './.'
                elif '1' in g and not '1' in t.genotype:
                    self.tracks[track]['renotype'] = './.'
                else:
                    self.tracks[track]['renotype'] = t.genotype
            else:
                self.tracks[track]['renotype'] = g

    def correct_genotypes(self):
        c = config.Configuration()
        for track in self.tracks:
            t = c.tracks[track]
            self.tracks[track]['genotype'] = t.genotype
            if self.tracks[track]['ultrahigh'] > 10:
                self.tracks[track]['genotype'] = './.'
            if float(self.tracks[track]['ultrahigh']) / len(self.tracks[track]['kmers']) > 0.5:
                self.tracks[track]['genotype'] = './.'
            if t['genotype'] == '0/0' and float(t['lp']) > 0.05 and t['renotype'] != '0/0':
                self.tracks[track]['genotype'] = '1/0'
        self.export_genotypes('corrected')



