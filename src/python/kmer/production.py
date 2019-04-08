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

from kmer import (
    bed,
    config,
    gapped,
    counter,
    dynamic,
    reduction,
    simulator,
    counttable,
    map_reduce,
    statistics,
    visualizer,
    programming,
)

from kmer.kmers import *
from kmer.commons import *
from kmer.chromosomes import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MixKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'MixKmersJob'
    _category = 'programming'
    _previous_job = None 
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = MixKmersJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.half_mers = {}
        self.depth_kmers = {}
        self.inner_kmers = {}
        self.gapped_kmers = {}
        self.junction_kmers = {}
        self.load_kmers()
        self.merge_kmers()
        #self.export_tracks()
        exit()

    def load_kmers(self):
        c = config.Configuration()
        job = reduction.FilterLociIndicatorKmersJob()
        with open(os.path.join(job.get_current_job_directory(), 'kmers.json'), 'r') as json_file:
            self.inner_kmers = json.load(json_file)
        self.load_junction_kmers()
        self.load_depth_of_coverage_kmers()
        #self.load_gapped_kmers()

    def merge_kmers(self):
        kmers = {}
        kmers['half_mers'] = self.half_mers
        kmers['depth_kmers'] = self.depth_kmers
        kmers['inner_kmers'] = self.inner_kmers
        kmers['gapped_kmers'] = self.gapped_kmers
        kmers['junction_kmers'] = self.junction_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump(kmers, json_file, indent = 4, sort_keys = True)

    def load_junction_kmers(self):
        c = config.Configuration()
        job = dynamic.CountJunctionKmersJob()
        with open(os.path.join(job.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                kmer = canonicalize(kmer)
                if kmer in self.inner_kmers:
                    self.inner_kmers.pop(kmer, None)
                else:
                    self.junction_kmers[kmer] = kmers[kmer]
        #with open(os.path.join(self.get_current_job_directory(), 'junction_kmers.json'), 'w') as json_file:
        #    json.dump(self.junction_kmers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
        #    json.dump(self.inner_kmers, json_file, indent = 4)

    def load_depth_of_coverage_kmers(self):
        n = 100000
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count == 1 and kmer.find('N') == -1:
                canon = canonicalize(kmer)
                if not canon in self.inner_kmers and not canon in self.junction_kmers:
                    self.depth_kmers[canon] = {'loci': {}, 'count': 0}
                    n -= 1
                    if n == 0:
                        break
        print('Counting', green(len(self.depth_kmers)), 'depth kmers')
        #with open(os.path.join(self.get_current_job_directory(), 'depth_kmers.json'), 'w') as json_file:
        #    json.dump(self.depth_kmers, json_file, sort_keys = True, indent = 4)
        self.unload_reference_counts_provider()

    def load_gapped_kmers(self):
        c = config.Configuration()
        job = gapped.CountUniqueGappedKmersJob()
        tracks = job.load_previous_job_results()
        for track in tracks:
            with open(os.path.join(job.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for kmer in kmers:
                    if kmers[kmer]['gap'] != -1:
                        left = kmer[:c.hsize]
                        right = kmer[-c.hsize:]
                        self.gapped_kmers[kmer] = kmers[kmer]
                        self.gapped_kmers[kmer]['count'] = 0
                        self.gapped_kmers[kmer]['doubt'] = 0
                        if not left in self.half_mers:
                            self.half_mers[left] = {}
                        self.half_mers[left][right] = kmer 
                        left = reverse_complement(left)
                        right = reverse_complement(right)
                        if not right in self.half_mers:
                            self.half_mers[right] = {}
                        self.half_mers[right][left] = kmer
        print('Counting', green(len(self.gapped_kmers)), 'gapped kmers')
        #with open(os.path.join(self.get_current_job_directory(), 'half_mers.json'), 'w') as json_file:
        #    json.dump(self.half_mers, json_file, indent = 4)
        #with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
        #    json.dump(self.gapped_kmers, json_file, indent = 4)

    def export_tracks(self):
        c = config.Configuration()
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.gapped_kmers:
            for track in self.gapped_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['gapped_kmers'][kmer] = self.gapped_kmers[kmer]
                self.tracks[track]['gapped_kmers'][kmer]['type'] = 'gapped'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'tracks.json'), 'w') as json_file:
            json.dump(self.tracks, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MixCounterJob(map_reduce.FirstGenotypingJob, counter.BaseExactCountingJob):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'MixCounterJob'
    _category = 'programming'
    _previous_job = None 
    _counter_mode = 4

    @staticmethod
    def launch(**kwargs):
        job = MixCounterJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.half_mers = {}
        self.depth_kmers = {}
        self.inner_kmers = {}
        self.gapped_kmers = {}
        self.junction_kmers = {}
        self.load_inner_kmers()
        self.load_gapped_kmers()
        self.round_robin()

    def load_inner_kmers(self):
        c = config.Configuration()
        job = reduction.CountLociIndicatorKmersJob()
        with open(os.path.join(job.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            self.inner_kmers = json.load(json_file)
        self.load_junction_kmers()
        #self.load_depth_of_coverage_kmers()
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
            #for kmer in self.depth_kmers:
            #    _kmers[kmer] = {}
            #    _kmers[kmer]['loci'] = {}
            json.dump(_kmers, json_file, indent = 4)

    def load_depth_of_coverage_kmers(self):
        n = 100000
        self.load_reference_counts_provider() 
        for kmer, count in self.reference_counts_provider.stream_kmers():
            if count == 1 and kmer.find('N') == -1:
                canon = canonicalize(kmer)
                if not kmer in self.inner_kmers and not reverse_complement(kmer) in self.inner_kmers:
                    self.depth_kmers[canon] = {'loci': {}}
                    n -= 1
                    if n == 0:
                        break
        print('Counting', green(len(self.depth_kmers)), 'depth kmers')
        self.unload_reference_counts_provider()

    def load_junction_kmers(self):
        c = config.Configuration()
        job = dynamic.CountJunctionKmersJob()
        with open(os.path.join(job.get_previous_job_directory(), 'kmers.json'), 'r') as json_file:
            kmers = json.load(json_file)
            for kmer in kmers:
                kmer = canonicalize(kmer)
                if kmer in self.inner_kmers:
                    self.inner_kmers.pop(kmer, None)
                else:
                    self.junction_kmers[kmer] = kmers[kmer]

    # TODO: Remove side from this
    # There shouldn't be any overlap between gapped kmers and other kmers because of the filtering applied to them
    def load_gapped_kmers(self):
        c = config.Configuration()
        job = gapped.CountUniqueGappedKmersJob()
        tracks = job.load_previous_job_results()
        for track in tracks:
            with open(os.path.join(job.get_previous_job_directory(), tracks[track]), 'r') as json_file:
                kmers = json.load(json_file)
                for side in ['inner', 'outer']:
                    for kmer in kmers[side]:
                        if kmers[side][kmer]['gap'] != -1:
                            left = kmer[:c.hsize]
                            right = kmer[-c.hsize:]
                            self.gapped_kmers[kmer] = kmers[side][kmer]
                            self.gapped_kmers[kmer]['count'] = 0
                            self.gapped_kmers[kmer]['doubt'] = 0
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

    def merge_count(self, kmer, tokens):
        if kmer in self.gapped_kmers:
            count = tokens[0] 
            self.gapped_kmers[kmer]['count'] += count
        else:
            count = tokens[0] 
            total = tokens[1] 
            canon = canonicalize(kmer)
            if canon in self.inner_kmers:
                self.inner_kmers[canon]['count'] += count / 2
                self.inner_kmers[canon]['total'] += total / 2
            else:
                self.junction_kmers[canon]['count'] += count / 2
                self.junction_kmers[canon]['total'] += total / 2

    def reduce(self):
        c = config.Configuration()
        self.merge_counts()
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
            json.dump(self.inner_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
            json.dump(self.gapped_kmers, json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'junction_kmers.json'), 'w') as json_file:
            json.dump(self.junction_kmers, json_file, indent = 4)
        self.tracks = {}
        for kmer in self.inner_kmers:
            for track in self.inner_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['inner_kmers'][kmer] = self.inner_kmers[kmer]
                self.tracks[track]['inner_kmers'][kmer]['type'] = 'inner'
        for kmer in self.gapped_kmers:
            for track in self.gapped_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['gapped_kmers'][kmer] = self.gapped_kmers[kmer]
                self.tracks[track]['gapped_kmers'][kmer]['type'] = 'gapped'
        for kmer in self.junction_kmers:
            for track in self.junction_kmers[kmer]['tracks']:
                if not track in self.tracks:
                    self.tracks[track] = {'inner_kmers': {}, 'gapped_kmers': {}, 'junction_kmers': {}}
                self.tracks[track]['junction_kmers'][kmer] = self.junction_kmers[kmer]
                self.tracks[track]['junction_kmers'][kmer]['type'] = 'junction'
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump({track: track + '.json' for track in self.tracks}, json_file, indent = 4)

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class MixIntegerProgrammingJob(programming.IntegerProgrammingJob):

    _name = 'MixIntegerProgrammingJob'
    _category = 'programming'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        c = config.Configuration()
        job = MixIntegerProgrammingJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def execute(self):
        c = config.Configuration()
        self.prepare()
        self.create_output_directories()
        self.gapped_kmers = []
        self.gapped_tracks = {}
        self.inner_kmers = []
        self.inner_tracks = {}
        #
        self.load_gapped_kmers()
        self.load_inner_kmers()
        self.merge_kmers()
        self.solve()
        exit()

    def load_gapped_kmers(self):
        print(cyan('=============================================================================================='))
        print(cyan('Loading gapped kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.gapped_kmers = payload['gapped_kmers']
                self.gapped_tracks = payload['tracks']
            return
        self.gapped_kmers_solver = gapped.GappedKmersIntegerProgrammingJob()
        self.gapped_kmers_solver.create_output_directories()
        self.gapped_kmers_solver.load_inputs()
        self.gapped_kmers_solver.distribute_workload()
        self.gapped_kmers_solver.wait_for_children()
        self.gapped_kmers = copy.deepcopy(self.gapped_kmers_solver.index_kmers())
        self.gapped_tracks = copy.deepcopy(self.gapped_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'gapped_kmers.json'), 'w') as json_file:
            json.dump({'gapped_kmers': self.gapped_kmers, 'tracks': self.gapped_tracks}, json_file, indent = 4)

    def load_inner_kmers(self):
        print(cyan('=============================================================================================='))
        print(cyan('Loading inner kmers...'))
        print(cyan('=============================================================================================='))
        if self.resume_from_reduce:
            with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'r') as json_file:
                payload = json.load(json_file)
                self.inner_kmers = payload['inner_kmers']
                self.inner_tracks = payload['tracks']
            return
        self.inner_kmers_solver = reduction.LociIndicatorKmersIntegerProgrammingJob()
        self.inner_kmers_solver.create_output_directories()
        self.inner_kmers_solver.load_inputs()
        self.inner_kmers_solver.distribute_workload()
        self.inner_kmers_solver.wait_for_children()
        self.inner_kmers = copy.deepcopy(self.inner_kmers_solver.index_kmers())
        self.inner_tracks = copy.deepcopy(self.inner_kmers_solver.index_tracks())
        with open(os.path.join(self.get_current_job_directory(), 'inner_kmers.json'), 'w') as json_file:
            json.dump({'inner_kmers': self.inner_kmers, 'tracks': self.inner_tracks}, json_file, indent = 4)

    def merge_kmers(self):
        c = config.Configuration()
        self.tracks = {}
        self.tracks.update(self.inner_tracks)
        self.tracks.update(self.gapped_tracks)
        n = 0
        tmp = sorted([t for t in self.tracks])
        self.lp_kmers = self.gapped_kmers + self.inner_kmers
        for track in tmp:
            self.tracks[track] = {'index': n,
                'inner_kmers': [index + len(self.gapped_kmers) for index in self.inner_tracks[track]['kmers']] if track in self.inner_tracks else [],
                'gapped_kmers': self.gapped_tracks[track]['kmers'] if track in self.gapped_tracks else [],
            }
            n += 1
            with open(os.path.join(self.get_current_job_directory(), 'kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)
        self.lp_kmers = self.inner_kmers + self.gapped_kmers
        with open(os.path.join(self.get_current_job_directory(), 'kmers.json'), 'w') as json_file:
            json.dump({'kmers': self.lp_kmers}, json_file, indent = 4)
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('Solving composite model on'), blue(len(self.lp_kmers)), cyan('kmers and '), blue(len(self.tracks)), cyan('tracks'))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        print(cyan('=============================================================================================='))
        for kmer in self.lp_kmers:
            r = 0
            for track in kmer['tracks']:
                r += kmer['tracks'][track]
            kmer['coverage'] = c.coverage if kmer['type'] == 'gapped' else kmer['coverage']
            kmer['residue'] = 0 if kmer['type'] == 'gapped' else kmer['reference'] - r
            kmer['count'] = min(kmer['count'], kmer['coverage'] * kmer['reference'])
            l = len(self.tracks[kmer['tracks'].keys()[0]]['inner_kmers'])
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

    def add_snp_linear_constraints(self, problem):
        problem.variables.add(names = ['i' + str(index) for index, kmer in enumerate(self.lp_kmers)], ub = [1.0] * len(self.lp_kmers))
        offset = len(self.tracks) + 2 * len(self.lp_kmers)
        for track in self.tracks:
            ind = [len(self.tracks) + 2 * len(self.lp_kmers) + index for index in self.tracks[track]['inner_kmers'] + self.tracks[track]['gapped_kmers']]
            problem.linear_constraints.add(
                lin_expr = [cplex.SparsePair(
                    ind = ind,
                    val = [1.0] * len(ind),
                )],
                rhs = [math.floor(0.1 * len(ind))],
                senses = ['L']
            )

    def solve(self):
        programming.IntegerProgrammingJob.solve(self)
        for track in self.tracks:
            self.tracks[track]['kmers'] = {'inner_kmers': {}, 'gapped_kmers': {}}
            self.tracks[track]['errors'] = {'inner_kmers': [], 'gapped_kmers': []}
            self.tracks[track]['indices'] = {'inner_kmers': [], 'gapped_kmers': []}
        for index, kmer in enumerate(self.lp_kmers):
            for track in kmer['tracks']:
                self.tracks[track]['kmers'][kmer['type'] + '_kmers'][kmer['kmer']] = kmer
                self.tracks[track]['errors'][kmer['type'] + '_kmers'].append(abs(self.errors[index]))
                self.tracks[track]['indices'][kmer['type'] + '_kmers'].append(index)
        for track in self.tracks:
            with open(os.path.join(self.get_current_job_directory(), 'kmers_' + track + '.json'), 'w') as json_file:
                json.dump(self.tracks[track], json_file, indent = 4)

    def calculate_confidence_score(self, track):
        score = 0
        for error in track['errors']['inner_kmers']:
            score += 1.0 / (1 + error)
        for error in track['errors']['gapped_kmers']:
            score += 10 / (1 + error)
        return score

# ============================================================================================================================ #
# Main
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    getattr(sys.modules[__name__], c.job).launch(resume_from_reduce = c.reduce)
