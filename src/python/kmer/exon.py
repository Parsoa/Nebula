# ============================================================================================================================ #
# ============================================================================================================================ #
# MapReduce job that extracts the kmers for the exons specified in a BED file
# ============================================================================================================================ #
# ============================================================================================================================ #

class ExtractExonKmersJob(map_reduce.Job):

    # ============================================================================================================================ #
    # job-specific stuff
    # ============================================================================================================================ #

    def parse_track_file(self, path):
        tracks = {}
        with open(path) as bed_file:
            line = bed_file.readline()
            while line:
                tokens = line.split()
                track = bed.BedTrack(tokens[-3], int(tokens[-2]), int(tokens[-1]))
                if not track.name in tracks:
                    tracks[track.name] = track
                line = bed_file.readline()
        return tracks

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = KmerNormalDistributionFittingJob(job_name = 'ExtractExonKmersJob_', previous_job_name = "", **kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def prepare(self):
        self.kmers = {}

    def check_cli_arguments(self, args):
        # --reference: the reference genome to use
        # --bed: the BED file to use
        pass

    def find_thread_count(self):
        c = config.Configuration()
        self.num_threads = c.max_threads

    def load_inputs(self):
        c = config.Configuration()
        # 
        tracks = self.parse_track_file(c.bed_file)
        for i in range(0, self.num_threads):
            self.batch[i] = {} # avoid overrding extra methods from MapReduce
        #
        index = 0
        for track in tracks:
            self.batch[index][track] = tracks[track]
            index += 1
            if index == self.num_threads:
                index = 0

    def transform(self, track, track_name):
        c = config.Configuration()
        seq = bed.extract_sequence(track)
        for kmer in get_all_kmers(seq, c.ksize):
            if not kmer in self.kmers:
                self.kmers[kmer] = 0
            self.kmers[kmer] += 1
        # we are not chaning the input track, storing results in a object property instead
        return track

    def output_batch(self, batch):
        with open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w') as json_file:
            json.dump(self.kmers, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
        exit()

    def reduce(self):
        c = config.Configuration()
        kmers = {}
        # merge all the kmer counts from previous steps
        for i in range(0, self.num_threads):
            print('batch', i)
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    for kmer in batch:
                        if not kmer in kmers:
                            kmers[kmer] = 0
                        kmers[kmer] += batch[kmer]
        with open(os.path.join(self.get_current_job_directory(), 'merge.json'), 'w') as json_file:
            json.dump(kmers, json_file, sort_keys = True, indent = 4, separators = (',', ': '))
            