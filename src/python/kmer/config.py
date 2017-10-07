class Configuration:

    class __impl:
        def __init__(self, ksize, khmer_table_size, khmer_num_tables, fastq_file):
            self.ksize = ksize
            self.khmer_table_size = khmer_table_size
            self.khmer_num_tables = khmer_num_tables
            self.fastq_file = fastq_file

        def kmer_size(self):
            return self.ksize

    __instance = None

    def __init__(self, ksize=None, khmer_table_size=None, khmer_num_tables=None, fastq_file=None):
        if Configuration.__instance is None:
            Configuration.__instance = Configuration.__impl(ksize, khmer_table_size, khmer_num_tables, fastq_file)

    def __getattr__(self, attr):
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        return setattr(self.__instance, attr, value)
 