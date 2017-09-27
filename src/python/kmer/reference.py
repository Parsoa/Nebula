import os
import pwd

import pybedtools

class ReferenceGenome:
    class __Impl:
        def __init__(self, path):
            self.path = path
            self.fasta = pybedtools.BedTool(path)

    __instance = None

    def __init__(self, path=None):
        if not ReferenceGenome.__instance:
            ReferenceGenome.__instance = ReferenceGenome.__Impl(path)
    def __getattr__(self, name):
        return getattr(self.__instance, name)
