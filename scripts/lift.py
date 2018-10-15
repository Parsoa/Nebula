#!/usr/bin/python
import os
import sys
import subprocess

if __name__ == '__main__':
    print os.getcwd()
    bed = sys.argv[1]
    name = bed[:bed.find('.bed')]
    FNULL = open(os.devnull, 'w')
    with open(name + '.convert.bed', 'w') as c:
        with open(bed) as f:
            for line in f.readlines():
                tokens = line.split()
                name = tokens[0] + '_' + tokens[1] + '_' + tokens[2] + '.bed'
                with open(name, 'w') as t:
                    t.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\t' + tokens[6] + '\n')
                command = '/home/pkhorsand/local/bin/liftOver ' + name + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg19ToHg38.over.chain res.bed un.bed'
                print command
                output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)
                with open('./res.bed', 'r') as r:
                    l = r.readline()
                print(l)
                tokens = tokens[:3] + l.split()
                tokens = [''] + tokens
                s = reduce(lambda x, y: x + '{:12}'.format(y), tokens)
                c.write(s + '\n')
                os.remove('./' + name)
        os.remove('./un.bed')
        os.remove('./res.bed')

