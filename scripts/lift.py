#!/usr/bin/python
import os
import sys
import subprocess

# usage
# lift.py <file.bed> <19|38>
# or
# lift.py chrom begin end

# output
# liftted coordinates, original coordinates, other fielda (if available)

def lift_track():
    with open('tmp.bed', 'w') as b:
        b.write(sys.argv[1] + '\t' + sys.argv[2] + '\t' + sys.argv[3] + '\n')
    command = '/home/pkhorsand/local/bin/liftOver ' + 'tmp.bed' + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg19ToHg38.over.chain res.bed un.bed'
    output = subprocess.call(command, shell = True)
    command = 'cat res.bed'
    output = subprocess.call(command, shell = True)

if __name__ == '__main__':
    print os.getcwd()
    bed = sys.argv[1]
    to = int(sys.argv[2]) if len(sys.argv) == 3 else 38
    if len(sys.argv) == 5:
        to = sys.argv[4]
        lift_track()
        exit()
    name = bed[:bed.find('.bed')]
    FNULL = open(os.devnull, 'w')
    with open(name + '.lift.bed', 'w') as c:
        with open(bed) as f:
            for line in f.readlines():
                tokens = line.split()
                if 'chrom' in tokens[0].lower():
                    c.write('CHROM\tBEGIN\tEND\tCHROM_OLD\tBEGIN_OLD\tEND_OLD\t' + reduce(lambda x, y: x + '\t' + y, tokens[3:]) + '\n')
                    continue
                name = tokens[0] + '_' + tokens[1] + '_' + tokens[2] + '.bed'
                with open(name, 'w') as t:
                    t.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\n')
                if to == 38:
                    command = '/home/pkhorsand/local/bin/liftOver ' + name + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg19ToHg38.over.chain res.bed un.bed'
                else:
                    command = '/home/pkhorsand/local/bin/liftOver ' + name + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg38ToHg19.over.chain res.bed un.bed'
                output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)
                with open('./res.bed', 'r') as r:
                    l = r.readline()
                l = l.split()
                if len(l) >= 3:
                    tokens = l + tokens
                    tokens = [''] + tokens
                    s = reduce(lambda x, y: x + '\t' + y, tokens)
                    c.write(s.strip() + '\n')
                os.remove('./' + name)
        os.remove('./un.bed')
        os.remove('./res.bed')

