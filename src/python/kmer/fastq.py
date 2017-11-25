
def parse_fastq(fp):
    name = None
    skip = False
    for line in fp:
        if skip: # quality line
            skip = False
            continue
        if line[0] == '@': # first line
            name = line[:-1] 
            continue
        if line[0] == '+': # third line
            skip = True
            continue
        seq = line[:-1] # ignore the EOL character
        yield seq, name