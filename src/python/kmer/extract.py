
chroms = {}

def search(kmer):
    for chrom in chroms:
        f = 0
        while f != -1:
            f = chroms[chrom].find(kmer, f + 1)
            if f != -1:
                print(chrom, f)
                print(chroms[chrom][f -32 : f + 64])

def extract_chromosomes(name, chromosomes):
    sequence = ''
    ref = open(name)
    line = ref.readline().lower().strip()
    found = False
    m = 0
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:].lower()
            if True:#chrom in chromosomes:
                print('extracting ' + chrom)
                while True:
                    line = ref.readline().lower().strip()
                    if line.startswith('>') or len(line) == 0:
                        print(len(sequence), 'bases')
                        yield sequence, chrom
                        sequence = ''
                        found = True
                        m += 1
                        #if m == len(chromosomes):
                        #    return
                        break
                    sequence += line.upper()
        # this is to avoid skipping the last line we read for the previous chromosome (header of next)
        if found:
            found = False
            continue
        line = ref.readline().lower().strip()
        if len(line) == 0 or not line:
            break

def extract_whole_genome(name):
    print('extracting whole genome')
    c = ['chr' + str(x) for x in range(1, 23)]
    c.append('chrx')
    c.append('chry')
    for seq, chrom in extract_chromosomes(name, c):
        chroms[chrom] = seq
    return chroms
