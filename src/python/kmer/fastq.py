
def parse_fastq(fp):
    name = None
    HEADER_LINE = 0
    SEQUENCE_LINE = 1
    THIRD_LINE = 2
    QUALITY_LINE = 3
    state = HEADER_LINE
    # need to skip invalid lines
    for line in fp:
        if state == HEADER_LINE and line[0] == '@':
            state = SEQUENCE_LINE
            name = line[:-1] 
            continue
        if state == SEQUENCE_LINE:
            state = THIRD_LINE
            seq = line[:-1] # ignore the EOL character
            yield seq, name
            continue
        if state == THIRD_LINE:
            state = QUALITY_LINE
            continue
        if state == QUALITY_LINE:
            state = HEADER_LINE
            continue