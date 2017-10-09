from kmer import bed

def main():
    bed.configure()
    bed.read_tracks_from_bed_file()

if __name__ == '__main__' :
    main()