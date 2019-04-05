#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd src/python
python -m kmer.main "$@"
