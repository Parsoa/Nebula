#!/bin/bash
cd /share/hormozdiarilab/Codes/NebulousSerendipity
cd src/python
#python create.py
dd if=/dev/zero of=output.json  bs=1M  count=24
dd if=/dev/zero of=/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/hormozdiarilab/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/hormozdiarilab/Codes/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/hormozdiarilab/Codes/NebulousSerendipity/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/hormozdiarilab/Codes/NebulousSerendipity/src/output.json  bs=1M  count=24
dd if=/dev/zero of=/share/hormozdiarilab/Codes/NebulousSerendipity/src/python/output.json  bs=1M  count=24
