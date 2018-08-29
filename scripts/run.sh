
#!/bin/bash
export P=$(pwd)
echo PWD $P
JOB=$(echo $P | awk -F/ '{ print $NF }')
echo JOB $JOB
/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/$JOB.sh $@
