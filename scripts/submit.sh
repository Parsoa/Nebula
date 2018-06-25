#!/bin/bash
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --job-name="NebulousSerendipity"
#SBATCH --ntasks=48
#SBATCH --time=12:00:00
#SBATCH --mail-user=pkhorsand@ucdavis.edu
#SBATCH --output=/share/hormozdiarilab/Codes/NebulousSerendipity/current.out
#SBATCH --mail-type=ALL
P=$(pwd)
echo $P
cd /share/hormozdiarilab/Codes/NebulousSerendipity
source venv2/bin/activate
cd $P
SRC=$1
shift
/share/hormozdiarilab/Codes/NebulousSerendipity/scripts/$SRC "$@"
