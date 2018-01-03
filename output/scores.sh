for var in $(ls *.bnd); do
	echo $var >> scores.bnd
	tail -n 2 "./$var" >> scores.bnd
done
