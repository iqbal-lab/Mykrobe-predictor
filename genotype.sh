echo "Counting KMERS" $2 $3
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -o "$1"_"$4".jf -F 2 -m $4 -c 3 -t 40 -s 10000000 <(zcat $2) <(zcat $3) 
echo "Comparing to PANEL" $1
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query "$1"_"$4".jf -s panel_k"$4".kmers | ./genotype.py $1 $4
