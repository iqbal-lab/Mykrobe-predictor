echo "Counting KMERS"
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -m 15 -c 3 -t 50 -s 10000000 <(zcat $1)
echo "Comparing to PANEL"
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query mer_counts.jf -s panel_k15.kmers | ./genotype.py
