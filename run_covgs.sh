#./get_unique_kmers.py > prelim_panel_k15_unique_kmers.fa

echo "Counting KMERS"
# /data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -F 2 -m 15 -c 3 -t 50 -s 10000000  /data2/users/phelim/data/staph/val/C00012831_R00000022.fq_1.fastq /data2/users/phelim/data/staph/val/C00012831_R00000022.fq_2.fastq
echo "Comparing to PANEL"
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query mer_counts.jf -s prelim_panel_k15_unique_kmers.text | ./genotype.py
