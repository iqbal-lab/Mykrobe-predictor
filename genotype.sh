
# rm *.fa
# rm *.kmer
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%70000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_k31.fasta &
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%400000==0){file=sprintf("panel_k31_%d.kmer",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_k31.kmers


echo "Counting KMERS" $2 $3
/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -o /tmp/"$1"_"$4".jf -F 2 -m $4 -c 3 -t 40 -s 10000000 <(zcat $2) <(zcat $3) 
echo "Comparing to PANEL" $1
parallel --gnu -j 13 /data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query /tmp/"$1"_"$4".jf -s {} ::: *kmer | ./genotype.py $1 $4
rm /tmp/*.jf &
