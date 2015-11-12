# rm *.fa
# rm *.kmer
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%70000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_list_k31.fasta &
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%400000==0){file=sprintf("panel_k31_%d.kmer",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_list_k31.kmers


# echo "Counting KMERS" $2 $3
# /data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -o /tmp/"$1"_"$4"_"$5".jf -F 2 -m $5 -c 3 -t 10 -s 10000000 <(zcat $2) <(zcat $3) 
# echo "Comparing to PANEL" $1
parallel --gnu -j 1 /data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query /tmp/"$1"_"$4"_"$5".jf -s {} ::: *kmer > /tmp/"$1"_"$4"_"$5".count
/home/phelimb/git/atlas/bin/genotype.py $1 $4 $5 /tmp/"$1"_"$4"_"$5".count && rm /tmp/"$1"_"$4"_"$5".jf && rm /tmp/"$1"_"$4"_"$5".count 