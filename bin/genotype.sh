# rm *.fa
# rm *.kmer
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%70000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_list_k31.fasta &
# awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%400000==0){file=sprintf("panel_k31_%d.kmer",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < panel_list_k31.kmers

#echo "Counting KMERS" $2 $3
#/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish count -o /tmp/"$1"_"$4"_"$5".jf -m $5 -c 3 -t 2 -s 10000000 <(zcat $2)
echo "Comparing to PANEL" $1
#/data2/users/phelim/tools/jellyfish-2.2.0/bin/jellyfish query /tmp/"$1"_"$4"_"$5".jf -s panel_tb_k31.kmers > /tmp/"$1"_"$4"_"$5".count
/home/phelimb/git/atlas-core/bin/genotype.py $1 $4 $5 /tmp/"$1"_"$4"_"$5".count --all #&& rm /tmp/"$1"_"$4"_"$5".jf && rm /tmp/"$1"_"$4"_"$5".count 