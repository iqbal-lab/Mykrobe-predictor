file=$1

readlink -f $file > /tmp/list.txt
echo /tmp/list.txt > /tmp/clist.txt
readlink -f prelim_panel_k15.fasta > /tmp/fasta_list.txt

/home/phelimb/git/cortex/bin/cortex_var_31_c1 --se_list  /tmp/list.txt --kmer_size 15 \
--mem_height 20 --mem_width 250 --align /tmp/fasta_list.txt,no \
--align_input_format LIST_OF_FASTA --max_read_len 10000