# mixed staph
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.staph --list /home/phelimb/git/mykrobe_paper/json_tests/staph/mixture.filelist \
--install_dir /home/phelimb/git/myKrobe-predictor/ > staph/mixed_staph.json & 
#  pure aureus
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.staph --list \
/Net/banyan/data2/users/zam/staph_amr/filelists/C00000784.filelist \
--install_dir /home/phelimb/git/myKrobe-predictor/ > staph/pure_aureus.json &

#  pure pure_non_aureus staph
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.staph --install_dir \
/home/phelimb/git/myKrobe-predictor/ \
--file /Net/banyan/data2/users/phelimb/staph/sra_download/Staphylococcus_epidermidis/SRR221652.fastq.gz >  staph/pure_non_aureus.json &

#  pure non staph
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.staph --install_dir \
/home/phelimb/git/myKrobe-predictor/ \
--file /data2/users/phelim/data/tb/val_2_bams/C00022264.fastq.gz > staph/non_staph.json & 

#NTM_mixed
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.tb --install_dir \
/home/phelimb/git/myKrobe-predictor/ --list \
/data2/users/phelim/ana/tb/species/mixed/mixture_filelists/C00009088_ERR350153/C00009088_ERR350153_40_20.txt > tb/MTBC_NTM_mixed.json & 
#/NON_MYCO
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.tb --install_dir \
/home/phelimb/git/myKrobe-predictor/ --file /Net/banyan/data2/users/phelimb/staph/sra_download/Staphylococcus_epidermidis/SRR221652.fastq.gz > tb/NON_MYCO.json & 
#pure_MTBC
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.tb --install_dir \
/home/phelimb/git/myKrobe-predictor/ --file /data2/users/phelim/data/tb/der/C00008823.fastq.gz > tb/pure_MTBC.json & 
#/pure_ntm
/home/phelimb/git/myKrobe-predictor/bin/Mykrobe.predictor.tb --install_dir \
/home/phelimb/git/myKrobe-predictor/ --file /data2/users/phelim/data/tb/NTMs/C00018214.bam > tb/pure_ntm.json &
wait
