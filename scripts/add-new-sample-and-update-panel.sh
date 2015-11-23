comid=$1
mkdir -p /data2/users/phelim/ana/tb/atlas/results/"$comid"
## Add variants
vcf=`ls  /data2/users/phelim/ana/tb/cortex/results/"$comid"/vcfs/*_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf`
~/git/atlas-core/bin/add-new-variants-to-database.py $vcf --sample $comid
## Genotype
ctx=`ls /data2/users/phelim/ana/tb/cortex/results/"$comid"/binaries/uncleaned/31/*.ctx`
~/git/atlas-core/bin/genotype.sh $comid $ctx tb 31 && ./bin/genotype.py  $comid /tmp/"$comid"_tb_31.fasta.colour_covgs
## Sensitivity
~/git/atlas-core/bin/sensitivity.py $comid tb 31 > /data2/users/phelim/ana/tb/atlas/results/"$comid"/sensitivity_before.txt
## Update Panel
~/git/atlas-core/bin/make-variant-panel.py /data2/users/phelim/ana/tb/stampy/ref/NC_000962.2.fasta 
## Regenotype
~/git/atlas-core/bin/genotype.sh $comid $ctx tb 31 && ./bin/genotype.py  $comid /tmp/"$comid"_tb_31.fasta.colour_covgs
## New sensitivity
~/git/atlas-core/bin/sensitivity.py $comid tb 31 > /data2/users/phelim/ana/tb/atlas/results/"$comid"/sensitivity_after.txt
