### ALTERNATIVE SPLICING ANALYSIS ###
## SUPPA ##

cd /path

#Generate events
python ./suppa.py generateEvents -i ./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_tr_exon.gtf -o ./code/AS_APA/output/output/APA_AS_corr/ORFanage_events -e SE SS MX RI FL -f ioe
