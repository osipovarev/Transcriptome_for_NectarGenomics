## Correct gene names to hg38 for the enrichment analysis

hg38_dict=~/Documents/LabDocs/NectarivoryProject/absrel/absrel_analysis_2024/galGal6_gene.hg38_gene_symbol.tsv

#########################################
#                                       #
#           ClusterProfiler             #
#                                       #
#########################################


## enrichGO in DESeq2 pairwise

for d in up down; \
do \
  for t in pectoralis liver heart duodenum palate;  \
  do  \
    for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
    do \
      renameToHLscaffolds.py -c 1 -a  March_2024.$d.$sp.$t.top_0.05.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.$d.$sp.$t.top_0.05.txt; \
      goenrich_genelist.R -w $(pwd) -g hg38.$d.$sp.$t.top_0.05.txt -u ../background_genes.$t.txt -o ClusterProfiler/March_2024.hg38.goenrich.$d.$sp.$t.top_0.05.tsv; echo -e "$d,$sp,$t"; \
    done; \
  done; \
done


## remove children in rank 2 
for f in $(ls rank2.goenrich.*tsv); do echo $f; golist=$(cut -f2 $f | tail +2 | tr '\n' ','); for g in $(cut -f2 $f | tail +2); do get_go_children.py -f $GOOBO -go $g -l $golist ; done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; wc -l to_exclude_go.lst; filter_annotation_with_list.py -b -c 2 -a $f -l to_exclude_go.lst > noChildren.$f; done;


## enrichGO in Intraspecies

for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; do echo $sp; for t in pectoralis liver; do echo $t; renameToHLscaffolds.py -c 1 -a up_${t}_genes_${sp}_not_nonnectar.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.up_${t}_genes_${sp}_not_nonnectar.txt; goenrich_genelist.R -w $(pwd) -g hg38.up_${t}_genes_${sp}_not_nonnectar.txt -o ClusterProfiler/hg38.goenrich.up_${t}_genes_${sp}_not_nonnectar.txt; done; done


for sp in common_swift zebra_finch cockatiel; do echo $sp; for t in pectoralis liver; do echo $t; renameToHLscaffolds.py -c 1 -a up_${t}_genes_${sp}_not_nectar.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.up_${t}_genes_${sp}_not_nectar.txt; goenrich_genelist.R -w $(pwd) -g hg38.up_${t}_genes_${sp}_not_nectar.txt -o ClusterProfiler/hg38.goenrich.up_${t}_genes_${sp}_not_nectar.txt; done; done



pairs = Annas_hummingbird.common_swift,
         New_Holland_honeyeater.zebra_finch,
         rainbow_lorikeet.cockatiel


#########################################
#                                       #
#                 Metascape             #
#                                       #
#########################################


## Pairwise DESeq2
d=up
d=down

for d in up down; \
do \
  for t in pectoralis liver heart duodenum; \
  do \
    for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
    do \
      printf $sp"\t"; cat $d.$sp.$t.top_0.05.txt | tr '\n' ','; echo; \
    done > for_metascape.$d.3lists.$t.txt; \
  done; \
done



## Remove children GO terms
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo

metaf=metascape.duodenum_down.rank3.tsv
golist=$(cut -f4 $metaf  | tail +2 | tr '\n' ',')

for g in $(cut -f4 $metaf | tail +2); do get_go_children.py -f $GOOBO -go $g -l $golist ; done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst

filter_annotation_with_list.py -b -c 4 -a $metaf -l to_exclude_go.lst > noChildren.$metaf



### Metascape enrichment 2024.01.17

for t in pectoralis liver heart duodenum; \
do \
 echo $t; \
 for d in up down; \
 do \
  f=checked.${d}_${t}.metascape_GO_enrichments.tsv; \
  golist=$(awk -F"\t" '$11 != ""{print }' $f | cut -f4 | tr '\n' ','); \
  for g in $(awk -F"\t" '$11 != ""{print }' $f | cut -f4); \
  do \
   get_go_children.py -f $GOOBO -go $g -l $golist; done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; \
   filter_annotation_with_list.py -c 4 -b -l to_exclude_go.lst -a <(awk -F"\t" '$11 != ""{print }' $f) > noChildren.$f; \
  done; \
 done




### Intraspecies DESeq2
for t in pectoralis liver; \
do 
  for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet;
  do
    printf $sp"\t";
    cat up_${t}_genes_${sp}_not_nonnectar.txt | tr '\n' ','; echo; \
  done > for_metascape_3lists.up_${t}.not_nonnectar.txt ;
done

for t in pectoralis liver; \
do \
  for sp in common_swift zebra_finch cockatiel; \
  do \
    printf $sp"\t"; \
    cat up_${t}_genes_${sp}_not_nectar.txt | tr '\n' ','; \
    echo; \
  done > for_metascape_3lists.up_${t}.not_nectar.txt ; \
done




#########################################
#				                    #
#              ENRICHR	           	#
#                                       #
#########################################

## Get GO terms for comparison
RNADIR=/Users/osipova/Documents/LabDocs/Bird_transcriptomics/Transcriptome_for_NectarGenomics/Enrichment_tests/
sp=Annas_hummingbird
sp=New_Holland_honeyeater
sp=rainbow_lorikeet
t=liver
#t=pectoralis
t=duodenum
t=heart


sublime ../DESeq2_results_Kallisto/up.$sp.$t.top_0.02.txt ../DESeq2_results_Kallisto/down.$sp.$t.top_0.02.txt
### Run enrichment analysis with Enrichr
# go to: https://maayanlab.cloud/Enrichr/
# submit gene lists!

up=$RNADIR/up.$t.$sp.tsv
down=$RNADIR/down.$t.$sp.tsv
mv ~/Downloads/BioPlanet_2019_table.txt $up
mv ~/Downloads/BioPlanet_2019_table.txt $down

cat $up $down | awk -F"\t" '{print $1"\t"$4"\t"$7}'  | awk -F"\t" '$2<.05 {print $1}' > GO_terms_to_look.$t.$sp.lst
# OR:
cat $up $down | awk -F"\t" '{print $1"\t"$3"\t"$7}'  | awk -F"\t" '$2<.01 {print $1}' > GO_terms_to_look.$t.$sp.lst

## Get adjusted p-value and Odds ratio for each term
get_stats_from_enrichr.sh GO_terms_to_look.$t.$sp.lst $up,$down 4 7 > enrichr_results.GO.$t.$sp.tsv
# term  pval1   odds1   pval2   odds2   pval3   odds3   pval4   odds4
# term  up_pval up_odds down_pval down_odds

## Run analysis with jupytyer-notebook:
jupytyer-notebook enrichment_analysis.ipynb

 

