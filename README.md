##########################################################
#                                                        #
#  Pipeline describing DGE analysis after quantification #
#                                                        #
##########################################################

```
RSCRIPTS=/Users/osipova/Documents/Scripts/rnaseq_tools/

info=ALL_info.tsv
analysis=kallisto
QUANTDIR=Kallisto_quant_ALL/
QUANTFILE=validated.abundance.tsv
QUANTPATH=/projects/project-eosipova/Bird_transcriptomics/
PREFF=''
SUFF=''
```

## Get Salmon/Kallisto quantification results from delta
```
for l in $(cut -f1 $info | grep -v lib);
do
 md ${QUANTDIR}${PREFF}${l}/;
 scp eosipova@10.10.50.202:${QUANTPATH}${QUANTDIR}${PREFF}${l}${SUFF}/${QUANTFILE} ${QUANTDIR}${PREFF}${l}/;
done
```

## Get one2one orthologs from delta
```
for db in HLcalAnn5 HLapuApu1 HLphyNov1 HLtriMol2 HLnymHol2 HLtaeGut4;
do
 one2ones=$db.one2ones.toga.my_anno.dict_GENES.tsv
 scp eosipova@10.10.50.202:/projects/project-eosipova/NectarivoryProject/Genome_annotation_2021/$db/$one2ones .
done
```


## Get gene level quantification from Salmon/Kallisto
```
info=ALL_info.tsv

for db in HLcalAnn5 HLapuApu1 HLphyNov1 HLtriMol2 HLnymHol2 HLtaeGut4;
do
 for l in $(grep $db $info | cut -f1);
 do
  Rscript $RSCRIPTS/get_gene_level_count.R  \
  ${QUANTDIR}${PREFF}$l/${QUANTFILE}  \
  $db.one2ones.toga.my_anno.dict_GENES.tsv  \
  $analysis \
  ${QUANTDIR}${PREFF}$l/gene_${QUANTFILE};
 done;
done
```



## Get genes that don't have much variation in effective length for each pair of species; also that are one2ones
```
### set 1
db1=HLcalAnn5
db2=HLapuApu1
l1=L73700
l2=L73723
sp1=Annas_hummingbird
sp2=common_swift

### set 2
db1=HLphyNov1
db2=HLtaeGut4
l1=L75541
l2=L73702
sp1=New_Holland_honeyeater
sp2=zebra_finch

### set 3
db1=HLtriMol2
db2=HLnymHol2
l1=L75544
l2=L75547
sp1=rainbow_lorikeet
sp2=cockatiel
```

## Get counts of one2ones in pairs 
```
one2ones=one2ones_${db1}_${db2}.lst

db1_ones=$db1.one2ones.toga.my_anno.dict_GENES.tsv
db2_ones=$db2.one2ones.toga.my_anno.dict_GENES.tsv
intersect_multiple_files.py -f <(cut -f2 $db1_ones) <(cut -f2 $db2_ones) | sort -u > $one2ones


intersect_multiple_files.py -f \
<(
  compare_values_two_anno.py \
  -d 100 \
  -a1 <(get_gene_stats_from_anno.py -c 3 -s all -i <(awk '{print $2"\t"$1}' $db1_ones) -a ${QUANTDIR}${l1}/${QUANTFILE}) \
  -a2 <(get_gene_stats_from_anno.py -c 3 -s all -i <(awk '{print $2"\t"$1}' $db2_ones) -a ${QUANTDIR}${l2}/${QUANTFILE}) \
  ) \
  $one2ones > same_len.$one2ones
```

## Filter quant files for genes with little variation in eff length and are one2ones in both species 
```
for l in $(grep "$db1\|$db2" $info | cut -f1);
do
 echo $l;
 filter_annotation_with_list.py \
  -l same_len.$one2ones \
  -a ${QUANTDIR}${PREFF}$l/gene_${QUANTFILE} \
  -c 1 > ${QUANTDIR}${PREFF}$l/same_len.gene_${QUANTFILE};
done
```


## Run DESeq2 analysis with R!

### NB!! sp1 = TEST species
### NB!! sp2 = REFERENCE species
### NB!! header in ALL_info:   lib sample  species type  db 

```
for tissue in gonads heart fat liver pectoralis pancreas kidney duodenum proventriculus tongue palate;
do
echo $tissue;
Rscript $RSCRIPTS/deseq2_expression_analysis.R \
  -w $(pwd) \
  -d $QUANTDIR \
  -c same_len.gene_${QUANTFILE} \
  -m $info \
  -t $tissue \
  -s1 $sp1 \
  -s2 $sp2 \
  -od deseq2_res.$tissue.$sp1.$sp2.tsv  \
  -on norm_counts.deseq2_res.$tissue.$sp1.$sp2.tsv;
done
```

## Fix headings
```
for f in $(ls deseq2_res.*.tsv);
do
 echo $f;
 sed 's/^baseMean/gene\tbaseMean/' $f > file;
 mv file $f;
done
```


## Run DESeq2 analysis for INTRAspecies, pectoralis VS liver
### NB!! header in ALL_info:   lib samples  type  species  db 
```
for sp in Annas_hummingbird common_swift cockatiel New_Holland_honeyeater zebra_finch rainbow_lorikeet;
do
 Rscript $RSCRIPTS/deseq2_expression_analysis.R \
 -w $(pwd) \
 -d Kallisto_quant_ALL/ \
 -c same_len.gene_${QUANTFILE} \
 -m ALL_info.tsv \
 -t $sp \
 -s1 pectoralis \
 -s2 liver \
 -od deseq2_res.$sp.pectoralis.liver.tsv \
 -on norm_counts.deseq2_res.$sp.pectoralis.liver.tsv;
done
```




