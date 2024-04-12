# Multifactorial DESeq2 analysis

### Set up
```
RSCRIPTS=/Users/osipova/Documents/Scripts/rnaseq_tools/

info=ALL_info.tsv
QUANTDIR=Kallisto_quant_ALL/
QUANTFILE=validated.abundance.tsv
QUANTPATH=/projects/project-eosipova/Bird_transcriptomics/

### which tissues will be considered as factors
tissues=liver,pectoralis


### set 1
sp1=Annas_hummingbird
sp2=common_swift

### set 2
sp1=New_Holland_honeyeater
sp2=zebra_finch

### set 3
sp1=rainbow_lorikeet
sp2=cockatiel
```

### run multifactorial DESeq2 analysis
```
Rscript $RSCRIPTS/deseq2_expression_analysis.R \
  -multi \
  -w $(pwd) \
  -d $QUANTDIR \
  -c same_len.gene_${QUANTFILE} \
  -m $info \
  -t $tissues \
  -s1 $sp1 \
  -s2 $sp2 \
  -od multi_deseq2_res.$sp1.tsv  \
  -on norm_counts.multi_deseq2_res.$sp1.tsv;
```


### Get genes of interest -> run GO enrichment analysis
```
for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	ref=ref.multi_deseq2_res.$sp.tsv; \
	interaction=interaction.multi_deseq2_res.$sp.tsv; \
	awk '($7>=0.01 && $3<=2 && $3>=-2){print $1}' $ref > nonsign.${ref%tsv}lst; \

	genes=$sp.multi_deseq2_res.lst; \
	intersect_multiple_files.py -f <(awk '($7<0.05 && $3> 1 ){print $1}' $interaction | grep -v baseMean ) nonsign.${ref%tsv}lst > up.$genes; \
	intersect_multiple_files.py -f <(awk '($7<0.05 && $3<-1 ){print $1}' $interaction | grep -v baseMean ) nonsign.${ref%tsv}lst > down.$genes; \
done
```

### Alternatively (using it now)
```
for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	ref=ref.multi_deseq2_res.$sp.tsv; \
	awk '($7<0.01 && ($3<=-2 || $3>=2)){print $1}' $ref > sign.${ref%tsv}lst; \
done

for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	interaction=interaction.multi_deseq2_res.$sp.tsv; \
	genes=$sp.multi_deseq2_res.lst; \
	filter_annotation_with_list.py -b -c 1 -l sign.ref.multi_deseq2_res.in_any.lst -a <(awk '($7<0.05 && $3>1 ){print $1}' $interaction | grep -v baseMean ) > up.$genes; \
	filter_annotation_with_list.py -b -c 1 -l sign.ref.multi_deseq2_res.in_any.lst -a <(awk '($7<0.05 && $3<-1 ){print $1}' $interaction | grep -v baseMean ) > down.$genes; \
done
```

```
hg38_dict=~/Documents/LabDocs/NectarivoryProject/absrel/absrel_analysis_2024/galGal6_gene.hg38_gene_symbol.tsv

for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	for d in up down; \
	do \
		echo $d; \
		genes=$sp.multi_deseq2_res.lst; \
		renameToHLscaffolds.py -c 1 -a $d.$genes -d <(sed 's/\t/,/' $hg38_dict) > hg38.$d.$genes; \ 
		goenrich_genelist.R -w $(pwd) -g hg38.$d.$genes -u ../background_genes.all_tissues.txt -o goenrich.hg38.$d.${genes%lst}tsv; \
	done; \
done
```


### get rank3 sets
```
for d in up down; \
do \
	for g in $(cut -f1 goenrich.hg38.down.*multi_deseq2_res.tsv | g -v ^ID | s | uniq -c | awk '$1>=3{print $2}'); \
	do \
		grep $g goenrich.hg38.$d.Annas_hummingbird.multi_deseq2_res.tsv ; \
	done | cut -f1,2 > 3_way_convergent_terms.multi_deseq2.${d}_pectoralis.tsv' \
done
```


### exclude children
```
GOOBO=/Users/osipova/Documents/LabDocs/GO_terms_genes/go.obo
golist=$(cut -f2 $f | tail +2 | tr '\n' ',');
f=3_way_convergent_terms.multi_deseq2.${d}_pectoralis.tsv

for g in $(cut -f2 $f | tail +2); \
do \
	get_go_children.py -f $GOOBO -go $g -l $golist ; \
done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; 

filter_annotation_with_list.py -b -c 2 -a $f -l to_exclude_go.lst > noChildren.$f;
```

