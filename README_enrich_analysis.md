# Funtional enrichment analysis with ClusterProfiler

### Correct gene names to hg38 for the enrichment analysis
```
hg38_dict=~/Documents/LabDocs/NectarivoryProject/absrel/absrel_analysis_2024/galGal6_gene.hg38_gene_symbol.tsv
```


## enrichGO in DESeq2 pairwise
```
for d in up down; \
do \
  for t in pectoralis liver heart duodenum palate;  \
  do  \
    for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
    do \
      renameToHLscaffolds.py -c 1 -a  March_2024.$d.$sp.$t.top_0.05.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.$d.$sp.$t.top_0.05.txt; \
      goenrich_genelist.R -w $(pwd) -g hg38.$d.$sp.$t.top_0.05.txt -u ../background_genes.$t.txt -o ClusterProfiler/March_2024.hg38.goenrich.$d.$sp.$t.top_0.05.tsv; \
      echo -e "$d,$sp,$t"; \
    done; \
  done; \
done
```


## remove children in rank 2 sets
```
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo
c=1

for f in $(ls rank2.goenrich.*tsv); \
do \
  echo $f; \
  golist=$(cut -f${c} $f | tail -n +2|  tr '\n' ','); \
  for g in $(echo $golist | tr ',' '\n');  \
  do \
    get_go_children.py -f $GOOBO -go $g -l $golist ; \
  done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; \
  wc -l to_exclude_go.lst; \
  filter_annotation_with_list.py -b -c 2 -a $f -l to_exclude_go.lst > noChildren.$f; \
done;
```

## enrichGO in DESeq2 pairwise genes 2-way convergent 
```
for d in up down; \
do \
  for t in pectoralis liver heart duodenum palate; \
  do \
    do cat $d.*.$t.top_0.05.txt | s | uniq -c | awk '$1>=2{print $2}' > $d.rank2.$t.top_0.05.txt ; \
    renameToHLscaffolds.py -c 1 -a $d.rank2.$t.top_0.05.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.$d.rank2.$t.top_0.05.txt; \
    goenrich_genelist.R -w $(pwd) -g hg38.$d.rank2.$t.top_0.05.txt -u ../background_genes.$t.txt -o ClusterProfiler/hg38.goenrich.$d.rank2.$t.top_0.05.txt; \
    echo -e "$d,$sp,$t"; \
  done; \
done
```


## enrichGO in INTRAspecies
### nectar
```
for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
  echo $sp; \
  for t in pectoralis liver; \
  do \
    echo $t; \
    
    renameToHLscaffolds.py -c 1 -a up_${t}_genes_${sp}_not_nonnectar.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.up_${t}_genes_${sp}_not_nonnectar.txt; \
    
    goenrich_genelist.R -w $(pwd) -g hg38.up_${t}_genes_${sp}_not_nonnectar.txt -o ClusterProfiler/hg38.goenrich.up_${t}_genes_${sp}_not_nonnectar.txt; \
  
  done; \
done
```

### nonnectar
```
for sp in common_swift zebra_finch cockatiel; \
do \
  echo $sp; \
  for t in pectoralis liver; \
  do \
    echo $t; \
    
    renameToHLscaffolds.py -c 1 -a up_${t}_genes_${sp}_not_nectar.txt -d  <(sed 's/\t/,/' $hg38_dict) > hg38.up_${t}_genes_${sp}_not_nectar.txt; \

    goenrich_genelist.R -w $(pwd) -g hg38.up_${t}_genes_${sp}_not_nectar.txt -o ClusterProfiler/hg38.goenrich.up_${t}_genes_${sp}_not_nectar.txt; \
  
  done; \
done

```
```
pairs = Annas_hummingbird.common_swift,
         New_Holland_honeyeater.zebra_finch,
         rainbow_lorikeet.cockatiel
```
