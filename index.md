[Zhanglab@Columbia](https://hanruizhang.github.io/zhanglab/) by [Hanrui Zhang](https://github.com/hanruizhang) [2019-07-01]

The material is modified from the [CBW](https://bioinformatics.ca/) workshop on [pathway and network analysis](https://bioinformaticsdotca.github.io/Pathways_2019).   

## Summary of the Workflow
1. **Process Data**: Derive a list of differentially expressed genes from RNA sequencing data. ([our workflow](https://hanruizhang.github.io/RNAseq-analysis-workflow/))
2. **Identify Pathways**: Identify enriched pathways using over-representation analysis or Gene Set Enrichment Analysis.
3. **Visualize**: Create an Enrichment Map displaying the landscape of pathways.
4. **Build the network**: [ReactomeFI](https://reactome.org/tools/reactome-fiviz#Download_and_Launch_ReactomeFIViz) - investigate and visualize functional interaction among genes in hit pathways.
5. **Predict gene function**: [GeneMANIA](https://genemania.org/) - predict the function of a gene or gene set.
6. **Discover the Regulons**: [iRegulon](http://iregulon.aertslab.org/) - sequence based discovery of the TF, the targets and the motifs/tracks from a set of genes.


## 1. Laptop set-up instruction
[https://bioinformaticsdotca.github.io/Pathways_laptop_setup_instructions](https://bioinformaticsdotca.github.io/Pathways_laptop_setup_instructions)

## 2. Reading materials and references
[https://bioinformaticsdotca.github.io/Pathways_2019_prework](https://bioinformaticsdotca.github.io/Pathways_2019_prework)

Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap
Reimand J, Isserlin R, Voisin V, Kucera M, Tannus-Lopes C, Rostamianfar A, Wadi L, Meyer M, Wong J, Xu C, Merico D, Bader GD. Nat Protoc. 2019 Feb;14(2):482-517. [PubMed Abstract](https://www.ncbi.nlm.nih.gov/pubmed/30664679). [Full-text can be downloaded here](http://baderlab.org/Publications#EM_2019).

## 3. Over-representation analysis and enrichment analysis
### 3.1 [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) for over-representation analysis: Using two lists of genes as the inputs

* Foreground genes should be the differentially expressed genes using different FC and FDR cutoff, e.g. Log2(FC)>0.585 & FDR<0.05.
* Background genes can be the "expressed" genes.
* Advanced option: Use FDR and adjust threshold as needed.
* Notes:
	* Data sources: Start with GO_Biological Processes and check "no electronic GO annotations". Can also include KEGG, Reactiome etc. 
	* Adjust term size to exclude general terms: the default is 10,000 and it is generally good to change to 1,000. If there are still a lot of results can reduce further. 
	* If there are ambiguous IDs, choose the one with the most GO terms, or the first one on the list if all are the same.  
	* The GEM is the generic enrichment file and it is formatted in a way that Enrichment map Cytoscape app can recognize.  It is missing some of the info that is found in the csv but you can use it directly with the Cytoscape app.   
	* The query URL results are not permenant   
	* When you have the large term sizes you will get more general terms coming up. With the previous version of g:profiler you were able to specify the min and max geneset size.  We used to recommend min of 3 and max of 300.  Unfortunately with the latest release of g:profiler you are not able to filter prior to searching by these thresholds.      
* [g:Convert](https://biit.cs.ut.ee/gprofiler/convert): 
	* Target name space: 
		* ENTREZGENE: Entrez gene symbol
		* ENTREZGENE_ACC: Entrez gene ID (unique, and recommended)
* [g:Orth](https://biit.cs.ut.ee/gprofiler/orth): Convert mouse and human ortholog
* [g:SNPense](https://biit.cs.ut.ee/gprofiler/snpense): Map human SNPs to genes

### 3.2 [GSEA](http://software.broadinstitute.org/gsea/index.jsp): Use all genes with fold change and p-value (recommended)
* GSEA analysis needs two files
	* A rank file (.rnk): using the codes below (and modify as needed) to generate a txt file then add ".rnk" to the end of the file name

```
# To prepare .rnk file
## Add a rank column using the following calculation, make sure to use p value, not padj.   
df$rank = -log10(df$pvalue) * sign(df$log2FoldChange)
## Order the rank by descending
df1 = df[order(df$rank, decreasing = TRUE),]
## Remove duplicate and subset GeneName and rank column
df1[which(duplicated(df1$SYMBOL)),]
df2 = df1[,c(9,10)]
df3 = (subset(df2, SYMBOL!="..." & SYMBOL!="..." & SYMBOL!="..."))
## Write to a .rnk file
write.table(df3, file=".txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

```

       
* A pathway definition file (.gmt): 
	* Can be downloaded from [http://baderlab.org/GeneSets](downloaded from http://baderlab.org/GeneSets).
	* Recommended file: **Human_ GOBP_ All_ Pathways_ no_ GO_ iea_ {Date}_{ID}.gmt**
* Number of permutation: always use 1000.
* Important tips: **set Enrichment Statistics to p2 if you want to add more weight on the most top up-regulated and top down-regulated.**. Can try P1.5 to see the difference in terms of enrichment plot.
* Click on "Success" to launch results.  

## 4. Network Visualization and Analysis with Cytoscape - Enrichment Map

* Add an expression table: The rationale for using the expression table is that when you use it you can see the more fine grain details of the sample data. Sometimes individual samples could have stronger signals than the others or it could be equal across the board.  It is just an additional feature that we have found people like to take advantage of. The gene names in the expression table and the rank file should match (they do not have to be in the same order just the identifiers need to match).
* Edge cutoff: we would like to keep it low unless the data look like a hairball.
* Enrichment map - Style - Chart data: Color by FDR or by datasets
* Circle represents number of genes in each enriched pathway; edge represents the overlap of genes between pathway; color of the circle can be adjusted to represent Q value.

#### * EnrichmentMap (GSEA and g:Profiler) step by step protocol: [https://www.biorxiv.org/content/early/2017/12/12/232835](https://www.biorxiv.org/content/early/2017/12/12/232835) and [https://github.com/BaderLab/Cytoscape_workflows](https://github.com/BaderLab/Cytoscape_workflows)


## Additional information


* Use [Biostars](https://www.biostars.org/) to post questions and search answers 



		

_This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license._

