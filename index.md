---
title: "Pathway and Network Analysis Workflow"
output:
    html_document:
        toc: false
---

[Zhanglab@Columbia](https://hanruizhang.github.io/zhanglab/) by [Hanrui Zhang](https://github.com/hanruizhang) [2019-07-01], updated on [2021-02-13].    

The material is modified from the [CBW](https://www.bioinformatics.ca/) workshop on pathway and network analysis [2020](https://bioinformaticsdotca.github.io/Pathways_2020). 

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
The protocol uses publicly available software packages (GSEA v.3.0 or higher, g:Profiler, Enrichment Map v.3.0 or higher, Cytoscape v.3.6.0 or higher) and custom R scripts that apply publicly available R packages (edgeR, Roast, Limma, Camera). Custom scripts are available in the Supplementary Protocols and at GitHub web sites [https://github.com/BaderLab/Cytoscape_workflows/tree/master/EnrichmentMapPipeline](https://github.com/BaderLab/Cytoscape_workflows/tree/master/EnrichmentMapPipeline) and [https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/index.html](https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/index.html).

## 3. Over-representation analysis and enrichment analysis
### 3.1 [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) for over-representation analysis: Using two lists of genes as the inputs

* **Select organism that matches input query gene list.**
* **Foreground genes**: Should be the differentially expressed genes using different FC and FDR cutoff, e.g. Log2(FC)>1.0 & FDR<0.01. You may download the practice data in .txt file [here](/data/M0-HMDM_M1-HMDM_Increase_FDR0.01_FC2.0.txt). 
	* The list include DE genes upregualted in human monocyte-derived macrophages (HMDM) stimualted with LPS+TFNg, i.e. M(LPS+IFNgamma), vs HMDM without stimulation (M0). The results are obtained by running the RNA-seq workflow also available on the lab web page [here](https://hanruizhang.github.io/RNAseq-analysis-workflow/).
	* It is recommended that for a list that can be ranked, e.g. a list of DE genes, the genes should be ranked and the box of **"Ordered query"** should be checked.
	* The ranking can be based on the same methods for GSEA analysis as described below. 
* **Background genes**: Can select the "Only annotated genes" upon clicking "Advanced Options". Or can use the "expressed" genes and you may download the practice data in .txt file [here](/data/M0-HMDM_Background.txt).
	* The "expressed" genes can be defined by customized criteria. For example:
		* Require that the sum of normalized counts for all samples is 10 or higher as **"expressed"** (this is considered the minimal pre-filtering applied before DESeq2 analysis to keep only rows that have at least 10 reads total).
		`rowSums(counts(dds)) >= 10`
		* Additional filtering can be used, e.g. for at least 3 samples have a count of 10 or higher.
		`rowSums(counts(dds) >= 10) >= 3`
		* For RNA-seq that is considered to be genome-wide coverage, it is usually just fine to use the "Only annotated genes" as the background, but it is better to run the analysis in both ways and compare the results.
* **Statistical threshold**: Click "Advanced option" and select FDR.
* **Data source**: It is recommended to start with "GO: Biological Process" and check "no electronic GO annotations", and "Reactome" for the initial analysis. Then can repeat the analysis with other or all gene sets/pathways included.
* **Run query**: If there are ambiguous IDs, choose the one with the most GO terms, or the first one on the list if all are the same. 
* **Results interpretation**:
	* Adjust term size to exclude general terms: the default is 10,000 and it is generally good to change to 1,000. If there are still a lot of results can reduce further. This is because large pathways are of limited interpretative value, whereas numerous small pathways decrease the statistical power because of excessive multiple testing.
		* **For enrichment map analysis, may try min = 3 and max = 250 to limit the results for more informative map.** 
		* With the previous version of g:profiler you were able to specify the min and max geneset size.  We used to recommend min of 3 and max of 300.  Unfortunately with the latest release of g:profiler you are not able to filter prior to searching by these thresholds.  
		* Keep a note for this filtering strategy, e.g. name the results folder including the min and max number.      
* **Save the results**: The query URL results are not permenant. The GEM is the generic enrichment file and it is formatted in a way that Enrichment map Cytoscape app can recognize.  It is missing some of the info that is found in the csv but you can use it directly with the Cytoscape app.
   
* [g:Convert](https://biit.cs.ut.ee/gprofiler/convert): 
	* Target name space: 
		* ENTREZGENE: Entrez gene symbol
		* ENTREZGENE_ACC: Entrez gene ID (unique, and recommended)
* [g:Orth](https://biit.cs.ut.ee/gprofiler/orth): Convert mouse and human ortholog
* [g:SNPense](https://biit.cs.ut.ee/gprofiler/snpense): Map human SNPs to genes

### 3.2 [GSEA](http://software.broadinstitute.org/gsea/index.jsp): Use all genes with fold change and p-value (recommended)
* Download the GSEA desktop app by clicking [here](http://software.broadinstitute.org/gsea/index.jsp) and then "Download". Login with your email and choose the Mac vs Windows version appropriate.
* GSEA analysis needs two files
	* A rank file (.rnk): using the codes below (and modify as needed) to generate a txt file then add ".rnk" to the end of the file name

```
# To prepare .rnk file
## Read the DESeq2 output
DE <- read.csv("../output/DESeq2.csv"), header = TRUE, sep = ",")

## Filter to remove all the NAs 
DEnoNA <- DE %>% filter(!is.na(SYMBOL) & !is.na(pvalue) & !is.na(padj))

## Filter to remove duplicated SYMBOLS
duplicate <- DEnoNA[which(duplicated(DEnoNA$SYMBOL)),]
duplicate_SYMBOL <- duplicate$SYMBOL
DEfinal <- DEnoNA[!grepl(paste(duplicate_SYMBOL, collapse = "|"), DEnoNA$SYMBOL),]

## Prepare the rnk file and save to the output folder
### Add a rank column using the following calculation, make sure to use p value, not padj.  
DEfinal$rank = -log10(DEfinal$pvalue) * sign(DEfinal$log2FoldChange)
### order by rank and subset SYMBOL and rank column
rnk = DEfinal[order(DEfinal$rank, decreasing = TRUE), 8:9]

### Write to a .rnk file
write.table(rnk, file="../output/DESeq2.rnk"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


```

       
* A pathway definition file (.gmt): 
	* Can be downloaded from [http://baderlab.org/GeneSets](downloaded from http://baderlab.org/GeneSets) and use the "current release". Recommended file: **Human_ GOBP_ All_ Pathways_ no_ GO_ iea_ {Date}_{ID}.gmt**
	* gmt file for mouse databases can also be downloaded from gprofiler, but GSEA results using the mouse databases only has GO ID but not description. 
	* Can also just use the gmt files already available in the GSEA destop app, but can only be for human data (human gene symbol). One can use g:orth to convert mouse gene symbol in the rank file to human orthologs and provide the rnk file for GSEA analysis.
* Number of permutation: always use 1000.
* Name the analysis with the gene sets used, and the min and max cutoff.
* Min = 15, max = 200 to include only smaller gene sets.
* Enrichment statistics: The default is "weighted". Can also **set Enrichment Statistics to p2 if you want to add more weight on the most top up-regulated and top down-regulated.**. Can try P1.5 to see the difference in terms of enrichment plot.
* Click "run" to start the analysis.
* Click on "Success" to launch results.  

## 4. Network Visualization and Analysis with Cytoscape - Enrichment Map
#### Network Visualization and Analysis with Cytoscape: Enrichment Map from g:Profiler results.
https://bioinformaticsdotca.github.io/Pathways_2019_Module3_Lab-EM_GProfiler    

#### Network Visualization and Analysis with Cytoscape: create an enrichment map from GSEA results.
https://bioinformaticsdotca.github.io/Pathways_2019_Module3_Lab-EM_GSEA    

**Notes:**
* Enrichment map documentation: [https://enrichmentmap.readthedocs.io/en/latest/](https://enrichmentmap.readthedocs.io/en/latest/)
* For GSEA results-based analysis, adding an expression table is recommended: The rationale for using the expression table is that when you use it you can see the more fine grain details of the sample data. Sometimes individual samples could have stronger signals than the others or it could be equal across the board.  It is just an additional feature that we have found people like to take advantage of. The gene names in the expression table and the rank file should match (they do not have to be in the same order just the identifiers need to match).
* The parameter for the enrichment map should be set up the way that the map is informative and representative (not like a hairball or miss any important information): The FDR cutoff for the number of nodes can be between 0.001 to 0.01; the edge connectivity should be between sparse to the middle. 
* Once the map is generated, we can further filter by FDR or Edge cutoff, which we would like to keep it low unless the data look like a hairball.
* Enrichment map - Style - Chart data: Color by FDR or by datasets
* Circle represents number of genes in each enriched pathway; edge represents the overlap of genes between pathway; color of the circle can be adjusted to represent Q value.


## 5. Additional information
* Other Bader lab Cytoscape workflows are available at [https://github.com/BaderLab/Cytoscape_workflows](https://github.com/BaderLab/Cytoscape_workflows)
* Use [Biostars](https://www.biostars.org/) to post questions and search answers
* This workflow does not cover the contents below, which can be found through workshop link above.
	* **Build the network** by [ReactomeFI](https://reactome.org/tools/reactome-fiviz#Download_and_Launch_ReactomeFIViz) - investigate and visualize functional interaction among genes in hit pathways.
	* **Predict gene function**: [GeneMANIA](https://genemania.org/) - predict the function of a gene or gene set.
	* **Discover the Regulons**: [iRegulon](http://iregulon.aertslab.org/) - sequence based discovery of the TF, the targets and the motifs/tracks from a set of genes.



		

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

