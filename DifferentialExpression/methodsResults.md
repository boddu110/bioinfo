-----

author: Jyothika Boddu output: html\_document: toc: true toc\_depth: 4
toc\_float: true dev: ‘svg’ md\_document: variant: gfm bibliography:
bibliography.ris nocite:
‘(<span class="citeproc-not-found" data-reference-id="*">**???**</span>)’
— \#\# Differential Expression Analysis

## Methods

The two main steps in performing differential expression analysis are to
estimate the relative abundance of transcripts, and to apply statistical
models to test for differential expression between treatment groups.
Estimating relative abundance is basically determining how many NGS
reads match a given gene within a genome.

Salmon, a tool for quantifying transcript abundance from RNA-seq reads
that is reliable and fast. Salmon is the first transcriptome-wide
quantifier to correct bias for fragment GC content. We conducted salmon
index and quantity in this module. The index is a framework that salmon
use to quasi-map RNA-seq reads during quantification. Salmon conducts
its inference using a descriptive and practical RNA-seq data model that
takes into account experimental properties and limitations commonly
observed in actual RNA-seq results. Next, it includes a collection of
target transcripts from de-novo assembly to km length 25 indexes. Last,
the data were matched and the salmon quant of the reports became
associated with all the samples.

The second method is to use tximport, You will need a table mapping
transcripts to genes to use the Salmon data in tximport. This allows
translate transcript-level abundance projections from specific
quantification pipelines into count-based computational inferencing
engines.

The third method is to use DESeq2, here we compose an R script with the
contents to import the Salmon alignments into DESeq2 and conduct the
study of differential expressions.

## Results

``` r
library(knitr)
annot<- read.csv("deAnnotated.csv",header=T)
kable(annot)
```

|  X | pathway      | ko        | log2FoldChange |      padj | Factor                        | description                               |
| -: | :----------- | :-------- | -------------: | --------: | :---------------------------- | :---------------------------------------- |
|  1 | path:ko03015 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | mRNA surveillance pathway                 |
|  2 | path:ko04071 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Sphingolipid signaling pathway            |
|  3 | path:ko04111 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Cell cycle - yeast                        |
|  4 | path:ko04144 | ko:K17920 |       1.118851 | 0.0000011 | Menthol\_Menthol\_vs\_Control | Endocytosis                               |
|  5 | path:ko04151 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | PI3K-Akt signaling pathway                |
|  6 | path:ko04152 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | AMPK signaling pathway                    |
|  7 | path:ko04261 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Adrenergic signaling in cardiomyocytes    |
|  8 | path:ko04390 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway                   |
|  9 | path:ko04391 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway - fly             |
| 10 | path:ko04530 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Tight junction                            |
| 11 | path:ko04728 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Dopaminergic synapse                      |
| 12 | path:ko05142 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Chagas disease (American trypanosomiasis) |
| 13 | path:ko05160 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Hepatitis C                               |
| 14 | path:ko05165 | ko:K04354 |     \-2.068930 | 0.0273593 | Menthol\_Menthol\_vs\_Control | Human papillomavirus infection            |

## References
