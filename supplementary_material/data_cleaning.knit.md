---
output:
  BiocStyle::html_document
---



# Import Data

We import the raw data from the SummarizedExperiment. It contains only one assay, in a matrix-like object with the genes in the row and the samples in the columns. The data is the counts for gene and sample.


```r
library(SummarizedExperiment)

se <- readRDS(file.path("rawCounts", "seLUSC.rds"))
se
```

```
class: RangedSummarizedExperiment 
dim: 20115 553 
metadata(5): experimentData annotation cancerTypeCode
  cancerTypeDescription objectCreationDate
assays(1): counts
rownames(20115): 1 2 ... 102724473 103091865
rowData names(3): symbol txlen txgc
colnames(553): TCGA.18.3406.01A.01R.0980.07
  TCGA.18.3407.01A.01R.0980.07 ... TCGA.90.7767.11A.01R.2125.07
  TCGA.92.7340.11A.01R.2045.07
colData names(549): type bcr_patient_uuid ...
  lymph_nodes_aortic_pos_by_ihc lymph_nodes_aortic_pos_total
```

The raw table of counts contains RNA-seq data from 20115 genes and 553 samples.

# Extract paired samples

We want to do a paired experiment. For this reason, from the 553 samples, we will keep those that share patient identification, which is in the factor _bcr\_patient\_barcode_.

We select those patient identification that are in normal tissue and in tumor tissue.


```r
normal <- data.frame(colData(se)[colData(se)$type == 'normal',])$bcr_patient_barcode
tumor <- data.frame(colData(se)[colData(se)$type == 'tumor',])$bcr_patient_barcode
length(normal)
```

```
[1] 51
```

```r
length(tumor)
```

```
[1] 502
```

We have 51 samples with normal tissue and 502 samples with tumoral tissue. Let's see which of the normal samples share patient identification with tumor samples.


```r
common_bcr_patient_barcode <- normal[normal %in% tumor]
length(common_bcr_patient_barcode)
```

```
[1] 47
```

47 patient identification are in normal samples and tumor samples. This means that for those patients, there were extraction for normal tissue and for tumoral tissue. 

Now, we filter the SummarizedExperiment object to keep only these patients.


```r
paired_seLUSC <- se[,colData(se)$bcr_patient_barcode %in% common_bcr_patient_barcode]
paired_seLUSC
```

```
class: RangedSummarizedExperiment 
dim: 20115 94 
metadata(5): experimentData annotation cancerTypeCode
  cancerTypeDescription objectCreationDate
assays(1): counts
rownames(20115): 1 2 ... 102724473 103091865
rowData names(3): symbol txlen txgc
colnames(94): TCGA.22.4593.01A.21R.1820.07
  TCGA.22.4609.01A.21R.2125.07 ... TCGA.90.7767.11A.01R.2125.07
  TCGA.92.7340.11A.01R.2045.07
colData names(549): type bcr_patient_uuid ...
  lymph_nodes_aortic_pos_by_ihc lymph_nodes_aortic_pos_total
```

Now, there are 20115 genes, the same as before because we have not filtered genes, and 94 samples, which are 47 normal samples and 47 tumoral samples, with the same patient identificacion.


```r
length(unique(paired_seLUSC$bcr_patient_barcode))
```

```
[1] 47
```

```r
length(levels(paired_seLUSC$bcr_patient_barcode))
```

```
[1] 495
```

Analyzing the length of the array that contains the patient identification and its levels, we see that, although we have discarded the patients that did not pass the filter of the paired design, the levels are still with all the patient identifications, 495.

When doing differential expression analysis, the levels are taken to build the design matrix. For this reason, as the leves are still the one in the raw data, we have to take out the levels that are not in the array of patient identifications.


```r
paired_seLUSC$bcr_patient_barcode <- droplevels(paired_seLUSC$bcr_patient_barcode)
length(unique(paired_seLUSC$bcr_patient_barcode))
```

```
[1] 47
```

```r
length(levels(paired_seLUSC$bcr_patient_barcode))
```

```
[1] 47
```

Now we have 47 patient identifications and also 47 levels. Now we can save this filtered SummarizedExperiment in an RDS file.


```r
saveRDS(paired_seLUSC, file.path("results", "paired_seLUSC.rds"))
```
