
# RNA-seq Analysis Using R: DESeq2 Workflow

This notebook demonstrates a step-by-step workflow for analyzing RNA-seq data using DESeq2. The analysis includes preprocessing the expression matrix, handling metadata, and identifying differentially expressed genes (DEGs).

---

## Libraries and Data Import
Load necessary libraries and read in the data.

```R
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)

# Raw counts data
rawCounts <- read.delim("GSE275290_raw_counts.tsv")
```

---

## Mapping Gene IDs to Symbols
Convert gene IDs to gene symbols using `org.Hs.eg.db`.

```R
gene_ids = as.character(rawCounts$GeneID)

gene_symbols = mapIds(
  org.Hs.eg.db, 
  keys = gene_ids, 
  keytype = "ENTREZID", 
  column = "SYMBOL", 
  multiVals = "first"
)

gene_df = data.frame(
  GeneID = names(gene_symbols), 
  GeneSymbol = as.vector(gene_symbols), 
  stringsAsFactors = FALSE
)

data = merge(rawCounts, gene_df, by = "GeneID", all.x = TRUE)

data = data %>% 
  dplyr::select(GeneSymbol, everything(), -GeneID)
```

---

## Preparing Metadata
Clean and process metadata.

```R
Allmetadata <- read.csv("Phenotable.csv")

metadata = Allmetadata %>% 
  dplyr::select(Sample.Name, treatment) %>% 
  rename(SampleID = Sample.Name, Condition = treatment)

metadata = metadata %>% 
  mutate(Condition = case_when(
    Condition == "FBZ for 48 hours" ~ "FBZ", 
    Condition == "DMSO; Vehicle" ~ "Control", 
    TRUE ~ Condition
  ))
```

---

## Preprocessing Expression Matrix
### Removing Missing Values
```R
sum(is.na(data))
data = na.omit(data)
```

### Handling Duplicates
```R
sum(duplicated(data$GeneSymbol))

exp.data = data %>% 
  dplyr::select(-GeneSymbol)

exp.data = apply(exp.data, 2, as.numeric)
exp.data.agg = aggregate(exp.data, by = list(data$GeneSymbol), FUN = mean)
rownames(exp.data.agg) = exp.data.agg$Group.1

exp.data.agg = exp.data.agg[-1]
exp.data.agg = exp.data.agg %>% 
  dplyr::select(-Group.1)
```

### Removing Zero-Variance Genes
```R
varRow = apply(exp.data.agg, 1, var, na.rm = TRUE)
constRow = (varRow == 0 | is.na(varRow))
exp.data.agg = exp.data.agg[!constRow, ]
```

---

## Differential Expression Analysis with DESeq2

### Preparing Data
Ensure metadata and expression data match.
```R
exp = exp.data.agg
all(colnames(exp) %in% metadata$SampleID)

metadata$Condition = factor(metadata$Condition, levels = c("Control", "FBZ"))
exp = exp[, metadata$SampleID]
exp = round(exp)
```

### Creating DESeq Dataset
```R
dds = DESeqDataSetFromMatrix(
  countData = exp, 
  colData = metadata, 
  design = ~ Condition
)

dds = dds[rowSums(counts(dds)) >= 10, ]
```

### Running DESeq2
```R
dds.run = DESeq(dds)
res = results(dds.run, contrast = c("Condition", "FBZ", "Control"), alpha = 0.05)
res = res[complete.cases(res), ]
```

---

## Results and Visualization
Summarize and visualize results.

```R
summary(res)
plotMA(res)

res.df = as.data.frame(res)
```

---

## Conclusion
This notebook walks through RNA-seq analysis using DESeq2, from preprocessing raw counts to identifying DEGs. The workflow ensures clean data and accurate results, making it a robust approach for RNA-seq analysis.
