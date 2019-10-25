---
title: "RnaSeqAnalysis"
author: "Luisa Santus"
date: "10/18/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Importing and formatting

Import Libraries: 
```{r import}
dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
