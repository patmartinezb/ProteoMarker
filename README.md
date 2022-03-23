## ProteoMarker: an easy-to-use tool to process, analyze and detect biomarkers from TMT proteomics data

ProteoMarker is a Shiny-based web application that enables the analysis of TMT proteomics data matrix (from MaxQuant, Proteome Discoverer or customized by the user), with the end goal of creating or selecting a representativa biomarker panel, capable of distinguising among conditions. This tool provides an interactive interface, easy to use for those users with no experience in the R programming language.

```{=html}
<img src="workflow.png" height="200" width="900"/>
```
## Sidebar tabs

1.  **Upload data**: Upload both the metadata and data matrix, and log2-transform the data (optional)

2.  **Quality control (QC) and filtering**: Includes sub-tabs for missing value plots, graphical overview of the data and filtering observations (proteins) based on the percentage of missing values present

3.  **Normalization**: Normalization of the data, allowing differente methods (optional)

4.  **Imputation**: Imputation of the normalizad data, allowing differente methods (optional)

5.  **Batch effect**: Includes sub-tabs for batch effect diagnosis and correction (optional)

6.  **Differential expression (DE) analysis**: Performs differential expression analysis, allowing different approaches depending on the number of groups to compare - diagnostic plots are included

7.  **Biomarker selection**: Includes sub-tabs for data partition, diagnostic plots and table for the logistic Lasso regression, and evaluation plots for the final model

8.  **Help**: Information of what variables the metadata and data matrix should include, and the rationale behind the different normalization and imputation methods

## **Requirements for the data matrix file**

Regardless of the source of the file (MaxQuant, Proteome Discoverer, etc.), the data matrix should only include a column for the **protein identifications** ("Protein.IDs" for MQ, and "Accession" for the rest), and the columns with the **actual data** (with their names matching the `key` column from the Metadata file).

## **Requirements for the metadata file**

For the app to work, the following variables have to be present:

-   **Run**: the MS run

-   **Mixture**: coincides with the MS run, and refers to the TMT batch

-   **TechRepMixture**: technical replicates (if any); in case there are none, use 1 for every sample

-   **Fraction**: fractions present (if any); in case there are none, use 1 for every sample

-   **Channel**: specific tag within the TMT batch

-   **Condition**: experimental design (i.e. control, treated)

-   **BioReplicate**: each sample within each condition

-   **key**: overlapping name of each sample with its name in the raw data matrix

Besides these mandatory columns, **other covariables can be added** (i.e. sex, age...).

> For the app to work, it is **mandatory** that the headers of the metadata file coincide with the ones stated above, and that the keys are the same in both the metadata and the data matrix.

The **reference channel** (if there is one), has to be marked as "Norm" in the Condition and BioReplicate variables of the metadata file.
