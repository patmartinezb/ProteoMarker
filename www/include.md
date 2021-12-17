
<style>

table{width:100%;}

</style>

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

Example:

|         |             |                    |              |             |               |                  |                                  |         |                |          |              |
|:-------:|:-----------:|:------------------:|:------------:|:-----------:|:-------------:|:----------------:|:--------------------------------:|:-------:|:--------------:|:--------:|:------------:|
| **Run** | **Mixture** | **TechRepMixture** | **Fraction** | **Channel** | **Condition** | **BioReplicate** |             **key**              | **Age** | **Inhibitors** | **Type** | **Bleeding** |
|  run1   |  Mixture1   |         1          |      1       |  channel.1  |     Norm      |       Norm       | Reporter.intensity.corrected.1.1 |   NA    |       NA       |    NA    |      NA      |
|  run1   |  Mixture1   |         1          |      1       |  channel.2  |    Control    |    Control_1     | Reporter.intensity.corrected.2.1 |   23    |    Control     | Control  |   Control    |
|  run1   |  Mixture1   |         1          |      1       |  channel.3  |    Control    |    Control_2     | Reporter.intensity.corrected.3.1 |   23    |    Control     | Control  |   Control    |
|  run1   |  Mixture1   |         1          |      1       |  channel.4  |    Control    |    Control_3     | Reporter.intensity.corrected.4.1 |   32    |    Control     | Control  |   Control    |
|  run1   |  Mixture1   |         1          |      1       |  channel.5  |     Hema      |      Hema_1      | Reporter.intensity.corrected.5.1 |   28    |       no       |   HAG    |     yes      |
|  run1   |  Mixture1   |         1          |      1       |  channel.6  |     Hema      |      Hema_2      | Reporter.intensity.corrected.6.1 |    8    |       no       |   HAG    |     yes      |

------------------------------------------------------------------------

## **Requirements for the data matrix file**

Regardless of the source of the file (MaxQuant, Proteome Discoverer, etc.), the data matrix should only include a column for the **protein identifications** ("Protein.IDs" for MQ, and "Accession" for the rest), and the columns with the **actual data** (with their names matching the key column from the Metadata file).

------------------------------------------------------------------------

## **Imputation methods**

All imputation methods were taken from the {imputeLCMD} and {impute} package. From their documentation:

-   **MinDet** - Imputation of left-censored missing data using a deterministic minimal value approach. Performs the imputation of left-censored missing data using a deterministic minimal value approach.Considering a peptide/protein expression data matrix with n columns corresponding to biological samples and p lines corresponding to peptides/proteins, for each sample (column), the missing entries are replaced with a minimal value observed in that sample. The minimal value observed is estimated as being the q-th quantile (e.g. q = 0.01) of the observed values in that sample.

-   **MinProb** - Imputation of left-censored missing data using stochastic minimal value approach. Performs the imputation of left-censored missing data by random draws from a Gaussian distribution centered in a minimal value. Considering a peptide/protein expression data matrix with n columns corresponding to biological samples and p lines corresponding to peptides/proteins, for each sample (column), the mean value of the Gaussian distribution is set to a minimal value observed in that sample. The minimal value observed is estimated as being the q-th quantile (e.g. q = 0.01) of the observed values in that sample. The standard deviation is estimated as the median of the peptide/protein-wise standard deviations. Note that when estimating the standard deviation of the Gaussian distribution, only the peptides/proteins which present more than 50% recorded values are considered.

-   **KNN** - A function to impute missing expression data, using nearest neighbor averaging.

-   **Min** - Manual function to impute missing expression data, using the minimum value of the data matrix.

------------------------------------------------------------------------

## **Normalization methods**

For the mean (meanNorm) and median (medianNorm) normalization methods, functions were created manually; and for the quantile normalization method (quantNorm) the normalize.quantiles function, from the {preprocessCore} package, was used.

------------------------------------------------------------------------

## **Batch effect**

Although the option to remove the batch effect (if there is any) is given (via de removeBatchEffect function, from the {limma} package), this is not recommended, but rather it is better to model said batch effect (including covariates in the model formula). This can be done with {limma, but to address batch effect and then use the Student's T test, one must use the bacth-corrected dataset.

Thus, {limma} can work with a batch-corrected dataset or modeling the batch (**not both**, so be careful when setting the parameters of the limma DA); while the Student's T test needs the batch-corrected dataset (if there is batch effect, if not, use the original dataset).
