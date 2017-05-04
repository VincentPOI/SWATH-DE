\name{DEprot}
\alias{DEprot}
\title{
 SWATH differential expression analysis
}
\description{
DEprot is used to perform a differential expression analysis on a set of SWATH mass spectrometry protein quantification data. The functions apply different type of normalization on the dataset and calculates the foldchange and the student's test p value for the differents samples present in the data. The function produce a Volcano Plot for each comparison as well as tables of the results.
}
\usage{
DEprot(data, nbCond = NULL, nbRep = NULL, condName = NULL, normalization = "SCALING",
       pvalue.treshold = 0.05, fc.treshold = 1)
}

\arguments{
  \item{data}{
a data.frame, the row names being the proteins names and the differents columns being the proteins quantification for each samples (replicates of a same conditions must be adjoining columns).
}
  \item{nbCond}{
number of conditions in the dataset.
}
  \item{nbRep}{
number of replicates per conditions in the dataset (all conditions must have the same number of replicates).
}
  \item{condName}{
a vector of strings representing the names of the differents conditons, if NULL or the size of the vector doesn't correspond to the number of conditions the names will be numbers from 1 to the number of conditions.
}
  \item{normalization}{
the normalization technique to use. "mean" for mean centering and scaling normalization, "quantile" for quantile normalization, "median" for median centering and scaling, "null" to not perform any normalisation in case the dataset is already normalized.
}
  \item{pvalue.treshold}{
p value treshold value to select the significamtly differentially expressed proteins.
}
  \item{fc.treshold}{
fold change treshold value to select the significamtly differentially expressed proteins.
}
}
\details{
The function takes in entry protein quantification data (issued from SWATH mass spectrometry), then it performs normalization on the dataset to reduce the noise caused by 
This package allows to perform differential expression analysis on protein quantification data.
The package was made for data issued from SWATH mass spectrometry used to quantify proteins in samples 
from differents conditions 
}
\value{
DEprot returns a list containing the different components :
\item{normdata}{
a data frame of the normalized data
}\item{resultsTable}{
a matrix of dataframes containing the protein name, the p value and the fold change for each comparison

}\item{DEprotein}{
a matrix of dataframes containing the differentialy expressed protein name, their p value and their fold change for each comparison.
}

}
\references{
exemple dataset :
Bjelosevic S, Pascovici D, Ping H, et.al. Quantitative age-specific variability of plasma proteins in healthy neonates, children and adults. Molecular & Cellular Proteomics (March 23 2017)
}
\examples{

data(SWATH.example.data)
results<-DEprot(SWATH.example.data,nbCond = 4,nbRep = 10, condName = c("control","patients",""),
                normalization = "scaling", fc.threshold = 0.5, pvalue.threshold = 0.05)

results$dataprocessed[1,2] # display the data frame of the results for the comparison of the two conditions

results$DEprotein[1,2] # display the significantly expressed proteins data frame for the comparison of the two conditions
}
