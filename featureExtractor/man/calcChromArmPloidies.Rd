% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcChromArmPloidies.main.R
\name{calcChromArmPloidies}
\alias{calcChromArmPloidies}
\title{Calculate overall chrom arm copy numbers}
\usage{
calcChromArmPloidies(
  cnv.file = NULL,
  cnv = NULL,
  out.file = NULL,
  sel.cols = NULL,
  mode = c("total_cn", "minor_cn"),
  min.rel.cum.segment.size = if (mode[1L] == "minor_cn") {     0.9 } else {     0.5 },
  max.rel.cum.segment.size.diff = 0.1,
  keep.chroms = paste0("chr", c(1:22, "X")),
  one.armed.chroms = c(13, 14, 15, 21, 22),
  verbose = F
)
}
\arguments{
\item{cnv.file}{Path to purple cnv file}

\item{cnv}{purple cnv file loaded as a dataframe}

\item{out.file}{Path to output file. If NULL, returns a named vector}

\item{sel.cols}{A character vector with the names: chrom, start, end, total_cn, major_cn,
minor_cn. The value corresponding to each name should refer to a column name in the txt file.
This is used to translate the column names in the txt file to the column names that the function
will use.}

\item{mode}{Can be total_cn (for determining arm CN changes) or minor_cn (for determining arm LOH)}

\item{min.rel.cum.segment.size}{If a chrom arm has a CN category that covers >0.5 (i.e 50\%;
default) of a chrom arm, this CN is the copy number of the arm}

\item{max.rel.cum.segment.size.diff}{This value (default 0.1) determines whether which CN
categories are considered to cover equal lengths of the chrom arm. For example, (by default) 2
CN categories covering 0.40 and 0.31 of a chrom arm are considered equally contributing. When
these CN categories have similar cumulative segment size as the one with the highest, if one of
these have the same CN as the genome CN, return the genome CN. Otherwise, simply return the one
with the highest segment support (as is done above).}

\item{keep.chroms}{Which chromosomes to keep?}

\item{one.armed.chroms}{Which chromosomes are considered to have only one arm?}

\item{verbose}{Show progress messages?}
}
\value{
A named vector of chrom arm copy numbers, or specified writes a table to out.file if
specified
}
\description{
This function first rounds copy numbers (CN) to integers so that CN segments can be
grouped together. Per chrom arm, the coverage of each CN category is calculated (i.e.
cumulative segment size). The chrom arm CN is (roughly) defined as the CN category with the
highest cumulative segment size
}
