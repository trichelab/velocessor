% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/izar_transform.R
\name{izar_transform}
\alias{izar_transform}
\title{generalizes the transformation from Izar et al 2020 for plate-seq data}
\usage{
izar_transform(txis, orig = "tpm", dupes = 10, pseudo = 1)
}
\arguments{
\item{txis}{a SingleCellExperiment}

\item{orig}{the name of the assay to transform ("tpm")}

\item{dupes}{an estimated duplication rate (default is 10, per Izar)}

\item{pseudo}{a pseudocount to add to (TPM/dupes) (default 1, ibid)}
}
\value{
\preformatted{       a SingleCellExperiment with assay `izar`
}
}
\description{
In order to compare the results of plate-seq and droplet-seq preparations,
Izar and coauthors proposed a transformation of transcripts per million
(TPM) for transcript \code{i} in cell \code{j} from a non-UMI plate-seq library as:
}
\details{
E\link{i,j} = log2( (TPM\link{i,j} / 10) + 1 )

Here we generalize this slightly by allowing a variable scaling and offset.
The goal is to make non-UMI plate-seq and UMI droplet-seq data comparable.
We caution the user that 'comparable' is a subjective term here.

Most plate-seq protocols do not add unique molecular indices (UMIs) to
each template fragment in a library (plate-seq protocols WITH UMIs include
Quartz-Seq2, SMART-Seq2, and STORM-UMI). In contrast, almost all droplet
protocols apply both a cell barcode (CB) and UMI barcode (UB) to each
template fragment, disambiguating whether two fragments that map to the
same reference sequence are from the same template molecule or not. In
order to compare the results of plate-seq and droplet-seq preparations,
Izar and colleagues proposed this transformation. One may additionally
apply methods such as \code{sctransform} to the resulting estimate, harmonize
the resulting matrix, or similar shenanigans. Alternatively, one may use
a UMI-enabled plate-seq preparation to elide this transformation.
}
\references{
Izar B, Tirosh I, Stover EH, et al. A single-cell landscape of high-grade
serous ovarian cancer. Nat Med 26, 1271–1279 (2020).
https://doi.org/10.1038/s41591-020-0926-0
}
