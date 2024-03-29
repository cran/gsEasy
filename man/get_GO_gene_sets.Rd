% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Interface.R
\name{get_GO_gene_sets}
\alias{get_GO_gene_sets}
\title{Create list of gene sets defined by GO term annotation}
\usage{
get_GO_gene_sets(
  GO_annotation_file,
  GO_file = "http://purl.obolibrary.org/obo/go.obo",
  min_genes = 15,
  max_genes = 500,
  verbose = TRUE
)
}
\arguments{
\item{GO_annotation_file}{File path of annotation file, which should contain a column of genes and a column of terms. Can be downloaded from at http://geneontology.org/gene-associations/gene_association.goa_human.gz.}

\item{GO_file}{File path of gene ontology.}

\item{min_genes}{Minimum number of genes in gene sets.}

\item{max_genes}{Maximum number of genes in gene sets.}

\item{verbose}{Print progress.}
}
\value{
List of character vectors of term IDs.
}
\description{
Note, this function takes several minutes to execute.
}
