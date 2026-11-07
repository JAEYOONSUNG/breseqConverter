\name{export_changed_sequences_to_genbank}
\alias{export_changed_sequences_to_genbank}
\title{Export GenBank Files Containing Reconstructed Sample Sequences}
\usage{
export_changed_sequences_to_genbank(
  reference_file,
  mutation_list,
  output_dir,
  samples = names(mutation_list),
  sample_labels = NULL,
  replicon_losses = NULL
)
}
\arguments{
\item{reference_file}{GenBank or FASTA reference, including multiple replicon records.}
\item{mutation_list}{Named per-sample tables of final mutation calls, with position,
mutation, and a replicon ID column. Display-only coverage rows are excluded.}
\item{output_dir}{Directory for per-sample, per-replicon GenBank and FASTA files.}
\item{samples}{Sample IDs in output order. A missing \code{WT} entry is an empty reference baseline.}
\item{sample_labels}{Optional named display labels for sample IDs.}
\item{replicon_losses}{Optional validated whole-replicon coverage-loss table with
sample, seq_id, evidence_origin, and depth_ratio columns. Evidence must be
\code{COVERAGE_RATIO} with a nonnegative ratio below 0.1.}
}
\value{
A list containing \code{manifest}, \code{feature_audit}, \code{edits}, and
\code{call_report} tables. The manifest includes absolute output paths, sequence
lengths, sequence MD5 checksums, and explicit loss or unresolved-call statuses.
}
\description{
Reconstructs nucleotide sequences from final mutation calls rather than adding
annotations to an unchanged reference. Reference features are lifted through
edits; fully deleted features are removed. Safe altered CDS translations using
translation tables 1 and 11 are recomputed, while invalid, partial, or exceptional
altered CDS features are conservatively flagged for reannotation without stale
translations or reference protein identifiers.
Applied edits receive colored \code{misc_feature} markers in the reconstructed
coordinates, with their original reference coordinates retained in notes.
Deletion markers identify an adjacent remaining base at the deletion junction;
they do not represent a positive-length deleted span in the new sequence.
}
\details{
Uncertain missing-coverage intervals are never applied as exact sequence edits.
Whole-replicon losses have manifest entries without empty sequence files.
Unsupported, overlapping, or reference-mismatching calls yield
\code{unresolved_calls} with no sequence file for that sample and replicon.
Supported feature locations include ranges, fuzzy endpoints, complement,
join, order, and adjacent between-base sites. Unsupported locations are omitted
and recorded in the feature audit. Output filenames end in
\code{_called_variants.gbk} to distinguish a reconstruction of finalized calls
from a complete independently assembled genome.
}
