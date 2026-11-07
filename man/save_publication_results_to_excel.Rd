\name{save_publication_results_to_excel}
\alias{save_publication_results_to_excel}
\title{Save Compact Publication Mutation and Evidence Tables}
\usage{
save_publication_results_to_excel(
  merged_results,
  output_file,
  sample_order = NULL,
  sample_labels = NULL,
  reference_file = NULL,
  gff_file = NULL,
  plot_evidence = NULL,
  annotation_overrides = NULL
)
}
\arguments{
\item{merged_results}{Named list of canonical merged mutation and evidence data frames.}
\item{output_file}{Path of the additional compact Excel workbook.}
\item{sample_order}{Optional exact chronological sample IDs. Empty reference baselines may be included.}
\item{sample_labels}{Optional named map from sample IDs to display labels.}
\item{reference_file}{Optional GenBank or FASTA reference for replicon labels and gene products. FASTA alone has no gene annotations.}
\item{gff_file}{Optional GFF file providing gene product annotations.}
\item{plot_evidence}{Optional audited display-evidence data frame with sample or sample_id, seq_id, start, end, size, evidence_origin, source_start, source_end, criterion, depth_ratio and raw_final_DEL_call. A separate sheet is created; these rows are never promoted to finalized calls. Optional appended and figure_label fields control displayed rows and their explicit labels.}
\item{annotation_overrides}{Optional data frame with sheet, seq_id, position, field and value columns. Only annotation columns in Mutations, Missing coverage or Plot-derived evidence may be changed. Mutation positions are numeric reference loci; missing-coverage positions use an original start-boundary string; plot-evidence positions use the displayed start coordinate. Every override must identify exactly one compact row.}
}
\value{Invisibly, normalized compact tables, sample order, presentation settings and audit mappings.}
\description{
Writes a compact workbook in addition to the canonical wide results. Mutation
rows retain the actually observed stage alleles. Pale fills mark reported calls;
bold marks new, changed or reappearing calls, without implying de novo mutation.
All populated content uses Arial 10 and the table headers begin on row seven.
}
\details{
Missing-coverage records are grouped only when their certain intervals overlap
on the same replicon and their reported gene contexts are compatible. Families
containing only evidence within a preceding finalized deletion are omitted from
this compact view; all original records remain in the companion wide workbook.
Uncertain interval boundaries and first-report statuses remain distinct from
finalized deletion calls. Reference annotations resolve gene lists where available.

Junction rows use both oriented endpoints when canonical paired keys are present.
Without those keys, every source endpoint is retained separately and the sheet
explicitly identifies unpaired endpoint mode. No pairs are inferred from adjacency.

An entirely empty WT stage is labelled Ref. as a reference-assembly baseline.
No missing mutation call is carried forward to a later stage. Annotation aliases
are explicit overrides and never modify calls, coordinates or the input data.
}
