# breseqConverter 1.1.0

* Reconstruct per-sample changed sequences from finalized calls and export
  annotated GenBank plus native SnapGene DNA through the one-file workflow.
  Audit unsupported calls, feature changes, whole-replicon losses and unresolved
  missing coverage explicitly.
* Publish a generic single-file entry point with no machine-specific data paths.
* Make `one_shot.R` run automatically in a fresh R session with bundled
  converter functions, reference preparation, per-replicon plots and Excel
  export. Preserve the previous script as `one_shot_legacy.R`.
* Add observed per-sample rings and explicit outside-to-inside sample order
  for the single-file workflow, preserving audited coverage display evidence.
* Use Cairo PDF output when available for mutation glyphs, and draw GC bars
  with their actual window widths, including short final windows.
* Restore the legacy GenomeDrawer GC palette: skew green/red and ratio
  blue/yellow. Add `gc_colors` overrides and derive colors from each track's
  plotted deviation, including when input color columns are stale.
* Include the existing standalone workflow with an explicit intergenic label
  toggle and GC color settings. Preserve mutation marks when hiding labels.
* Use subtle stage colors for every observed call in opt-in Excel exports,
  including retained calls. Use bold for new or changed calls, regular weight
  for revised MC evidence, and Arial 10 pt throughout.
* Add an optional intergenic label filter to the Circos workflows. The filter
  preserves mutation marks and data and is disabled by default.
* Add optional stage-specific Excel highlighting, readable headers, display
  labels, and explicit sample order. Formatted exports retain canonical headers
  for subsequent loading by the package.
* Preserve distinct replicons and full sample IDs when comparing and exporting
  stage-specific calls.
* Fix the classification helper name collision with the GenBank exporter,
  comma-formatted deletion lengths, and classification of calls first reported
  in later sample columns. Keep GenBank classification and colors unchanged.
* Parse numeric coordinates from oriented junction endpoints without changing
  the reported endpoint strings used to identify junction events.
