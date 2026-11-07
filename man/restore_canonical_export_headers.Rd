\name{restore_canonical_export_headers}
\alias{restore_canonical_export_headers}
\title{Restore Canonical Headers from a Presentation Export}
\usage{restore_canonical_export_headers(file, sheet, data)}
\arguments{
\item{file}{Path to the Excel workbook.}
\item{sheet}{Worksheet name or index.}
\item{data}{Data frame read without dropping empty columns.}
}
\value{The data frame with canonical headers and sample-order metadata, or unchanged if no metadata exists.}
\description{Internal loader for clean-header workbook round trips. It validates display headers before restoring their canonical names.}
\keyword{internal}
