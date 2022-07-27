\name{edgeRUsersGuide}
\alias{edgeRUsersGuide}

\title{View edgeR User's Guide}
\description{Finds the location of the edgeR User's Guide and optionally opens it.}

\usage{
edgeRUsersGuide(view=TRUE)
}

\arguments{
  \item{view}{logical, should the document be opened using the default PDF document reader?}
}

\value{Character string giving the file location.
If \code{view=TRUE}, the PDF document reader is started and the User's Guide is opened, as a side effect.}

\details{
The function \code{vignette("edgeR")} will find the short edgeR Vignette which describes how to obtain the edgeR User's Guide.
The User's Guide is not itself a true vignette because it is not automatically generated using \code{\link{Sweave}} during the package build process.
This means that it cannot be found using \code{vignette}, hence the need for this special function.

If the operating system is other than Windows, then the PDF viewer used is that given by \code{Sys.getenv("R_PDFVIEWER")}.
The PDF viewer can be changed using \code{Sys.putenv(R_PDFVIEWER=)}.
}

\seealso{
\code{\link{system}}
}

\author{Gordon Smyth}

\examples{
# To get the location:
edgeRUsersGuide(view=FALSE)
# To open in pdf viewer:
\dontrun{edgeRUsersGuide()}
}

\concept{Documentation}
