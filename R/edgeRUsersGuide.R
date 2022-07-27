edgeRUsersGuide <- function(view=TRUE)
#	Find and optionally view edgeR User's Guide
#	Gordon Smyth
#	23 May 2012.
{
	f <- system.file("doc","edgeRUsersGuide.pdf",package="edgeR")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}
