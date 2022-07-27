#  ZZZ.R

.onAttach <- function(libname, pkgname)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	14 Sep 2009. Last modified 9 July 2011.
{
	if(.Platform$OS.type=="windows" && .Platform$GUI=="Rgui" ) {
		winMenuAddItem("Vignettes","edgeR","shell.exec(system.file(\"doc\",\"edgeRUsersGuide.pdf\",package=\"edgeR\"))")
	}
}
