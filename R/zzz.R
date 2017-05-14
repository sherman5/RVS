THISPKG <- "RVsharing"
.onAttach <- function(libname, pkgname)
{
    version <- packageDescription("RVsharing", fields="Version")
	packageStartupMessage(paste("Welcome to RVsharing version ", version,
        "\n", sep = "" ))
}

