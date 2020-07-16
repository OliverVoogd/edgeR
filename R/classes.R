require(methods)

# S4 classes

setClass("DGEExact",
representation("list")
)

setClass("DGEList",
representation("list")
)

setClass("DGEGLM",
representation("list")
)

setClass("DGELRT",
representation("list")
)

setClass("TopTags",
representation("list")
)

# Set inheritance
# The LargeDataObject class is set in limma and provides a show method

setIs("DGEList","LargeDataObject")
setIs("DGEExact","LargeDataObject")
setIs("DGEGLM","LargeDataObject")
setIs("DGELRT","LargeDataObject")

# Show method

setMethod("show", "TopTags", function(object) {
	if(object$test=="exact") {
		cat("Comparison of groups: ",paste(rev(object$comparison),collapse="-"),"\n")
	} else {
		cat("Coefficient: ",object$comparison,"\n")
	}
	print(object$table)
})
