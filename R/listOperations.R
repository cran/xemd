listOperations <- function() {
    lapply(sort(getNamespaceExports("xemd")), function(name) {
        list(name = name,
             arguments = names(formals(getExportedValue("xemd", name))))
    })
}
