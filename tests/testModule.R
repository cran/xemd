library(xemd)

# Contract of the module description and the operation catalogue.
Module <- describeModule()
stopifnot(identical(names(Module), c("package", "version")), identical(Module$package, "xemd"),
          identical(Module$version, as.character(utils::packageVersion("xemd"))))
Operations <- listOperations()
Exports <- sort(getNamespaceExports("xemd"))
stopifnot(length(Operations) == length(Exports),
          identical(vapply(Operations, `[[`, character(1L), "name"), Exports),
          all(vapply(Operations, function(x) identical(names(x), c("name", "arguments")), logical(1L))),
          all(vapply(seq_along(Operations), function(i) {
              identical(Operations[[i]]$arguments, names(formals(get(Exports[i], envir = asNamespace("xemd")))))
          }, logical(1L))))
