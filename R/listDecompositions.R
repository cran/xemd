# Methods that decomposeTSL accepts: the scalar decompositions, applied to
# each OCID of a group, and the coupled ones, applied to the OCID of a group
# together. transformSST, the analysis operations, resiftModes and
# estimateMaskFrequency have no record operation.
listDecompositions <- function() {
    c("decomposeEMD", "decomposeEEMD", "decomposeCEEMD", "decomposeCEEMDAN", "decomposeICEEMDAN",
      "decomposeMaskedEMD", "decomposeIteratedMaskedEMD", "decomposeVMD", "decomposeSVMD",
      "decomposeLMD", "decomposeIF", "decomposeEWT", "decomposeTVFEMD",
      "decomposeBEMD", "decomposeMEMD", "decomposeNAMEMD", "decomposeMVMD")
}
