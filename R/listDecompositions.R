# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
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
