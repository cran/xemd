# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
listOperations <- function() {
    lapply(sort(getNamespaceExports("xemd")), function(name) {
        list(name = name,
             arguments = names(formals(getExportedValue("xemd", name))))
    })
}
