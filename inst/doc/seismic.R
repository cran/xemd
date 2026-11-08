## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----table, message = FALSE---------------------------------------------------
library(xemd)
library(data.table)

Time <- seq(0, 9, length.out = 3000)
Channels <- cbind(sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
                      0.5 * Time,
                  sin(pi * Time), sin(2 * pi * Time))
IDs <- c("AT", "VT", "DT", "AI", "CAV")
OCIDs <- c("H1", "H2", "UP")
TSL <- rbindlist(lapply(seq_along(IDs), function(i) {
    rbindlist(lapply(seq_along(OCIDs), function(j) {
        data.table(RecordID = "R1", DIR = OCIDs[j], ID = IDs[i],
                   OCID = OCIDs[j], t = Time,
                   s = Channels[, (i + j - 2) %% 3 + 1])
    }))
}))
TSL[, .N, by = .(RecordID, ID, OCID)]

## ----methods------------------------------------------------------------------
listDecompositions()

## ----scalar-------------------------------------------------------------------
Input <- copy(TSL)
OUT <- decomposeTSL(TSL, method = "decomposeEMD",
                    control = list(sifting = list(stopRule = "doubleThreshold")),
                    by = "RecordID")
names(OUT)
OUT$termination
OUT$diagnostics[ID == "AT"]
identical(Input, TSL)

## ----modes--------------------------------------------------------------------
OUT$modes[, .N, by = .(OCID, ID, mode)][order(OCID, ID, mode)][1:6]

## ----coupled------------------------------------------------------------------
TimeV <- seq_len(1000) / 1000
SeriesV <- cbind(cos(2 * pi * 2 * TimeV) + 0.25 * cos(2 * pi * 24 * TimeV) +
                     0.0625 * cos(2 * pi * 288 * TimeV),
                 cos(2 * pi * 2 * TimeV), 0.25 * cos(2 * pi * 24 * TimeV))
Record <- data.table(ID = "AT", OCID = rep(OCIDs, each = 1000),
                     t = rep(TimeV, 3), s = as.vector(SeriesV))
Control <- list(initialCenters = (0:2) / (6 / 1000),
                alpha = 2000 * (1 / 1000)^2 / (8 * pi^2), tau = 0.5,
                fixedZeroCenter = FALSE, tolerance = 1e-6,
                maxIterations = 98)
Coupled <- decomposeTSL(Record, method = "decomposeMVMD", control = Control,
                        by = character())
names(Coupled)
Coupled$centers

## ----regularization-----------------------------------------------------------
Displaced <- Time
Displaced[1501] <- Time[1501] + 9 / 2999 / 2
Irregular <- rbindlist(lapply(c("AT", "VT"), function(id) {
    data.table(ID = id, OCID = rep(c("H1", "H2"), each = 3000),
               t = rep(Displaced, 2),
               s = c(sin(pi * Displaced) + sin(2 * pi * Displaced) +
                         sin(6 * pi * Displaced) + 0.5 * Displaced,
                     sin(pi * Displaced)))
}))
Regularized <- decomposeTSL(Irregular, method = "decomposeEMD",
                            control = list(sifting = list(stopRule = "doubleThreshold")),
                            by = character())
Regularized$signal[, .(samples = .N, first = min(t), last = max(t)),
                   by = .(OCID, ID)]

