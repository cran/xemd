## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----emd----------------------------------------------------------------------
library(xemd)

Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = Inf)
EMD <- decomposeEMD(Signal, sifting = Sifting)
str(EMD)

## ----reconstruction-----------------------------------------------------------
max(abs(rowSums(EMD$modes) + EMD$residue - Signal))

## ----assess-------------------------------------------------------------------
assessIMF(cbind(EMD$modes, EMD$residue))

## ----vmd----------------------------------------------------------------------
SamplingInterval <- 1 / 1000
TimeV <- seq_len(1000) * SamplingInterval
SignalV <- cos(2 * pi * 2 * TimeV) + 0.25 * cos(2 * pi * 24 * TimeV) +
    0.0625 * cos(2 * pi * 288 * TimeV)
VMD <- decomposeVMD(SignalV, samplingInterval = SamplingInterval,
                    initialCenters = (0:2) / (6 * SamplingInterval),
                    alpha = 2000 * SamplingInterval^2 / (8 * pi^2),
                    tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6,
                    maxIterations = 98)
VMD$centers
VMD$stopReason

## ----attributes---------------------------------------------------------------
Attributes <- computeInstantaneousAttributes(EMD$modes,
                                             samplingInterval = 9 / 2999)
apply(Attributes$frequency, 2, stats::median)

## ----spectrum-----------------------------------------------------------------
Spectrum <- buildHilbertSpectrum(Attributes$frequency, Attributes$amplitude,
                                 edges = c(1, 2, 3, 4, 5), exponent = 2,
                                 samplingInterval = 9 / 2999)
Spectrum$marginal

