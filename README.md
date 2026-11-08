## xemd

`xemd` provides R functions for signal decomposition and mode analysis. The
library and the `xemd` command-line interface are installed together from the
repository root with the product installer.

In empirical mode decomposition (EMD), the signal's extrema define the local
scales used to extract components. Each extracted component is subtracted
from the remaining signal, leaving a final residue. Hilbert analysis describes
each component through its analytic signal: its modulus gives instantaneous
amplitude, and the derivative of its continuous phase divided by $2\pi$ gives
instantaneous frequency where that phase is defined and differentiable. The
[introduction](https://kanameishi.github.io/xemd/articles/xemd.html) develops these concepts alongside an
EMD example and a variational mode decomposition (VMD) example.

A decomposition gives a representation of the observed signal. Its modes
do not automatically identify separate physical sources or isolate noise.
The [methods guide](https://kanameishi.github.io/xemd/articles/methods.html) explains how the methods define
their components; the [interpretation guide](https://kanameishi.github.io/xemd/articles/validation.html)
distinguishes numerical convergence, reconstruction, mode diagnostics and
evidence for retaining a component.

## Guides

| Task | Documentation |
| --- | --- |
| Understand modes, residues and instantaneous frequency | [Introduction](https://kanameishi.github.io/xemd/articles/xemd.html) |
| Choose a method and understand its controls | [Methods and controls](https://kanameishi.github.io/xemd/articles/methods.html) |
| Organize seismic records and reconstruct a selection | [Seismic records](https://kanameishi.github.io/xemd/articles/seismic.html) |
| Interpret diagnostics, noise comparisons and selection | [Validation and interpretation](https://kanameishi.github.io/xemd/articles/validation.html) |
| Check exact arguments, units and returned fields | [Function reference](https://kanameishi.github.io/xemd/reference/index.html) |
| Find the published work on the methods | [Bibliography](https://kanameishi.github.io/xemd/articles/references.html) |

## EMD example

This signal combines sinusoidal oscillations with a trend. The sifting
controls specify the envelope-centering test, the boundary treatment and
the iteration limit; they are the configuration of this example.

```r
library(xemd)

Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = Inf)
OUT <- decomposeEMD(Signal, sifting = Sifting)
```

`OUT$modes` contains the extracted components and `OUT$residue` contains
what remains. Their sum reconstructs the input up to rounding. This additive
identity does not establish whether an extracted mode recovers one of the
original sinusoids. `OUT$converged` reports acceptance by the chosen sifting
rule, while `OUT$stopReason` reports why extraction ended. The
[worked introduction](https://kanameishi.github.io/xemd/articles/xemd.html#empirical-mode-decomposition)
shows those outputs, a reconstruction check and the IMF diagnostics.
