## xemd 0.1.6

- Licence changed from GPL-3 to AGPL (>= 3). The new `LICENSE.note` states
  the copyright holder's reservation of text and data mining, which is not a
  licence term; the code adapted from third parties keeps the licences listed
  in `NOTICE`.
- The conditions that name a series carry it in their `unit` field instead of
  a prefix of the message: the error of `decomposeTSL()` when the channels of
  a group do not share one time grid or a series has fewer than two samples,
  now of class `xemdTimeGrid`, its warning about a non-uniform time grid, and
  the `xemdZeroEnergy` error of `measureModesTSL()`.
- In the command-line interface, an error inside a step of `run` or `batch`,
  a method of `compare` or the reading of an input keeps its class and exit
  status and locates itself in the fields `phase`, `unit` and `path`, which
  `--verbose` shows.
- `xemd doctor` diagnoses without repairing, as the family rules: one JSON
  report on standard output, exit 0 only when every required check passes,
  and `doctor --help` listing the checks. The `--json` and `--install` flags
  are gone; `xemd repair [--source DIR]` installs the unavailable R packages
  or reinstalls the library from a checkout. `xemd --version` prints the
  package with its version and then the installed build commit.
- The maintainer scripts of `inst/` (`build.R`, `package.R`, `install.R`,
  `rhub.R`, `cran-check.R`, `cran-submit.R`, `rhub-check.R`) no longer ship
  in the package tarball; they stay in the repository for the maintainer.

## xemd 0.1.5

- The conditions the library signals carry the classes of the family console:
  `xemdZeroEnergy`, `xemdInsufficientModes` and `xemdZeroEnvelope` inherit
  from `moduleFailure`, and the warning of `decomposeTSL()` about a
  non-uniform time grid is a `moduleWarning`. Callers that caught the own
  classes keep working.
- The command-line interface reads its arguments from the launcher and leaves
  the rendering of every message, warning and error to the installer's console
  renderer: an error is written as `xemd: <message>`, exit status 2 for a
  usage error (class `xemdUsage`) and 1 otherwise; the launcher accepts
  `--verbose` and `--quiet`. The package code is otherwise unchanged.

## xemd 0.1.4

- The maintainer scripts of `inst/` (`cran-check.R`, `cran-release.R`,
  `cran-submit.R`, `rhub-check.R`) travel in the source package, as in the
  other modules of the family; the package code is unchanged.
- The introductory vignette is titled "Introduction", and the vignettes,
  the articles, the README and two help pages write the package name in
  code font or with `\pkg{}` instead of as a plain word.

## xemd 0.1.3

- `decomposeEEMD()`, `decomposeCEEMD()`, `decomposeMaskedEMD()` and
  `decomposeBEMD()` no longer return a mode whose sifting lost the extrema
  needed for its envelopes. Each realization of EEMD and CEEMD stops before
  that candidate, masked EMD discards the stage when one of its phases
  loses them, and BEMD discards the incomplete vector candidate. The
  content of the attempt stays in the residue, the modes extracted before
  it are unchanged, and the divisors, phases and directions are preserved.
  A candidate that reaches the iteration limit is kept and flagged
  `converged = FALSE`, as the reference implementations do.
- Those four decompositions return `termination`, a data frame that
  records the attempted mode index, the subtractions of the attempt and its
  `stopReason` (`terminalResidue`, `modeLimit` or `insufficientExtrema`,
  and `iterationLimit` for a phase of masked EMD listed with the one that
  failed), with the realization and its sign, the phase or the first
  direction without an envelope where the method has them.
  `decomposeTSL()` carries the table with the keys of the series.
- The tests gain a script that checks, on internal inputs recorded from the
  decomposition of two horizontal components of a strong-motion record,
  that an interrupted attempt leaves the residue unchanged and reports its
  cause.

## xemd 0.1.2

- The package carries eight contract test scripts on toy cases. The suite
  exercises the exported functions and checks output structures,
  reconstruction identities, preservation of the caller's random stream
  with explicit noise, and rejection of invalid inputs. The conformance
  tests shipped in version 0.1.1 are no longer part of the package.
- The examples of `decomposeCEEMD()`, `decomposeCEEMDAN()`,
  `decomposeEEMD()`, `decomposeICEEMDAN()` and `resiftModes()` use 1000
  samples instead of 3000; `decomposeMEMD()` uses 256 instead of 512.

## xemd 0.1.1

- Code adapted from pytvfemd and from curve_ext of the Synchrosqueezing
  Toolbox, which version 0.1.0 did not declare in `DESCRIPTION`, `NOTICE` or
  the file headers, is credited; executable code is unchanged. `DESCRIPTION`
  adds Stefano Bianchi and Alessandro Longo, authors of pytvfemd, and
  Jianfeng Lu, original author of curve_ext, and gives the maintainer the
  role `trl`.
- `NOTICE` adds pytvfemd 0.3.14, with its copyright and license notice and
  the seven files that adapt its code for `decomposeTVFEMD()`, and lists in
  the SynchWave entry the two files of `extractSSTModes()` that adapt
  curve_ext and curve_ext_multi. Those nine files carry a header with the
  source, its authors, its license and the date of modification.
- The six files adapted from SynchWave keep LGPL (>= 2) as the license of
  their source and state that `xemd` distributes them under GPL-3, as section 3
  of the GNU Library General Public License version 2 permits.
- The help pages of the methods whose article `inst/REFERENCES.bib` lists
  cite it, with its DOI where the entry gives one, and the pages of
  `decomposeTVFEMD()`, `transformSST()` and `extractSSTModes()` name the
  adapted code and its license.

## xemd 0.1.0

- `decomposeEWT()` takes `detection`, a list that selects the boundary
  detector of the EWT authors' code (`locmax`, `locmaxmin` or `locmaxminf`),
  the logarithm and the regularization of the spectrum, and with
  `bandCount = NULL` counts the bands by the relative threshold of Gilles
  (2013); `detection = NULL` keeps the previous partition unchanged.
  `extractSSTModes()` traces ridges of a synchrosqueezed transform by the
  penalized dynamic programme of the Synchrosqueezing Toolbox and
  reconstructs one component per ridge from a band around it.
  `assessFractionalNoise()` compares the energy of each mode with the
  fractional Gaussian noise reference and threshold tabulated by Flandrin,
  Gonçalves and Rilling (2005).
- The help page of each function that uses code adapted from a public
  library names the library, its author and its license, as `NOTICE` does.
- Decompositions are named `decompose<METHOD>`: EMD; EEMD, CEEMD, CEEMDAN
  and ICEEMDAN, with `resiftModes()`; masked EMD, with
  `estimateMaskFrequency()`, and iterated masking EMD; VMD, MVMD and SVMD;
  BEMD, MEMD and NA-MEMD; LMD and iterative filtering; EWT; TVF-EMD. The
  synchrosqueezed transform is `transformSST()`.
- Decompositions receive `signal` without a time vector and return lists
  with the samples in the first dimension: `modes` and, where the method
  defines them, `residue` or `fitResidual`, `iterations`, `converged`,
  `stopReason` and the fields of its family. Functions whose formulation
  depends on the time scale receive `samplingInterval` in seconds, and
  frequencies are in Hz. Noise-assisted methods receive their noise
  realizations and draw no random numbers. Undefined formulas stop with the
  condition classes `xemdZeroEnergy`, `xemdZeroEnvelope` and
  `xemdInsufficientModes`.
- Analysis of modes: `computeInstantaneousAttributes()`,
  `buildHilbertSpectrum()`, `assessIMF()`, `assessWhiteNoise()`,
  `measureDispersion()`, `measureOrthogonality()`, `measurePMSI()`,
  `measureQRF()`, `measureFCF()`, `measureDTW()` and
  `selectModesEnergy()`.
- `decomposeTSL()` decomposes the series of a long table of records by
  record, quantity (`ID`) and component (`OCID`), regularizes grids that
  are not uniform with a warning, and returns homogeneous tables for every
  method. It replaces `TSL2IMF()`, retired without an alias.
  `listDecompositions()` returns the methods it accepts.
- The command-line interface, the skill and the operations catalogue
  `inst/api/operations.json` follow `decomposeTSL()`: `xemd decompose` passes the JSON `control` object as the
  named list of arguments of the method, with nested objects as lists,
  arrays as vectors and numbers as doubles, and writes the list of tables to
  an RDS file; `xemd methods` prints `listDecompositions()`. Warnings, such
  as the regularization of a time grid, reach standard error as they occur.
- `xemd doctor` reports, one line per check (`--json` for programs), whether R meets the minimum of `xemd`, `xemd`
  loads, its R packages and jsonlite are available, the installation receipt
  names the version that loads, the first `xemd` on `PATH` is the launcher of
  the same prefix and the working directory is writable; a failed check sets
  exit status 1. `xemd doctor --install` installs the unavailable R packages
  from CRAN and checks again; with `--source DIR` it also reinstalls the `xemd`
  library from that checkout, without sudo, whether or not its version
  differs. The report names no library or file path.
  `describeModule()` still returns the package and its version.
- `measureModesTSL()`, `selectModesTSL()` and `reconstructTSL()` work on the
  tables of `decomposeTSL()` for every method: the metrics of each mode
  (frequency and bandwidth in Hz, energy, energy fraction and mean period)
  and of each series (number of modes and energy fractions of the remainder
  and of the cross products); the selection of modes with a declared
  criterion (explicit indices, a minimum energy fraction or a closed
  frequency band, to keep or to remove); and the sum of the retained modes,
  with or without the remainder. The command-line interface adds the verbs
  `compare`, `metrics`, `select`, `reconstruct`, `export`, `run` and
  `batch`, and every verb answers `--help` with its keys: `run` chains
  steps on one input and writes one RDS file per step and an index, `batch`
  runs the same steps on several inputs and stacks their tables, and
  `export` writes the tables as CSV files, long or with one column per
  mode. A series of `t` and `s` is read as one record, and a method that
  draws noise requires `seed`, which the output keeps. The verbs after
  `compare` work method by method and stack their tables with the stacking
  keys first: `method`, or `source` and then `method` after `batch`. A
  criterion `modes` whose entries mix `method` with none is a usage error.
  `xemd doctor --install --source DIR` reinstalls the library whether or not
  the version differs.
