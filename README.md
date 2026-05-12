# PySOC - Calculation of Spin-Orbit Coupling

PySOC is a Python and Fortran program for calculating spin-orbit coupling (SOC)
between singlet ground or excited states and triplet excited states from
linear-response excited-state calculations.

PySOC currently supports excited-state data from:

- Gaussian 09/16 TDDFT and TDA calculations
- DFTB+ TD-DFTB calculations
- ORCA TDDFT/TDA calculations

## Features

- Evaluation of SOC matrix elements between singlet and triplet states.
- Text-table and CSV command-line output.
- Gaussian TDDFT/TDA interface using `.log`/`.rwf` data.
- DFTB+ TD-DFTB interface using transition-vector and fitted-basis data.
- ORCA TDDFT/TDA interface using printed orbital data plus full TD vectors from
  `orca_2json`.
- One-electron, full Breit-Pauli, and effective-charge (`Zeff`) SOC modes.
- Cartesian integral backend through MolSOC and `soc_td`; ORCA spherical
  harmonic AO data are reordered and converted internally for the supported
  shells.

## Dependencies

Install the Python dependencies:

```console
$ pip install -r requirements.txt
```

or equivalently:

```console
$ pip install cclib tabulate periodictable scipy
```

PySOC also needs its Fortran backend executables:

- `pysoc/data/bin/molsoc`
- `pysoc/data/bin/soc_td`

The repository includes these paths. If you recompile either Fortran backend,
copy the newly built executable back into `pysoc/data/bin` before testing or
running `pysoc`.

Program-specific external tools:

- Gaussian calculations require Gaussian's `rwfdump` executable in `PATH`.
- ORCA calculations should have `orca_2json` in `PATH` for complete TD vectors.
- ORCA calculations also need an `.mkl` file, or `orca_2mkl` in `PATH` so PySOC
  can generate it from the `.gbw` file.

## Installation

Add the repository to `PYTHONPATH`:

```console
$ export PYTHONPATH="$PYTHONPATH:/path/to/pysoc"
```

Add the `bin` directory to `PATH`:

```console
$ export PATH="$PATH:/path/to/pysoc/bin"
```

The `bin/pysoc` command is a symlink to `pysoc/program/main.py`.

## General Requirements

Before running PySOC, first run the electronic-structure excited-state
calculation in Gaussian, DFTB+, or ORCA.

> [!NOTE]
> Shells higher than the f shell are not currently supported by the full PySOC
> workflow. Use basis sets without g or higher angular-momentum shells.

The command-line interface can guess the program from the main file extension:

- `.log` -> Gaussian
- `.xyz` -> DFTB+
- `.out` -> ORCA

If a file extension is ambiguous, pass `--program` explicitly. This is important
for Gaussian jobs whose main output file is named `.out`, because `.out` is
guessed as ORCA.

## Gaussian Calculations

Gaussian TDDFT and TDA calculations are supported for Gaussian 09 and 16. Other
Gaussian methods may work if the required `.rwf` sections are present, but they
are not the tested path.

PySOC requires:

- The Gaussian output file, normally `.log`.
- The Gaussian read-write file, `.rwf`.
- Gaussian's `rwfdump` executable in `PATH`.
- Cartesian d and f functions in Gaussian output ordering.
- Printed basis-set information in the output file.

Gaussian normally deletes the `.rwf` file at the end of a job. Keep it by adding
an `%rwf` line to the Gaussian input:

```text
%rwf=formaldehyde.rwf
```

Use `6D 10F GFInput` so PySOC receives Cartesian d/f functions and the basis set
is printed:

```text
%mem=1GB
%chk=formaldehyde.chk
%rwf=formaldehyde.rwf
# TD(50-50,nstates=5) wB97XD/TZVP 6D 10F GFInput

formaldehyde

0 1
C         -0.131829      -0.000001      -0.000286
O          1.065288       0.000001       0.000090
H         -0.718439       0.939705       0.000097
H         -0.718441      -0.939705       0.000136
```

Run PySOC with:

```console
$ pysoc formaldehyde.log --rwf_file formaldehyde.rwf
```

If `--rwf_file` is omitted, PySOC looks for a file with the same basename as the
Gaussian output and the `.rwf` suffix.

## ORCA Calculations

ORCA support is designed for TDDFT/TDA calculations with singlet and triplet
roots. ORCA uses spherical harmonic AO functions, while MolSOC and `soc_td`
operate internally with Cartesian Gaussian functions. PySOC handles the
spherical ordering and Cartesian conversion for the supported shells, using the
GTO basis data printed by ORCA.

PySOC requires:

- The ORCA `.out` file.
- The ORCA `.gbw` file.
- The ORCA `.cis` file for TD-vector export through `orca_2json`.
- The ORCA `.property.txt` file when produced by ORCA.
- The ORCA `.mkl` file, or `orca_2mkl` in `PATH` so PySOC can create it.
- `orca_2json` in `PATH` for complete TD vectors.

Recommended ORCA input options:

```text
! wB97X-D3 TZVP

%tddft
    nroots 10
    triplets true
end

%output
    Print[ P_Basis ] 2
    Print[ P_Overlap ] 1
    Print[ P_MOs ] 1
end

* xyz 0 1
C         -0.131829      -0.000001      -0.000286
O          1.065288       0.000001       0.000090
H         -0.718439       0.939705       0.000097
H         -0.718441      -0.939705       0.000136
*
```

`Print[ P_Basis ] 2` is strongly recommended because it lets PySOC use ORCA's
raw GTO basis block. If that block is missing, PySOC falls back to basis data in
the `.mkl` file and prints a warning because the contraction normalization may
not be the same as MolSOC expects.

Run PySOC with:

```console
$ pysoc ch2o_orca.out --program ORCA -s 5 -t 5
```

When ORCA support is working in the full-vector path, PySOC automatically writes
or updates `ch2o_orca.json.conf`, runs:

```console
$ orca_2json ch2o_orca.gbw -json
```

and reads the complete `X` and `Y` TD vectors from `ch2o_orca.json`. For TDA
roots, ORCA does not provide `Y`; PySOC correctly uses `X` for both `X+Y` and
`X-Y`.

If `orca_2json` or the JSON TD-vector data is unavailable, PySOC falls back to
the sparse CI amplitudes printed in the `.out` file and emits a warning. That
fallback is useful for debugging, but it is not recommended for quantitative SOC
values because ORCA usually prints only amplitudes above a threshold.

Request more ORCA roots than you plan to use in PySOC. Root ordering near the
highest requested Davidson root can be unstable; PySOC warns when selected roots
include the highest root printed by ORCA.

## DFTB+ Calculations

> [!NOTE]
> This section is under revision.

DFTB+ TD-DFTB calculations are supported through the original PySOC/TD-DFTB
interface.

> [!NOTE]
> PySOC's SOC backend uses Gaussian-type orbital data. DFTB+ uses Slater-type
> orbitals, so the parameter set used for the DFTB+ calculation must be fitted
> to GTOs for use with PySOC. PySOC includes a fitted set for `mio-1-1`, and an
> alternative fitted basis directory can be supplied with `--fitted_basis`.

An example `dftb_in.hsd`:

```text
Geometry = GenFormat {
  <<< "ch2o.gen"
}

Driver = {}

Hamiltonian = DFTB {
  SCC = Yes
  SCCTolerance = 1e-10
  MaxAngularMomentum = {
    H  = "s"
    C  = "p"
    O  = "p"
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/path/to/mio-1-1/"
    Separator = "-"
    Suffix = ".skf"
  }
  LinearResponse {
    NrOfExcitations = 10
    StateOfInterest = 0
    Symmetry = both
    HubbardDerivatives {
      H  = 0.347100   0.491900
      C  = 0.341975   0.387425
      O  = 0.467490   0.523300
    }
    WriteTransitions = Yes
    WriteTransitionDipole = Yes
    WriteXplusY = Yes
  }
}

Options {
  WriteEigenvectors = Yes
  WriteHS = No
}

ParserOptions {
  ParserVersion = 4
}
```

After the first DFTB+ calculation, set `WriteHS = Yes` and run the calculation a
second time so the Hamiltonian and overlap data are written.

Run PySOC with the geometry file:

```console
$ pysoc formaldehyde.xyz --program "DFTB+"
```

## Usage

Basic command:

```console
$ pysoc CALC_FILE
```

Examples:

```console
$ pysoc examples/ch2o_gaussian/gaussian.log -s 5 -t 5
$ pysoc examples/ch2o_orca/ch2o_orca.out --program ORCA -s 5 -t 5
$ pysoc examples/ch2o_tddftb/ch2o.xyz --program "DFTB+" -s 5 -t 5
```

Common options:

| Option | Meaning |
| --- | --- |
| `-p`, `--program` | Input program: `Gaussian`, `DFTB+`, or `ORCA`. |
| `-s`, `--singlets` | Number of singlet excited states to include. Defaults to all parsed states. |
| `-t`, `--triplets` | Number of triplet excited states to include. Defaults to all parsed states. |
| `-T`, `--calculation` | SOC integral mode: `auto`, `one`, `two`, or `zeff`. |
| `-S`, `--SOC_scale` | Scale factor for `Zeff` mode. |
| `-n`, `--no_ground` | Exclude the singlet ground state from the SOC table. |
| `-C`, `--CI_threshold` | Threshold for CI coefficients passed to `soc_td`. |
| `-c`, `--CSV` | Print comma-separated output instead of a formatted table. |
| `-o`, `--output` | Directory where intermediate files should be kept. |
| `--rwf_file` | Gaussian `.rwf` file. |
| `--fitted_basis` | DFTB+ fitted GTO basis directory. |

At present, `-s` and `-t` select the first N singlet and triplet states. They do
not select arbitrary state lists.

`--calculation auto` is the default. It uses `zeff` when MolSOC has effective
charges for all atoms in the molecule; otherwise it falls back to one-electron
SOC and prints a warning.

## Output

The default output is a table:

```console
$ pysoc examples/ch2o_gaussian/gaussian.log -s 5 -t 5
 Singlet    Triplet     RSS (cm-1)    +1 (cm-1)    0 (cm-1)    -1 (cm-1)
---------  ---------  ------------  -----------  ----------  -----------
  S(0)       T(1)          60.7419      42.9510      0.0077      42.9510
  S(0)       T(2)           0.0194       0.0137      0.0000       0.0137
  S(0)       T(3)          10.6234       0.0133     10.6234       0.0133
  S(0)       T(4)          59.8873      42.3467      0.0012      42.3467
  S(0)       T(5)          11.7332       0.0013     11.7332       0.0013
```

Each line is one singlet-triplet pair. `+1`, `0`, and `-1` are the triplet
sub-state couplings. `RSS` is the root-sum-square of those three components.

CSV output:

```console
$ pysoc examples/ch2o_gaussian/gaussian.log -s 5 -t 5 -c
Singlet,Triplet,RSS (cm-1),+1 (cm-1),0 (cm-1),-1 (cm-1)
S(0),T(1),60.7419154885397,42.95102,0.00769,42.95102
S(0),T(2),0.019431296920174937,0.01374,1e-05,0.01374
```

Use shell redirection to write CSV output to a file:

```console
$ pysoc examples/ch2o_gaussian/gaussian.log -s 5 -t 5 -c > SOC.csv
```

## Intermediate Files

By default, PySOC writes intermediate files to a temporary directory and removes
them at the end of the run. Use `-o` to keep them:

```console
$ pysoc examples/ch2o_orca/ch2o_orca.out --program ORCA -s 5 -t 5 -o ./pysoc_debug
```

Useful intermediate files include:

- `molsoc.inp`
- `molsoc_basis`
- `molsoc_overlap.dat`
- `molsoc_dipole.dat`
- `soint`
- `mo_ene.dat`
- `mo_coeff.dat`
- `ao_overlap.dat`
- `ci_coeff.dat`
- `soc_td_input.dat`
- `soc_out.dat`

For ORCA, the `orca_2json` cache files are written next to the ORCA output
files, not into the PySOC intermediate directory:

- `basename.json.conf`
- `basename.json`
- `basename.JSON.bibtex`

## Known Limitations

- Shells above f are not supported by the full workflow.
- ORCA quantitative SOC calculations need the `orca_2json` full-vector path.
  Sparse CI amplitudes printed in `.out` are not complete enough for production
  SOC values.
- ORCA support has been developed and tested primarily for TDDFT/TDA closed-shell
  examples.
- DFTB+ support requires a fitted GTO basis compatible with the DFTB+ parameter
  set.
- The command line currently selects the first N singlet/triplet states, not an
  arbitrary set of roots.

## Additional Documentation

1. Short tutorial: [doc/pysoc.pdf](doc/pysoc.pdf)
2. Chinese tutorial by sobereva Lu Tian: http://bbs.keinsci.com/thread-9442-1-1.html
3. Gaussian 16 tutorial by ggdh: http://bbs.keinsci.com/thread-19813-1-1.html

## Reference and Citation

Evaluation of Spin-Orbit Couplings with Linear-Response Time-Dependent Density
Functional Methods

Xing Gao, Shuming Bai, Daniele Fazzi, Thomas Niehaus, Mario Barbatti, and Walter
Thiel

J. Chem. Theory Comput., 2017, 13 (2), pp 515-524

DOI: 10.1021/acs.jctc.6b00915

If you use ORCA-generated data, also cite ORCA and any ORCA methods used in the
underlying electronic-structure calculation.

## Authorship

PySOC was originally written by Xing Gao et al. for Python 2.x.

[MolSOC](molsoc/) is used for the calculation of atomic integrals and was
originally written by Sandro Giuseppe Chiodo et al. (Computer Physics
Communications 185 (2014) 676-683).

PySOC was rewritten for Python 3.x by Oliver S. Lee.

Support for ORCA was added by Diego J. Alonso de Armiño.

