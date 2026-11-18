# GROMACS + MDAnalysis Trajectory Workbench

A local Python/Tkinter GUI for **GROMACS trajectory processing**, **MDAnalysis structural analysis**, **residue–residue contact maps**, and **VMD visualization**. Standard analyses can run in batches and export plotted and tabular results.

> **Research software status:** The MDAnalysis analysis paths have a synthetic-data smoke test. Real GROMACS/VMD workflows, especially the version-dependent SASA and hydrogen-bond commands, have **not** been comprehensively validated. Verify numerical results against reference calculations before publication.

## Features

| Feature | Method / tool | Notes |
| --- | --- | --- |
| PBC processing | GROMACS `gmx trjconv` | Whole → nojump → molecule center → rotation/translation fitting |
| RMSD | MDAnalysis coordinates / Kabsch fit | Reference = first sampled frame; selected group also used for fitting |
| Radius of gyration | MDAnalysis + NumPy | Mass-weighted when masses are available, otherwise equal-mass fallback |
| RMSF (atom) | MDAnalysis + NumPy | After fitting selected group to first sampled frame |
| RMSF (residue) | MDAnalysis + NumPy | **Mean atomic/bead RMSF per residue**, *not* residue-COM RMSF |
| Minimum distance | MDAnalysis | Group A ↔ B |
| Center-of-mass distance | MDAnalysis | Group A ↔ B; group must be physically whole |
| Contact count | MDAnalysis | Atom/bead pair count within cutoff |
| Residue contact map | MDAnalysis | Fraction of sampled frames in contact; matrix and plot exports |
| PCA | NumPy SVD | Cartesian coordinates after alignment; max 2,000 selected atoms and 750 frames |
| SASA | GROMACS | Native version-dependent command; verify in local installation |
| Hydrogen bonds | GROMACS | Native version-dependent command; atomistic interpretation only |
| Visualization | VMD | Optional GUI launch |

## Installation

**Requirements:** Python 3.10+, a graphical desktop session, Python Tkinter, GROMACS for PBC/SASA/H-bonds, and optionally VMD.

Ubuntu/Debian example:

```bash
sudo apt update
sudo apt install python3-tk python3-venv
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Install GROMACS and VMD separately using your institutional package manager, HPC modules, or official distribution. Ensure `gmx` is on your `PATH`, or enter the executable location in the GUI. Tkinter is part of many Python installations but can require separate installation.

Launch:

```bash
python trajectory_gui.py
```

## Basic workflow

1. **Input:** Supply a `.tpr`, `.xtc`, optional `.ndx`, and an output folder in the main window.
2. **PBC correction (optional):** Load index groups, choose **centering**, **fitting**, and **final output** groups, and select **Run all four steps**. This generates corrected XTC stages and a matching final GRO.
3. **Analysis:** Select **Analysis workbench…** to choose analyses, a matching topology/structure and trajectory, MDAnalysis selections, frame range, stride, and (where appropriate) a contact cutoff.
4. **Batch run:** Click **Run selected analyses**. Results are written into a **single `Analysis_Results/` directory** under the selected output folder. Files are prefixed with analysis name, selection(s), and a unique timestamp/microsecond run ID to avoid overwriting earlier DNA/protein analyses.
5. **Visualize:** Plot previews appear in the workbench. The main GUI also has **Contact map…** for a dedicated contact-map window and optional VMD launching.

**Already PBC corrected?** Skip the processing step and choose your matching structure/trajectory directly in the analysis workbench. The *minimum-image PBC* checkbox affects distance/contact computations; it is **not** trajectory repair.

## Selections and groups

MDAnalysis selections use the [MDAnalysis selection language](https://docs.mdanalysis.org/stable/documentation_pages/selections.html), for example:

| Purpose | Selection A | Selection B |
| --- | --- | --- |
| Whole protein RMSD, Rg or RMSF | `protein` | not used |
| Protein Cα atom RMSF | `name CA` | not used |
| DNA phosphate-bead RMSF (SIRAH example) | `name PX` | not used |
| DNA P atom RMSF | `nucleic and name P` | not used |
| Protein–DNA contact map | `protein` | `nucleic` |
| Intramolecular residues 1–100 | `protein and resid 1:100` | `protein and resid 1:100` |
| Two components of a complex | `segid A` | `segid B` |

These expressions only work when atom names, residue names, and segment IDs match the loaded topology. `nucleic` may not recognize custom coarse-grained residue names. Contact-map pairs are aggregated to residues based on the selected atoms/beads.

**GROMACS selections are separate.** The existing `.ndx` file supports GROMACS group selection for processing/native analyses, but is **not automatically converted** into MDAnalysis selections in this release. The GROMACS group fields accept group names. The code does not yet expose separate alignment and analysis selections for MDAnalysis RMSD/RMSF.

For complexes, center and fit using scientifically appropriate groups, keep all relevant interacting molecules in the final trajectory, and verify topology–trajectory atom order. **Do not apply original full-system index numbering blindly to an atom-subset trajectory.**

## Output files

For repeated runs, the shared directory retains independent outputs:

```text
<selected-output-folder>/Analysis_Results/
├── run_YYYYMMDD_HHMMSS_ffffff.json
├── RMSD__nucleic__<run-id>__rmsd.csv
├── RMSD__nucleic__<run-id>__rmsd.png
├── RMSD__nucleic__<run-id>__rmsd.pdf
├── Radius_of_gyration__nucleic__<run-id>__rg.csv
├── RMSD__protein__<another-run-id>__rmsd.csv
├── Residue_contact_map__...__contact_map.png
└── summary_<run-id>.json
```

Outputs vary by analysis. Results include CSV data and PNG/PDF plots (plus NumPy matrices for contact maps and PCA). Metadata is saved with each run. Results and raw trajectories are intentionally excluded from Git via `.gitignore`.

## Scientific definitions and limitations

- **RMSD:** Unweighted least-squares fit of Selection A onto the first sampled frame, followed by RMSD of Selection A. Cannot currently fit on a different component and measure another.
- **Rg:** Mass-weighted radius of gyration when valid masses exist; otherwise equal-weight fallback. Inspect metadata for `equal_mass_fallback`.
- **RMSF:** Fit Selection A to the first sampled frame, calculate per-atom fluctuation around its time-averaged coordinate. Residue-wise mode averages constituent atom RMSFs, not residue-center fluctuations.
- **Contacts:** For each residue pair, occupancy is the number of sampled frames with at least one selected atom/bead pair within the cutoff, divided by the number of analyzed frames. Atom/bead **contact count** is a separate measure.
- **PBC:** Coordinate-based Rg/RMSD/RMSF/PCA require physically whole, appropriately aligned input structures. The minimum-image toggle only affects applicable distance/contact methods.
- **Coarse graining:** Select appropriate beads and cutoffs for the force field. A 4.5 Å atomistic-style cutoff is not automatically meaningful for SIRAH.
- **SASA / H-bonds:** GROMACS version-dependent flags and selection protocols may require adjustment. Conventional atomistic hydrogen-bond definitions should not be applied unchanged to coarse-grained models.
- **PCA:** Cartesian, unweighted PCA on sampled aligned coordinates; not a replica ensemble or dihedral PCA.

## Testing

The included synthetic-data smoke test covers the nine MDAnalysis analysis options and checks that output metadata and plots are produced:

```bash
python test_analysis.py
```

This is **not** a test of GROMACS processing, GROMACS SASA/H-bonds, VMD, or a complete desktop interaction. Compare at least one real trajectory against established GROMACS/MDAnalysis calculations before interpreting new systems.

## Repository contents

```text
trajectory_gui.py       Main GUI: PBC processing and VMD launch
analysis_workbench.py   Multi-analysis GUI and safe result publishing
analysis_engine.py      Numeric MDAnalysis analyses; GROMACS adapters
contact_map_gui.py      Dedicated contact-map GUI and residue occupancy
requirements.txt       Python dependencies
test_analysis.py        Synthetic-trajectory smoke test
.github/workflows/      Python syntax checks and optional smoke test
.gitignore             Excludes molecular trajectories and generated results
```

## Citing underlying software

This repository builds on [GROMACS](https://www.gromacs.org/), [MDAnalysis](https://www.mdanalysis.org/), [NumPy](https://numpy.org/), [Matplotlib](https://matplotlib.org/), and [VMD](https://www.ks.uiuc.edu/Research/vmd/). When using the tool in a publication, cite the relevant upstream software and methods as well as any repository release/DOI you publish.

## License

**No license is granted by this repository until a `LICENSE` file is added.** Choose and add a license only after confirming rights to distribute all source code and any third-party components. GitHub's repository license picker can help select an appropriate open-source license.
