# Contributing

Please open an issue describing the use case, inputs (without private or unpublished simulation data), expected outcome, actual outcome, and GROMACS/MDAnalysis/Python versions.

For changes, use a feature branch, retain output-name collision protection, and run:

```bash
python -m compileall -q trajectory_gui.py analysis_workbench.py analysis_engine.py contact_map_gui.py test_analysis.py
python test_analysis.py
```

The synthetic test is not a substitute for physical/numerical validation. Add a regression test and state assumptions about PBC, fitting groups, atom/bead selections, units, and contact definitions whenever changing analysis code.

Avoid committing trajectories, unpublished structures, credentials, cluster paths, or large result files.
