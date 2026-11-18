"""Synthetic, fast, non-GROMACS integration smoke test."""
import tempfile
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from analysis_engine import run_analysis


def test_all():
    root = Path(tempfile.mkdtemp(prefix='gmx_workbench_test_'))
    u = mda.Universe.empty(6, n_residues=3, atom_resindex=np.array([0,0,1,1,2,2]), trajectory=True)
    u.add_TopologyAttr('names', ['C','N','C','N','C','N'])
    u.add_TopologyAttr('resnames', ['ALA','GLY','LYS'])
    u.add_TopologyAttr('resids', [1,2,3])
    u.add_TopologyAttr('masses', [12.,14.,12.,14.,12.,14.])
    pdb = root/'test.pdb'; dcd = root/'test.dcd'
    coords = np.array([[0,0,0],[0,0,1],[4,0,0],[4,0,1],[9,0,0],[9,0,1]], dtype=float)
    u.atoms.positions = coords
    with mda.Writer(str(pdb), n_atoms=6) as w:
        w.write(u.atoms)
    with mda.Writer(str(dcd), n_atoms=6) as w:
        for i in range(6):
            positions = coords.copy()
            positions[4:,0] = 9 + i
            positions[2:4,1] = 0.1*i
            u.atoms.positions = positions
            w.write(u.atoms)
    cfg = dict(topology=str(pdb), trajectory=str(dcd), selection_a='all',
               selection_b='resid 3', start=0, stop=None, stride=1,
               cutoff=5.0, neighbors=0, pbc=False)
    for label in ('RMSD', 'Radius of gyration', 'RMSF by residue', 'RMSF by atom',
                  'Group COM distance', 'Minimum distance', 'Contact count',
                  'Residue contact map', 'PCA'):
        dest = root / label.replace(' ', '_')
        cfg['selection_a'] = 'resid 1:2' if label in ('Minimum distance', 'Contact count') else 'all'
        res = run_analysis(label, cfg, dest)
        assert (dest / 'metadata.json').exists(), label
        assert (dest / res['plot']).exists(), label
        print(label, res)
    return root

if __name__ == '__main__':
    print('RESULT:', test_all())
