"""Analysis back-end for the GROMACS/MDAnalysis workbench.

Units: MDAnalysis positions Angstrom, time ps; RMSD/Rg/RMSF/distance Angstrom.
Independent of Tk; all outputs are written to a distinct per-run directory.
"""
import csv
import json
import re
import subprocess
import threading
from pathlib import Path
import numpy as np


ANALYSES = ('RMSD', 'Radius of gyration', 'RMSF by residue', 'RMSF by atom',
            'Group COM distance', 'Minimum distance', 'Contact count',
            'Residue contact map', 'PCA', 'SASA (GROMACS)', 'Hydrogen bonds (GROMACS)')


def _frames(u, cfg):
    start, stride = int(cfg['start']), int(cfg['stride'])
    stop = int(cfg['stop']) if cfg.get('stop') not in (None, '') else len(u.trajectory)
    if start < 0 or stride < 1 or not 0 <= start < stop <= len(u.trajectory):
        raise ValueError(f'Invalid frames: use 0 <= start < stop <= {len(u.trajectory)} and stride >= 1.')
    return range(start, stop, stride)


def _check_cancel(cancel):
    if cancel is not None and cancel.is_set():
        raise RuntimeError('Analysis cancelled.')


def _selection(u, expr):
    atoms = u.select_atoms(expr)
    if not len(atoms):
        raise ValueError(f'Selection contains no atoms: {expr!r}')
    return atoms


def _kabsch(mobile, reference):
    """Return aligned mobile (Nx3) onto reference (Nx3)."""
    m0, r0 = mobile.mean(axis=0), reference.mean(axis=0)
    x, y = mobile - m0, reference - r0
    U, _, Vt = np.linalg.svd(x.T @ y)
    adjust = np.eye(3)
    adjust[-1, -1] = np.linalg.det(U @ Vt)
    return x @ (U @ adjust @ Vt) + r0


def _align_all(coords, fit_indices, fit_reference):
    mobile = coords[fit_indices]
    center = mobile.mean(axis=0)
    centered = mobile - center
    target = fit_reference
    U, _, Vt = np.linalg.svd(centered.T @ (target - target.mean(axis=0)))
    adjust = np.eye(3)
    adjust[-1, -1] = np.linalg.det(U @ Vt)
    return (coords - center) @ (U @ adjust @ Vt) + target.mean(axis=0)


def _save_series(folder, key, x, y, xlabel, ylabel, xlabels=None, extra=None):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    folder = Path(folder)
    with (folder / f'{key}.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow([xlabel, ylabel])
        if xlabels is None:
            writer.writerows(zip(x, y))
        else:
            writer.writerows(zip(xlabels, y))
    fig = Figure(figsize=(8, 4.7), dpi=110, layout='constrained')
    ax = fig.add_subplot(111)
    ax.plot(x, y, '-', lw=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(key.replace('_', ' '))
    ax.grid(alpha=0.25)
    fig.savefig(folder / f'{key}.png', dpi=220)
    fig.savefig(folder / f'{key}.pdf')
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig).draw()
    report = {'output': f'{key}.csv', 'plot': f'{key}.png', 'n_points': len(y),
              'mean': float(np.mean(y)), 'std': float(np.std(y))}
    if extra:
        report.update(extra)
    return report


def _get_masses(atoms):
    try:
        masses = atoms.masses
        if (not np.all(np.isfinite(masses))) or np.any(masses <= 0):
            raise ValueError('Invalid masses')
        return masses, False
    except (ValueError, AttributeError, Exception):
        return np.ones(len(atoms)), True


def _gmx(cfg, analysis, folder, cancel=None):
    """Use GROMACS native implementations where MDAnalysis is not reliable.

    GROMACS selections are index-group names, not MDAnalysis selection strings.
    """
    gmx = cfg.get('gmx', 'gmx').strip() or 'gmx'
    tpr = cfg.get('tpr', '').strip()
    ndx = cfg.get('ndx', '').strip()
    if not tpr or not Path(tpr).is_file():
        raise ValueError('GROMACS SASA/H-bond analysis requires a valid TPR.')
    flags = [gmx]
    if analysis == 'SASA (GROMACS)':
        flags += ['sasa', '-s', tpr, '-f', cfg['trajectory'], '-o', str(folder / 'sasa.xvg')]
        selection_text = cfg.get('gmx_selection', 'Protein') + '\n'
    else:
        flags += ['hbond', '-s', tpr, '-f', cfg['trajectory'], '-num', str(folder / 'hbond.xvg')]
        selection_text = cfg.get('gmx_selection', 'Protein') + '\n' + cfg.get('gmx_selection_b', 'Protein') + '\n'
    if ndx:
        flags += ['-n', ndx]
    begin, end = cfg.get('begin_ps', ''), cfg.get('end_ps', '')
    if begin != '':
        flags += ['-b', str(float(begin))]
    if end != '':
        flags += ['-e', str(float(end))]
    if cfg.get('dt_ps', '') != '':
        flags += ['-dt', str(float(cfg['dt_ps']))]
    try:
        completed = subprocess.run(flags, input=selection_text, capture_output=True, text=True,
                                   errors='replace', timeout=3600, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f'GROMACS not found: {gmx}') from exc
    (folder / ('sasa_gromacs_command.json' if analysis == 'SASA (GROMACS)' else 'hbond_gromacs_command.json')).write_text(json.dumps({'args': flags, 'stdin': selection_text}, indent=2))
    (folder / ('sasa_gromacs.log' if analysis == 'SASA (GROMACS)' else 'hbond_gromacs.log')).write_text(completed.stdout + '\n' + completed.stderr)
    if completed.returncode:
        raise RuntimeError(f'GROMACS exited with code {completed.returncode}. See gromacs.log. Check version-specific selection prompts.')
    path = folder / ('sasa.xvg' if analysis == 'SASA (GROMACS)' else 'hbond.xvg')
    if not path.exists():
        raise RuntimeError('GROMACS finished without expected XVG output. See gromacs.log.')
    rows = []
    for line in path.read_text(errors='replace').splitlines():
        if not line.strip() or line.lstrip().startswith(('#', '@')):
            continue
        fields = line.split()
        try:
            rows.append([float(v) for v in fields])
        except ValueError:
            pass
    if not rows:
        raise RuntimeError('No numeric XVG rows found.')
    data = np.asarray(rows)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError('XVG needs at least two numeric columns.')
    units = 'nm²' if analysis == 'SASA (GROMACS)' else 'count'
    result = _save_series(folder, 'sasa' if analysis.startswith('SASA') else 'hydrogen_bonds',
                          data[:, 0], data[:, 1], 'Time (ps)', units)
    result['note'] = 'Native GROMACS analysis; selections and command interface depend on GROMACS version.'
    return result


def run_analysis(analysis, cfg, folder, progress=None, cancel=None):
    """Execute one analysis, write CSV/PNG/PDF and metadata, and return statistics."""
    if analysis not in ANALYSES:
        raise ValueError(f'Unsupported analysis: {analysis}')
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    try:
        if analysis.endswith('(GROMACS)'):
            result = _gmx(cfg, analysis, folder, cancel)
        elif analysis == 'Residue contact map':
            from contact_map_gui import residue_contacts, save_contact_results
            out = residue_contacts(cfg['topology'], cfg['trajectory'], cfg['selection_a'], cfg['selection_b'],
                                   cutoff=float(cfg['cutoff']), start=int(cfg['start']),
                                   stop=int(cfg['stop']) if cfg.get('stop') else None,
                                   stride=int(cfg['stride']), pbc=bool(cfg.get('pbc')),
                                   exclude_neighbors=int(cfg.get('neighbors', 0)),
                                   progress=progress, cancel=cancel)
            save_contact_results(out, folder, cfg)
            result = {'plot': 'contact_map.png', 'output': 'occupancy.csv', 'n_frames': out['n_frames'],
                      'mean_occupancy': float(np.mean(out['matrix']))}
        else:
            import MDAnalysis as mda
            from MDAnalysis.lib.distances import capped_distance, distance_array
            u = mda.Universe(cfg['topology'], cfg['trajectory'])
            a = _selection(u, cfg['selection_a'])
            b = None
            if analysis in ('Group COM distance', 'Minimum distance', 'Contact count'):
                b = _selection(u, cfg['selection_b'])
                if analysis in ('Minimum distance', 'Contact count') and np.intersect1d(a.indices, b.indices).size:
                    raise ValueError('Minimum distance/contact count require disjoint selections to avoid trivial self-contacts.')
            frame_ids = _frames(u, cfg)
            times, values = [], []
            first = u.trajectory[frame_ids[0]]
            refpos = a.positions.astype(float).copy()
            # For RMSD and RMSF, fit on the selected atoms; ref is first selected frame.
            if analysis == 'RMSD':
                for k, i in enumerate(frame_ids):
                    _check_cancel(cancel)
                    ts = u.trajectory[i]
                    aligned = _kabsch(a.positions.astype(float), refpos)
                    times.append(float(ts.time))
                    values.append(float(np.sqrt(np.mean(np.sum((aligned-refpos)**2, axis=1)))))
                    if progress:
                        progress(k+1, len(frame_ids))
                result = _save_series(folder, 'rmsd', times, values, 'Time (ps)', 'RMSD (Å)')
            elif analysis == 'Radius of gyration':
                masses, fallback = _get_masses(a)
                for k, i in enumerate(frame_ids):
                    _check_cancel(cancel)
                    ts = u.trajectory[i]
                    coords = a.positions.astype(float)
                    center = np.average(coords, axis=0, weights=masses)
                    rg = np.sqrt(np.average(np.sum((coords-center)**2, axis=1), weights=masses))
                    times.append(float(ts.time)); values.append(float(rg))
                    if progress: progress(k+1, len(frame_ids))
                result = _save_series(folder, 'rg', times, values, 'Time (ps)', 'Rg (Å)',
                                      extra={'equal_mass_fallback': fallback})
            elif analysis.startswith('RMSF'):
                # Align on selection A then calculate atomic fluctuations about time-average positions.
                mean = np.zeros_like(refpos)
                accum = np.zeros(len(a), dtype=float)
                n = 0
                for k, i in enumerate(frame_ids):
                    _check_cancel(cancel)
                    u.trajectory[i]
                    p = _kabsch(a.positions.astype(float), refpos)
                    n += 1
                    delta = p - mean
                    mean += delta / n
                    accum += np.sum(delta * (p - mean), axis=1)
                    if progress: progress(k+1, len(frame_ids))
                rmsf = np.sqrt(np.maximum(accum / n, 0))
                if analysis == 'RMSF by residue':
                    ids = a.resindices
                    unique = np.unique(ids)
                    vals = [float(np.mean(rmsf[ids == rid])) for rid in unique]
                    labels = [f'{u.residues[rid].segid}:{u.residues[rid].resname}{u.residues[rid].resid}' for rid in unique]
                    result = _save_series(folder, 'rmsf_residue', np.arange(1, len(vals)+1), vals,
                                          'Residue (ordinal)', 'Mean atomic RMSF (Å)', xlabels=labels)
                    result['note'] = 'Residue-wise mean of selected atom/bead RMSFs (not residue COM RMSF).'
                else:
                    result = _save_series(folder, 'rmsf_atom', np.arange(1, len(a)+1), rmsf,
                                          'Selected atom (ordinal)', 'RMSF (Å)',
                                          xlabels=[f'{u.residues[idx].resid}:{name}' for idx, name in zip(a.resindices, a.names)])
            elif analysis in ('Group COM distance', 'Minimum distance', 'Contact count'):
                cutoff = float(cfg['cutoff'])
                if cutoff <= 0: raise ValueError('Contact cutoff must be > 0 Å.')
                for k, i in enumerate(frame_ids):
                    _check_cancel(cancel)
                    ts = u.trajectory[i]
                    box = ts.dimensions if cfg.get('pbc') else None
                    if cfg.get('pbc') and (box is None or np.any(np.asarray(box[:3]) <= 0)):
                        raise ValueError('PBC requested but box dimensions are invalid.')
                    if analysis == 'Group COM distance':
                        m1, _ = _get_masses(a); m2, _ = _get_masses(b)
                        aa = np.average(a.positions, weights=m1, axis=0)
                        bb = np.average(b.positions, weights=m2, axis=0)
                        v = float(distance_array(aa[None, :], bb[None, :], box=box)[0, 0])
                    elif analysis == 'Minimum distance':
                        # distance_array has potentially quadratic memory; chunk on A.
                        v = min(float(np.min(distance_array(a.positions[j:j+256], b.positions, box=box)))
                                for j in range(0, len(a), 256))
                    else:
                        pairs = capped_distance(a.positions, b.positions, max_cutoff=cutoff,
                                                box=box, return_distances=False)
                        # Unique atom pairs. Identical groups count unordered pairs once, excluding self-pairs.
                        if np.array_equal(a.indices, b.indices):
                            pairs = pairs[pairs[:, 0] < pairs[:, 1]]
                        v = len(pairs)
                    times.append(float(ts.time)); values.append(v)
                    if progress: progress(k+1, len(frame_ids))
                name = {'Group COM distance':'com_distance','Minimum distance':'minimum_distance',
                        'Contact count':'contact_count'}[analysis]
                result = _save_series(folder, name, times, values, 'Time (ps)',
                                      'Number of atom/bead pairs' if analysis == 'Contact count' else 'Distance (Å)')
                if analysis == 'Group COM distance' and cfg.get('pbc'):
                    result['note'] = 'Distance between group COMs with minimum-image convention; group coordinates must be made whole first.'
            elif analysis == 'PCA':
                if len(a) > 2000 or len(frame_ids) > 750:
                    raise ValueError('PCA safety limit: <=2000 selected atoms and <=750 frames. Narrow selection or increase stride.')
                coords = []
                for k, i in enumerate(frame_ids):
                    _check_cancel(cancel)
                    u.trajectory[i]
                    coords.append(_kabsch(a.positions.astype(float), refpos).ravel())
                    if progress: progress(k+1, len(frame_ids))
                data = np.asarray(coords)
                if len(data) < 3: raise ValueError('PCA requires at least three frames.')
                centered = data - data.mean(axis=0)
                _, singular, vt = np.linalg.svd(centered, full_matrices=False)
                projected = centered @ vt[:2].T
                variance = singular**2 / (len(data)-1)
                explained = variance / variance.sum() if variance.sum() else np.zeros_like(variance)
                with (folder / 'pca_projection.csv').open('w', newline='') as f:
                    writer = csv.writer(f); writer.writerow(['time_ps', 'PC1_A', 'PC2_A'])
                    writer.writerows(zip([float(u.trajectory[i].time) for i in frame_ids], projected[:, 0], projected[:, 1]))
                np.save(folder / 'pca_components.npy', vt[:2])
                np.savetxt(folder / 'pca_variance.csv', np.column_stack((np.arange(1,len(explained)+1), explained)),
                           delimiter=',', header='component,variance_fraction', comments='')
                from matplotlib.figure import Figure
                fig = Figure(figsize=(7, 5), layout='constrained'); ax = fig.add_subplot(111)
                ax.scatter(projected[:,0], projected[:,1], c=np.arange(len(projected)), cmap='viridis', s=8)
                ax.set(xlabel='PC1 (Å)', ylabel='PC2 (Å)', title='PCA projection (aligned selection)')
                fig.savefig(folder / 'pca.png', dpi=220); fig.savefig(folder / 'pca.pdf')
                result = {'output': 'pca_projection.csv', 'plot': 'pca.png',
                          'pc1_fraction': float(explained[0]), 'pc2_fraction': float(explained[1])}
            else:
                raise ValueError('Unimplemented analysis')
        (folder / ('metadata_' + re.sub(r'[^a-z0-9]+', '_', analysis.lower()).strip('_') + '.json')).write_text(json.dumps({'analysis': analysis, 'settings': cfg,
                                                         'result': result}, indent=2, default=str))
        return result
    except Exception:
        (folder / ('FAILED_' + re.sub(r'[^a-z0-9]+', '_', analysis.lower()).strip('_') + '.txt')).write_text('Analysis failed. No result should be interpreted from this folder.\n')
        raise
