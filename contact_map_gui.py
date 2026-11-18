#!/usr/bin/env python3
"""MDAnalysis residue-contact occupancy maps for the GROMACS trajectory GUI.

The calculation is intentionally separate from Tkinter so it can be tested
or reused in notebooks. Distances and cutoffs are in Angstrom.
"""
import csv
import json
from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


def residue_contacts(topology, trajectory, selection_a, selection_b, cutoff=4.5,
                     start=0, stop=None, stride=1, pbc=False, exclude_neighbors=0,
                     progress=None, cancel=None):
    """Return a fraction-of-frames contact map for two static atom selections.

    A residue pair contacts if ANY selected atom pair is <= cutoff Angstrom.
    Per-frame duplicates are deduplicated. 'stop' is exclusive. The optional
    neighbor exclusion applies only if the two selections are identical.
    """
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

    if cutoff <= 0 or stride < 1 or start < 0:
        raise ValueError('Cutoff must be > 0, stride >= 1 and start >= 0.')
    u = mda.Universe(str(topology), str(trajectory))
    a = u.select_atoms(selection_a)
    b = u.select_atoms(selection_b)
    if len(a) == 0 or len(b) == 0:
        raise ValueError('Both selections must contain at least one atom/bead.')
    if stop is None:
        stop = len(u.trajectory)
    if not 0 <= start < stop <= len(u.trajectory):
        raise ValueError(f'Frame bounds must satisfy 0 <= start < stop <= {len(u.trajectory)}.')

    ra = a.residues
    rb = b.residues
    na, nb = len(ra), len(rb)
    if na * nb > 4_000_000:
        raise ValueError(f'Map is too large ({na} x {nb}). Narrow the selections.')
    # MDAnalysis residue indices uniquely distinguish residues even where
    # resids repeat across chains or segmentation differs.
    amap = np.searchsorted(ra.resindices, a.resindices)
    bmap = np.searchsorted(rb.resindices, b.resindices)
    counts = np.zeros((na, nb), dtype=np.int32)
    frames = list(range(start, stop, stride))
    same = np.array_equal(a.indices, b.indices)
    times = []
    n_contact_pairs_per_frame = []

    for k, frame in enumerate(frames):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('Contact-map calculation cancelled.')
        ts = u.trajectory[frame]
        if pbc:
            box = ts.dimensions
            if box is None or len(box) != 6 or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0):
                raise ValueError(f'Frame {frame} has no valid periodic box. Disable PBC or use a trajectory with box information.')
        else:
            box = None
        pairs = capped_distance(a.positions, b.positions, max_cutoff=cutoff,
                                box=box, return_distances=False)
        if pairs.size:
            ii, jj = amap[pairs[:, 0]], bmap[pairs[:, 1]]
            if same:
                # Exclude same-residue contacts and optionally local sequence neighbors.
                mask = np.abs(ii - jj) > exclude_neighbors
                ii, jj = ii[mask], jj[mask]
            flat = np.unique(ii.astype(np.int64) * nb + jj.astype(np.int64))
            counts.flat[flat] += 1
            n_contact_pairs_per_frame.append(len(flat))
        else:
            n_contact_pairs_per_frame.append(0)
        times.append(float(ts.time))
        if progress is not None and (k % max(1, len(frames) // 100) == 0 or k == len(frames) - 1):
            progress(k + 1, len(frames))

    def labels(residues):
        out = []
        for residue in residues:
            seg = str(residue.segid).strip()
            label = f'{seg}:' if seg else ''
            label += f'{residue.resname}{residue.resid}'
            out.append(label)
        return out

    return {
        'matrix': counts.astype(float) / len(frames),
        'counts': counts, 'labels_a': labels(ra), 'labels_b': labels(rb),
        'resindices_a': ra.resindices.tolist(), 'resindices_b': rb.resindices.tolist(),
        'times_ps': times, 'contacts_per_frame': n_contact_pairs_per_frame,
        'n_frames': len(frames), 'n_atoms_a': len(a), 'n_atoms_b': len(b),
        'same_selection': bool(same),
    }


def save_contact_results(results, folder, metadata):
    import numpy as np
    from matplotlib.figure import Figure
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    matrix = results['matrix']
    np.save(folder / 'occupancy.npy', matrix)
    with (folder / 'occupancy.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Residue A / Residue B', *results['labels_b']])
        for name, row in zip(results['labels_a'], matrix):
            writer.writerow([name, *[f'{v:.6f}' for v in row]])
    with (folder / 'contact_counts.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Residue A / Residue B', *results['labels_b']])
        for name, row in zip(results['labels_a'], results['counts']):
            writer.writerow([name, *row.tolist()])
    with (folder / 'contacts_vs_time.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time_ps', 'contacting_residue_pairs'])
        writer.writerows(zip(results['times_ps'], results['contacts_per_frame']))
    info = {**metadata, **{key: results[key] for key in
            ('n_frames', 'n_atoms_a', 'n_atoms_b', 'same_selection',
             'labels_a', 'labels_b', 'resindices_a', 'resindices_b')}}
    (folder / 'metadata.json').write_text(json.dumps(info, indent=2))
    fig = make_contact_figure(results)
    fig.savefig(folder / 'contact_map.png', dpi=250, bbox_inches='tight')
    fig.savefig(folder / 'contact_map.pdf', bbox_inches='tight')
    return folder


def make_contact_figure(results):
    from matplotlib.figure import Figure
    mat = results['matrix']
    fig = Figure(figsize=(8, 6), dpi=100, layout='constrained')
    ax = fig.add_subplot(111)
    image = ax.imshow(mat, origin='lower', aspect='auto', cmap='viridis',
                      interpolation='nearest', vmin=0, vmax=1)
    ax.set_xlabel('Selection B — residue index')
    ax.set_ylabel('Selection A — residue index')
    ax.set_title(f'Residue contact occupancy ({results["n_frames"]} frames)')
    for axis, labels in ((ax.xaxis, results['labels_b']), (ax.yaxis, results['labels_a'])):
        n = len(labels)
        if n <= 25:
            locs = list(range(n))
        else:
            import numpy as np
            locs = sorted(set(np.linspace(0, n - 1, min(12, n), dtype=int).tolist()))
        axis.set_ticks(locs)
        axis.set_ticklabels([labels[i] for i in locs], rotation=90 if axis is ax.xaxis else 0,
                            fontsize=7)
    fig.colorbar(image, ax=ax, label='Fraction of frames', shrink=0.85)
    return fig


class ContactMapWindow:
    def __init__(self, master, structure='', trajectory=''):
        self.win = tk.Toplevel(master)
        self.win.title('MDAnalysis — residue contact map')
        self.win.geometry('1050x850')
        self.win.minsize(760, 580)
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.running = False
        self.canvas = None
        self.figure = None
        self.values = {}
        defaults = dict(topology=structure, trajectory=trajectory, selection_a='protein',
                        selection_b='nucleic', cutoff='4.5', start='0', stop='',
                        stride='10', exclude_neighbors='3', folder='', pbc=False)
        for key, val in defaults.items():
            self.values[key] = tk.BooleanVar(value=val) if key == 'pbc' else tk.StringVar(value=val)

        main = ttk.Frame(self.win, padding=12)
        main.pack(fill='both', expand=True)
        main.columnconfigure(1, weight=1)
        inputs = [('topology', 'Topology / structure (matching GRO, PDB, TPR)'),
                  ('trajectory', 'Trajectory (XTC, DCD, ...)'),
                  ('selection_a', 'Selection A (MDAnalysis syntax)'),
                  ('selection_b', 'Selection B (MDAnalysis syntax)'),
                  ('cutoff', 'Contact cutoff (Å)'),
                  ('start', 'Start frame (0-based)'),
                  ('stop', 'Stop frame (exclusive; blank = end)'),
                  ('stride', 'Analyze every Nth frame'),
                  ('exclude_neighbors', 'Same-selection: exclude ±N residue neighbors'),
                  ('folder', 'Output folder (blank = next to trajectory)')]
        for i, (key, label) in enumerate(inputs):
            ttk.Label(main, text=label).grid(row=i, column=0, sticky='w', padx=(0, 10), pady=3)
            ttk.Entry(main, textvariable=self.values[key]).grid(row=i, column=1, sticky='ew', pady=3)
            if key in ('topology', 'trajectory', 'folder'):
                ttk.Button(main, text='Browse', command=lambda k=key: self.browse(k)).grid(
                    row=i, column=2, padx=(8, 0))
        row = len(inputs)
        ttk.Checkbutton(main, text='Apply minimum-image PBC distances (requires a valid box)',
                        variable=self.values['pbc']).grid(row=row, column=0, columnspan=3, sticky='w', pady=6)
        row += 1
        ttk.Label(main, text='Contact = at least one selected atom/bead pair within cutoff in a frame. '
                  'Select a matched structure and trajectory; group labels depend on topology.').grid(
                  row=row, column=0, columnspan=3, sticky='w', pady=(0, 6))
        row += 1
        bar = ttk.Frame(main)
        bar.grid(row=row, column=0, columnspan=3, sticky='ew', pady=4)
        self.run_button = ttk.Button(bar, text='Calculate contact map', command=self.run)
        self.run_button.pack(side='left')
        self.stop_button = ttk.Button(bar, text='Cancel', command=self.stop, state='disabled')
        self.stop_button.pack(side='left', padx=8)
        self.progress = ttk.Progressbar(bar, maximum=100, length=200)
        self.progress.pack(side='left', padx=10, fill='x', expand=True)
        self.status = tk.StringVar(value='Ready')
        ttk.Label(main, textvariable=self.status).grid(row=row + 1, column=0, columnspan=3, sticky='w')
        self.plot_frame = ttk.Frame(main)
        self.plot_frame.grid(row=row + 2, column=0, columnspan=3, sticky='nsew', pady=8)
        main.rowconfigure(row + 2, weight=1)
        self.win.after(100, self.poll)

    def browse(self, key):
        if key == 'folder':
            value = filedialog.askdirectory(parent=self.win)
        else:
            value = filedialog.askopenfilename(parent=self.win)
        if value:
            self.values[key].set(value)

    def run(self):
        if self.running:
            return
        cfg = {k: var.get() for k, var in self.values.items()}
        try:
            for key in ('topology', 'trajectory'):
                if not Path(cfg[key]).is_file():
                    raise ValueError(f'Choose an existing {key} file.')
            cfg['cutoff'] = float(cfg['cutoff'])
            cfg['start'] = int(cfg['start'])
            cfg['stop'] = int(cfg['stop']) if cfg['stop'].strip() else None
            cfg['stride'] = int(cfg['stride'])
            cfg['exclude_neighbors'] = int(cfg['exclude_neighbors'])
            if cfg['cutoff'] <= 0 or cfg['start'] < 0 or cfg['stride'] < 1 or cfg['exclude_neighbors'] < 0:
                raise ValueError('Cutoff must be positive; frame start and neighbor exclusion nonnegative; stride >= 1.')
            if not cfg['selection_a'].strip() or not cfg['selection_b'].strip():
                raise ValueError('Both atom selections are required.')
        except Exception as e:
            messagebox.showerror('Invalid input', str(e), parent=self.win)
            return
        self.running = True
        self.cancel.clear()
        self.run_button.configure(state='disabled')
        self.stop_button.configure(state='normal')
        self.progress['value'] = 0
        self.status.set('Calculating contacts...')

        def worker():
            try:
                results = residue_contacts(cfg['topology'], cfg['trajectory'],
                    cfg['selection_a'], cfg['selection_b'], cutoff=cfg['cutoff'],
                    start=cfg['start'], stop=cfg['stop'], stride=cfg['stride'],
                    pbc=cfg['pbc'], exclude_neighbors=cfg['exclude_neighbors'],
                    cancel=self.cancel,
                    progress=lambda n, total: self.events.put(('progress', 100 * n / total)))
                if self.cancel.is_set():
                    raise RuntimeError('Contact-map calculation cancelled.')
                parent = Path(cfg['folder']).expanduser() if cfg['folder'] else Path(cfg['trajectory']).parent
                from datetime import datetime
                destination = parent / ('contact_map_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
                folder = save_contact_results(results, destination, cfg)
                self.events.put(('done', (results, str(folder))))
            except Exception as exc:
                self.events.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    def stop(self):
        self.cancel.set()
        self.status.set('Cancelling after the current frame...')

    def poll(self):
        if not self.win.winfo_exists():
            return
        for _ in range(100):
            try:
                kind, data = self.events.get_nowait()
            except queue.Empty:
                break
            if kind == 'progress':
                self.progress['value'] = data
            else:
                self.running = False
                self.run_button.configure(state='normal')
                self.stop_button.configure(state='disabled')
                if kind == 'error':
                    self.status.set('Failed: ' + data)
                    messagebox.showerror('Contact map', data, parent=self.win)
                elif kind == 'done':
                    results, folder = data
                    self.status.set(f'Saved {results["n_frames"]} frames to {folder}')
                    self.show_plot(results)
        self.win.after(100, self.poll)

    def show_plot(self, results):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        self.figure = make_contact_figure(results)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        NavigationToolbar2Tk(self.canvas, self.plot_frame).update()
