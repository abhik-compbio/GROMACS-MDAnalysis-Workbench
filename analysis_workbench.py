"""Analysis launcher, concurrent-safe Tk UI. Heavy calculations never run on UI thread."""
import json
import queue
import re
import tempfile
import threading
from datetime import datetime
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from analysis_engine import ANALYSES, run_analysis


def _safe_tag(value, limit=42):
    """Readable, filename-safe selection tag (e.g. nucleic, protein, name_PX)."""
    tag = re.sub(r'[^A-Za-z0-9]+', '_', str(value)).strip('_')
    return (tag[:limit].rstrip('_') or 'selection')


def publish_analysis_outputs(staging, destination, label):
    """Promote a completed analysis atomically into a shared, flat directory.

    Every file carries an analysis/selection/run prefix; nothing is overwritten.
    Uses a unique staging path, so failed analyses cannot destroy prior results.
    """
    staging, destination = Path(staging), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    mappings = {}
    for source in sorted(staging.iterdir()):
        if not source.is_file():
            continue
        desired = destination / f'{label}__{source.name}'
        # Defensive collision guard even if two runs receive identical timestamps.
        counter = 2
        while desired.exists():
            desired = destination / f'{label}_{counter}__{source.name}'
            counter += 1
        source.replace(desired)
        mappings[source.name] = desired.name
    return mappings


class AnalysisWorkbench:
    def __init__(self, master, topology='', trajectory='', tpr='', ndx='', gmx='gmx'):
        self.win = tk.Toplevel(master)
        self.win.title('Trajectory analysis workbench')
        self.win.geometry('1150x880')
        self.win.minsize(850, 650)
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.running = False
        self.plot_canvas = None
        self.fig = None
        self.vars = {}
        defaults = {
            'topology': topology or tpr,
            'trajectory': trajectory,
            'tpr': tpr,
            'ndx': ndx,
            'gmx': gmx or 'gmx',
            'selection_a': 'protein',
            'selection_b': 'nucleic',
            'cutoff': '4.5',
            'start': '0', 'stop': '', 'stride': '10', 'neighbors': '3',
            'pbc': False,
            'folder': str(Path(trajectory).parent) if trajectory else '',
            'gmx_selection': 'Protein', 'gmx_selection_b': 'Protein',
            'begin_ps': '', 'end_ps': '', 'dt_ps': ''
        }
        for k, val in defaults.items():
            self.vars[k] = tk.BooleanVar(value=val) if isinstance(val, bool) else tk.StringVar(value=val)
        shell = ttk.Frame(self.win, padding=12)
        shell.pack(fill='both', expand=True)
        shell.columnconfigure(0, weight=2)
        shell.columnconfigure(1, weight=3)
        shell.rowconfigure(1, weight=1)
        desc = ttk.Label(shell, text='GROMACS + MDAnalysis | Select multiple analyses and click Run selected',
                         font=('Sans', 12, 'bold'))
        desc.grid(row=0, column=0, columnspan=2, sticky='w', pady=(0, 10))
        left = ttk.Frame(shell)
        left.grid(row=1, column=0, sticky='nsew', padx=(0, 12))
        left.columnconfigure(1, weight=1)
        self._form(left)
        right = ttk.Frame(shell)
        right.grid(row=1, column=1, sticky='nsew')
        right.columnconfigure(0, weight=1)
        right.rowconfigure(2, weight=1)
        ttk.Label(right, text='Analyses (multiselect)', font=('Sans', 11, 'bold')).grid(row=0, column=0, sticky='w')
        listbox_frame = ttk.Frame(right)
        listbox_frame.grid(row=1, column=0, sticky='ew', pady=(4, 8))
        self.analysis_box = tk.Listbox(listbox_frame, height=len(ANALYSES), selectmode='extended', exportselection=False)
        for label in ANALYSES:
            self.analysis_box.insert('end', label)
        self.analysis_box.pack(fill='x')
        for i in (0, 1, 2, 7):
            self.analysis_box.selection_set(i)
        preview = ttk.LabelFrame(right, text='Results (PNG preview)', padding=6)
        preview.grid(row=2, column=0, sticky='nsew')
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(1, weight=1)
        self.results_box = ttk.Combobox(preview, state='readonly')
        self.results_box.grid(row=0, column=0, sticky='ew')
        self.results_box.bind('<<ComboboxSelected>>', self.show_result)
        self.plot_holder = ttk.Frame(preview)
        self.plot_holder.grid(row=1, column=0, sticky='nsew')
        bottom = ttk.Frame(shell)
        bottom.grid(row=2, column=0, columnspan=2, sticky='ew', pady=(10, 4))
        self.run_button = ttk.Button(bottom, text='Run selected analyses', command=self.run)
        self.run_button.pack(side='left')
        self.stop_button = ttk.Button(bottom, text='Cancel', command=self.stop, state='disabled')
        self.stop_button.pack(side='left', padx=8)
        self.progress = ttk.Progressbar(bottom, length=200)
        self.progress.pack(side='left', padx=8, fill='x', expand=True)
        self.status = tk.StringVar(value='Ready. GROMACS commands require an installed gmx executable.')
        ttk.Label(shell, textvariable=self.status, wraplength=950).grid(row=3, column=0, columnspan=2, sticky='w')
        self.log = tk.Text(shell, height=7, wrap='word')
        self.log.grid(row=4, column=0, columnspan=2, sticky='nsew', pady=(8, 0))
        shell.rowconfigure(4, weight=0)
        self.completed = {}
        self.win.after(100, self.poll)

    def _form(self, parent):
        rows = [
            ('topology', 'Structure/topology (GRO/PDB/TPR)'),
            ('trajectory', 'Trajectory (XTC/TRR/DCD)'),
            ('tpr', 'GROMACS TPR (SASA/H-bond)'),
            ('ndx', 'GROMACS NDX (optional)'),
            ('gmx', 'GROMACS executable'),
            ('selection_a', 'MDAnalysis selection A'),
            ('selection_b', 'MDAnalysis selection B'),
            ('cutoff', 'Contact cutoff (Å)'),
            ('start', 'Start frame (inclusive)'),
            ('stop', 'Stop frame (exclusive; blank=end)'),
            ('stride', 'Frame stride'),
            ('neighbors', 'Map: exclude ±N neighboring residues'),
            ('folder', 'Output parent folder'),
            ('gmx_selection', 'GROMACS group 1 (SASA/H-bond)'),
            ('gmx_selection_b', 'GROMACS group 2 (H-bond)'),
            ('begin_ps', 'GROMACS start time (ps, optional)'),
            ('end_ps', 'GROMACS end time (ps, optional)'),
            ('dt_ps', 'GROMACS time step (ps, optional)')]
        for i, (key, label) in enumerate(rows):
            ttk.Label(parent, text=label).grid(row=i, column=0, sticky='w', pady=2)
            ttk.Entry(parent, textvariable=self.vars[key]).grid(row=i, column=1, sticky='ew', padx=4, pady=2)
            if key in ('topology', 'trajectory', 'tpr', 'ndx', 'folder'):
                ttk.Button(parent, text='…', width=3, command=lambda k=key: self.browse(k)).grid(row=i, column=2)
        ttk.Checkbutton(parent, text='Use minimum-image PBC for distances/contacts',
                        variable=self.vars['pbc']).grid(row=len(rows), column=0, columnspan=3, sticky='w', pady=5)
        ttk.Label(parent, text='MDAnalysis uses selections such as "protein", "nucleic", "segid A", '
                  '"resid 1:50". The GROMACS group fields require index-group names.\n'
                  'RMSD/RMSF: selection A is also the fitting group; reference = first analyzed frame.\n'
                  'Contact maps: cutoff is Å and output is fraction of frames.\n'
                  'Input trajectories should already be made whole for COM/Rg/RMSD; '
                  'PBC distance mode does not repair broken molecules.', wraplength=450,
                  justify='left').grid(row=len(rows)+1, column=0, columnspan=3, sticky='w', pady=8)

    def browse(self, key):
        if self.running:
            return
        path = filedialog.askdirectory(parent=self.win) if key == 'folder' else filedialog.askopenfilename(parent=self.win)
        if path:
            self.vars[key].set(path)

    def _settings(self, selections):
        cfg = {k: v.get() for k, v in self.vars.items()}
        for key in ('trajectory', 'topology'):
            if any(not a.endswith('(GROMACS)') for a in selections) and not Path(cfg[key]).is_file():
                raise ValueError(f'Invalid {key}: {cfg[key]}')
        if any(a.endswith('(GROMACS)') for a in selections):
            if not Path(cfg['trajectory']).is_file() or not Path(cfg['tpr']).is_file():
                raise ValueError('GROMACS analyses need valid trajectory and TPR files.')
        try:
            cfg['start'] = int(cfg['start'])
            cfg['stride'] = int(cfg['stride'])
            cfg['neighbors'] = int(cfg['neighbors'])
            cfg['cutoff'] = float(cfg['cutoff'])
            cfg['stop'] = int(cfg['stop']) if cfg['stop'] else None
        except ValueError as e:
            raise ValueError('Invalid frame range, stride, neighbor exclusion or cutoff.') from e
        if cfg['start'] < 0 or cfg['stride'] < 1 or cfg['neighbors'] < 0 or cfg['cutoff'] <= 0:
            raise ValueError('Start/neighbors >= 0, stride >= 1, cutoff > 0.')
        for selection in ('selection_a', 'selection_b'):
            if not cfg[selection].strip():
                raise ValueError(f'{selection} cannot be blank.')
        folder = Path(cfg['folder']).expanduser() if cfg['folder'] else Path(cfg['trajectory']).parent
        folder.mkdir(parents=True, exist_ok=True)
        return cfg, folder

    def run(self):
        if self.running:
            return
        selected = [ANALYSES[i] for i in self.analysis_box.curselection()]
        if not selected:
            messagebox.showerror('Analyses', 'Select at least one analysis.', parent=self.win)
            return
        try:
            cfg, parent = self._settings(selected)
        except Exception as e:
            messagebox.showerror('Invalid settings', str(e), parent=self.win)
            return
        # One permanent results directory; distinct run IDs prevent data loss.
        run_dir = parent / 'Analysis_Results'
        run_dir.mkdir(parents=True, exist_ok=True)
        run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        (run_dir / f'run_{run_id}.json').write_text(
            json.dumps({'run_id': run_id, 'selected': selected, 'settings': cfg}, indent=2))
        # Keep existing completed previews instead of clearing them on every run.
        self.running = True
        self.cancel.clear()
        self.run_button.configure(state='disabled')
        self.stop_button.configure(state='normal')
        self.progress['value'] = 0
        self.progress['maximum'] = len(selected) * 100
        self.status.set(f'Running {len(selected)} analyses; output: {run_dir}')
        self.log.insert('end', f'\nRun folder: {run_dir}\n')
        def worker():
            summary = {}
            for idx, analysis in enumerate(selected):
                if self.cancel.is_set():
                    summary[analysis] = {'status': 'cancelled'}
                    break
                self.events.put(('status', f'{idx+1}/{len(selected)}: {analysis}'))
                is_pair = analysis in ('Residue contact map', 'Group COM distance',
                                       'Minimum distance', 'Contact count', 'Hydrogen bonds (GROMACS)')
                group_a = cfg['gmx_selection'] if analysis.endswith('(GROMACS)') else cfg['selection_a']
                group_b = cfg['gmx_selection_b'] if analysis.endswith('(GROMACS)') else cfg['selection_b']
                group_tag = _safe_tag(group_a) + ('_vs_' + _safe_tag(group_b) if is_pair else '')
                kind_tag = _safe_tag(analysis)
                label = f'{kind_tag}__{group_tag}__{run_id}'
                try:
                    # Do not write unprefixed rmsd.csv / rg.csv into a shared folder:
                    # stage separately, then move successful outputs with unique names.
                    with tempfile.TemporaryDirectory(prefix='.analysis_tmp_', dir=run_dir) as staging:
                        result = run_analysis(analysis, cfg, staging,
                                              progress=lambda n, total, j=idx: self.events.put(('progress', (j+n/total)*100)),
                                              cancel=self.cancel)
                        mapping = publish_analysis_outputs(staging, run_dir, label)
                    result = dict(result)
                    for key in ('plot', 'output'):
                        if key in result and result[key] in mapping:
                            result[key] = mapping[result[key]]
                    # Keep metadata's advertised output filenames synchronized with
                    # the promoted, uniquely prefixed files.
                    metadata_key = 'metadata_' + _safe_tag(analysis).lower() + '.json'
                    metadata_file = run_dir / mapping.get(metadata_key, '') if mapping.get(metadata_key) else None
                    if metadata_file and metadata_file.is_file():
                        meta = json.loads(metadata_file.read_text())
                        meta['result'] = result
                        meta['output_files'] = mapping
                        metadata_file.write_text(json.dumps(meta, indent=2, default=str))
                    summary[analysis] = {'status':'ok', 'result':result, 'file_prefix':label}
                    preview_label = f'{analysis} | {group_tag} | {run_id}'
                    self.events.put(('completed', (preview_label, str(run_dir / result['plot']), result)))
                except Exception as e:
                    summary[analysis] = {'status':'error', 'message':str(e)}
                    self.events.put(('log', f'ERROR {analysis}: {e}\n'))
                self.events.put(('progress', (idx+1)*100))
            (run_dir / f'summary_{run_id}.json').write_text(json.dumps(summary, indent=2))
            self.events.put(('done', str(run_dir)))
        threading.Thread(target=worker, daemon=True).start()

    def stop(self):
        self.cancel.set()
        self.status.set('Cancellation requested; current operation may need to finish.')

    def poll(self):
        if not self.win.winfo_exists():
            return
        for _ in range(100):
            try:
                kind, info = self.events.get_nowait()
            except queue.Empty:
                break
            if kind == 'progress':
                self.progress['value'] = info
            elif kind == 'status':
                self.status.set(info)
            elif kind == 'log':
                self.log.insert('end', info); self.log.see('end')
            elif kind == 'completed':
                label, filename, result = info
                self.completed[label] = filename
                self.results_box.configure(values=list(self.completed))
                self.results_box.set(label)
                self.show_result()
                self.log.insert('end', f'OK {label}: {result}\n'); self.log.see('end')
            elif kind == 'done':
                self.running = False
                self.run_button.configure(state='normal')
                self.stop_button.configure(state='disabled')
                self.status.set('Finished. Results: ' + info)
                self.log.insert('end', 'Saved unique run summary; earlier files retained.\n'); self.log.see('end')
        self.win.after(100, self.poll)

    def show_result(self, event=None):
        label = self.results_box.get()
        if not label or label not in self.completed:
            return
        try:
            from PIL import Image, ImageTk
            for child in self.plot_holder.winfo_children():
                child.destroy()
            img = Image.open(self.completed[label])
            img.thumbnail((600, 470))
            self._image = ImageTk.PhotoImage(img)
            tk.Label(self.plot_holder, image=self._image).pack(anchor='center', expand=True)
        except ImportError:
            ttk.Label(self.plot_holder, text='Install Pillow for in-GUI image previews.\nPNG/PDF are saved.').pack()
        except Exception as e:
            self.log.insert('end', 'Preview unavailable: ' + str(e) + '\n')
