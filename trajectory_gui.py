#!/usr/bin/env python3
"""Local GROMACS trajectory processor. Python standard library only."""
import json
import os
from pathlib import Path
import queue
import re
import shlex
import signal
import subprocess
import tempfile
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from datetime import datetime


def read_groups(path):
    groups = []
    for line in Path(path).read_text().splitlines():
        match = re.match(r'^\s*\[\s*(.*?)\s*\]', line)
        if match:
            groups.append(match.group(1))
    if not groups:
        raise ValueError('No index groups found.')
    return groups


def pipeline(executable, tpr, xtc, ndx, outdir, stem, center, fit, output):
    """Keep all atoms through intermediate stages so TPR/index IDs remain valid."""
    base = [executable, 'trjconv', '-s', str(tpr), '-n', str(ndx)]
    paths = [Path(outdir) / (stem + suffix + '.xtc') for suffix in
             ('.whole', '.whole.nojump', '.whole.nojump.mol', '.fit')]
    extra = [(['-pbc', 'whole'], ['System']),
             (['-pbc', 'nojump'], ['System']),
             (['-pbc', 'mol', '-center'], [center, 'System']),
             (['-fit', 'rot+trans'], [fit, output])]
    steps = []
    source = Path(xtc)
    for target, (flags, selections) in zip(paths, extra):
        steps.append((base + ['-f', str(source), '-o', str(target)] + flags,
                      '\n'.join(selections) + '\n', target))
        source = target
    return steps


def tcl_literal(value):
    # binary format H* works with older Tcl shipped in VMD 1.9.3.
    # Hex encoding prevents Tcl substitution in filenames and selections.
    return '[encoding convertfrom utf-8 [binary format H* {' + str(value).encode('utf-8').hex() + '}]]'


def vmd_script(structure, trajectory, stride=1, toolkit='', selection='all'):
    if int(stride) < 1:
        raise ValueError('Frame stride must be at least 1.')
    text = 'set mid [mol new ' + tcl_literal(structure) + ' waitfor all]\n'
    text += 'mol addfile ' + tcl_literal(trajectory) + ' type xtc step ' + str(int(stride)) + ' waitfor all molid $mid\n'
    text += 'mol delrep 0 $mid\nmol representation VDW 0.35 12\nmol color ResName\n'
    text += 'mol selection ' + tcl_literal(selection) + '\nmol addrep $mid\n'
    text += 'display projection Orthographic\ncolor Display Background white\naxes location Off\n'
    if toolkit:
        text += 'source ' + tcl_literal(toolkit) + '\n'
    text += 'animate goto 0\ndisplay resetview\n'
    return text


class App:
    def __init__(self, root):
        self.root = root
        root.title('GROMACS trajectory processor')
        root.geometry('1050x950')
        root.minsize(800, 650)
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.process = None
        self.busy = False
        self.groups = []
        self.loaded = None
        self.generated_index = None
        self.temp = tempfile.TemporaryDirectory(prefix='gmx-gui-')
        self.v = {k: tk.StringVar() for k in ('gmx', 'tpr', 'xtc', 'ndx', 'out', 'stem', 'vmd', 'structure', 'view_xtc', 'toolkit', 'stride', 'selection')}
        self.v['gmx'].set('gmx')
        self.v['vmd'].set('vmd')
        self.v['stride'].set('10')
        self.v['selection'].set('all')
        self.last_outputs = {}
        style = ttk.Style()
        style.configure('Title.TLabel', font=('Sans', 20, 'bold'))
        body = ttk.Frame(root, padding=20)
        body.pack(fill='both', expand=True)
        ttk.Label(body, text='Trajectory processing', style='Title.TLabel').pack(anchor='w')
        ttk.Label(body, text='Whole molecules → no jumps → molecule centering → rotational / translational fit').pack(anchor='w', pady=(4, 18))
        files = ttk.Frame(body)
        files.pack(fill='x')
        files.columnconfigure(1, weight=1)
        labels = [('gmx', 'GROMACS executable'), ('tpr', 'Reference TPR'), ('xtc', 'Input trajectory XTC'),
                  ('ndx', 'Index NDX (optional)'), ('out', 'Output parent folder'), ('stem', 'Output filename prefix')]
        for row, (key, label) in enumerate(labels):
            ttk.Label(files, text=label).grid(row=row, column=0, sticky='w', padx=(0, 15), pady=5)
            ttk.Entry(files, textvariable=self.v[key]).grid(row=row, column=1, sticky='ew', pady=5)
            if key != 'stem':
                ttk.Button(files, text='Browse', command=lambda k=key: self.browse(k)).grid(row=row, column=2, padx=(8, 0))
        self.load_button = ttk.Button(body, text='Load groups from TPR / index', command=self.load_groups)
        self.load_button.pack(anchor='w', pady=12)
        groups = ttk.LabelFrame(body, text='Group selections', padding=10)
        groups.pack(fill='x')
        self.boxes = []
        for col, label in enumerate(('Centering group', 'Fit group', 'Final output group')):
            groups.columnconfigure(col, weight=1)
            ttk.Label(groups, text=label).grid(row=0, column=col, sticky='w')
            box = ttk.Combobox(groups, state='readonly')
            box.grid(row=1, column=col, sticky='ew', padx=(0, 10), pady=6)
            self.boxes.append(box)
        ttk.Label(body, text='Intermediate files retain System. Final output can contain your chosen group.\nChoose a stable core for fitting if you want to retain relative nucleosome motion.').pack(anchor='w', pady=10)
        actions = ttk.Frame(body)
        actions.pack(fill='x')
        self.run_button = ttk.Button(actions, text='Run all four steps', command=self.run)
        self.run_button.pack(side='left')
        self.stop_button = ttk.Button(actions, text='Stop', command=self.stop, state='disabled')
        self.stop_button.pack(side='left', padx=10)
        self.status = tk.StringVar(value='Select files, then load groups.')
        ttk.Label(actions, textvariable=self.status).pack(side='left', padx=10)
        self.progress = ttk.Progressbar(body, maximum=4)
        self.progress.pack(fill='x', pady=12)
        viewer = ttk.LabelFrame(body, text='VMD visualization (opens a separate VMD window)', padding=10)
        viewer.pack(fill='x', pady=(0, 10))
        viewer.columnconfigure(1, weight=1)
        for row, (key, label) in enumerate((('vmd', 'VMD executable'), ('structure', 'Matching GRO / PDB / PSF'),
                                            ('view_xtc', 'Trajectory to visualize'), ('toolkit', 'SIRAH toolkit TCL (optional)'))):
            ttk.Label(viewer, text=label).grid(row=row, column=0, sticky='w', padx=(0, 10), pady=3)
            ttk.Entry(viewer, textvariable=self.v[key]).grid(row=row, column=1, sticky='ew')
            ttk.Button(viewer, text='Browse', command=lambda k=key: self.browse(k)).grid(row=row, column=2, padx=5)
        options = ttk.Frame(viewer)
        options.grid(row=4, column=0, columnspan=3, sticky='ew', pady=6)
        ttk.Label(options, text='Load every Nth frame:').pack(side='left')
        ttk.Entry(options, width=5, textvariable=self.v['stride']).pack(side='left', padx=5)
        ttk.Label(options, text='Atom selection:').pack(side='left', padx=5)
        ttk.Entry(options, width=23, textvariable=self.v['selection']).pack(side='left')
        ttk.Button(options, text='Open in VMD', command=self.open_vmd).pack(side='left', padx=10)
        ttk.Button(options, text='Use completed final output', command=self.use_final).pack(side='left')
        analysis = ttk.LabelFrame(body, text='Trajectory analysis (MDAnalysis)', padding=10)
        analysis.pack(fill='x', pady=(0, 10))
        ttk.Label(analysis, text='Residue–residue contact occupancy maps (atomistic or coarse-grained).').pack(side='left', padx=(0, 12))
        ttk.Button(analysis, text='Contact map…', command=self.open_contact_map).pack(side='left')
        ttk.Button(analysis, text='Analysis workbench…', command=self.open_workbench).pack(side='left', padx=12)
        self.log = tk.Text(body, height=10, wrap='word', font=('Monospace', 10))
        self.log.pack(fill='both', expand=True)
        root.after(100, self.poll)
        root.protocol('WM_DELETE_WINDOW', self.close)

    def open_workbench(self):
        try:
            from analysis_workbench import AnalysisWorkbench
        except ImportError as exc:
            messagebox.showerror('Analysis module missing', str(exc))
            return
        structure = self.v['structure'].get().strip() or self.v['tpr'].get().strip()
        trajectory = self.v['view_xtc'].get().strip() or self.v['xtc'].get().strip()
        if self.last_outputs:
            structure = self.last_outputs['structure']
            trajectory = self.last_outputs['trajectory']
        AnalysisWorkbench(self.root, topology=structure, trajectory=trajectory,
                          tpr=self.v['tpr'].get().strip(), ndx=self.v['ndx'].get().strip(),
                          gmx=self.v['gmx'].get().strip())

    def open_contact_map(self):
        try:
            from contact_map_gui import ContactMapWindow
        except ImportError as exc:
            messagebox.showerror('Contact-map module missing', str(exc))
            return
        structure = self.v['structure'].get().strip()
        trajectory = self.v['view_xtc'].get().strip()
        if self.last_outputs:
            structure = self.last_outputs['structure']
            trajectory = self.last_outputs['trajectory']
        if not structure and self.v['tpr'].get().strip():
            structure = self.v['tpr'].get().strip()
        if not trajectory and self.v['xtc'].get().strip():
            trajectory = self.v['xtc'].get().strip()
        ContactMapWindow(self.root, structure=structure, trajectory=trajectory)

    def browse(self, key):
        if self.busy:
            return
        if key == 'out':
            value = filedialog.askdirectory()
        else:
            types = {'tpr': [('TPR', '*.tpr')], 'xtc': [('XTC', '*.xtc')], 'ndx': [('Index', '*.ndx')], 'structure': [('Structure', '*.gro *.pdb *.psf')], 'view_xtc': [('XTC', '*.xtc')], 'toolkit': [('Tcl script', '*.tcl')]}
            value = filedialog.askopenfilename(filetypes=types.get(key, [('Executable', '*')]))
        if value:
            self.v[key].set(value)
            if key == 'xtc':
                self.v['stem'].set(Path(value).stem)
                if not self.v['out'].get():
                    self.v['out'].set(str(Path(value).parent))

    def use_final(self):
        if not self.last_outputs:
            messagebox.showinfo('No completed output', 'Complete a processing run first, or browse existing structure and trajectory files.')
            return
        self.v['structure'].set(self.last_outputs['structure'])
        self.v['view_xtc'].set(self.last_outputs['trajectory'])

    def open_vmd(self):
        cfg = self.snapshot()
        try:
            for key in ('structure', 'view_xtc'):
                if not Path(cfg[key]).is_file():
                    raise ValueError('Select an existing matching structure and XTC trajectory.')
            if cfg['toolkit'] and not Path(cfg['toolkit']).is_file():
                raise ValueError('The SIRAH toolkit file does not exist.')
            script = vmd_script(str(Path(cfg['structure']).resolve()), str(Path(cfg['view_xtc']).resolve()),
                                int(cfg['stride']), str(Path(cfg['toolkit']).resolve()) if cfg['toolkit'] else '', cfg['selection'])
            folder = Path(tempfile.mkdtemp(prefix='vmd_view_', dir=Path(cfg['view_xtc']).resolve().parent))
            script_path = folder / 'view.tcl'
            script_path.write_text(script)
            # Keep launch scripts/logs beside the trajectory; an active VMD session survives GUI closure.
            args = [cfg['vmd'], '-e', str(script_path)]
            log = (folder / 'vmd.log').open('w')
            try:
                proc = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT,
                                        start_new_session=(os.name == 'posix'))
            finally:
                log.close()
            self.events.put(('log', '\nVMD launched: ' + shlex.join(args) + '\nVMD log: ' + str(folder / 'vmd.log') + '\n'))
            def monitor():
                code = proc.wait()
                if code:
                    self.events.put(('vmd_error', 'VMD exited with code ' + str(code) + '. See ' + str(folder / 'vmd.log')))
            threading.Thread(target=monitor, daemon=True).start()
        except Exception as exc:
            messagebox.showerror('VMD could not launch', str(exc))

    def snapshot(self):
        return {k: v.get().strip() for k, v in self.v.items()}

    def set_busy(self, value):
        self.busy = value
        self.run_button.configure(state='disabled' if value else 'normal')
        self.load_button.configure(state='disabled' if value else 'normal')
        self.stop_button.configure(state='normal' if value else 'disabled')

    def execute(self, args, stdin, logfile=None):
        if self.cancel.is_set():
            raise RuntimeError('Stopped by user.')
        self.events.put(('log', '$ ' + shlex.join(args) + '\n'))
        self.events.put(('log', 'Selections: ' + repr(stdin) + '\n'))
        env = dict(os.environ, GMX_MAXBACKUP='-1')
        process = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, errors='replace',
                                   start_new_session=(os.name == 'posix'), env=env)
        self.process = process
        if self.cancel.is_set():
            self.terminate(process)
        try:
            process.stdin.write(stdin)
            process.stdin.close()
            for line in process.stdout:
                self.events.put(('log', line))
                if logfile:
                    logfile.write(line)
                    logfile.flush()
            code = process.wait()
            if self.cancel.is_set():
                raise RuntimeError('Stopped by user. Partial outputs remain in the run folder.')
            if code:
                raise RuntimeError(f'GROMACS exited with code {code}. See log above.')
        finally:
            if process.poll() is None:
                self.terminate(process)
                process.wait()
            self.process = None

    def load_groups(self):
        cfg = self.snapshot()
        if not Path(cfg['tpr']).is_file():
            messagebox.showerror('Missing TPR', 'Select a valid TPR file.')
            return
        if cfg['ndx'] and not Path(cfg['ndx']).is_file():
            messagebox.showerror('Missing index', 'The selected index file does not exist.')
            return
        self.cancel.clear()
        self.set_busy(True)
        self.status.set('Loading index groups…')
        def worker():
            try:
                target = Path(self.temp.name) / ('groups_' + datetime.now().strftime('%H%M%S%f') + '.ndx')
                args = [cfg['gmx'], 'make_ndx', '-f', cfg['tpr'], '-o', str(target)]
                if cfg['ndx']:
                    args += ['-n', cfg['ndx']]
                self.execute(args, 'q\n')
                names = read_groups(target)
                if 'System' not in names or names.count('System') != 1:
                    raise ValueError('The index must contain one complete System group.')
                if len(names) != len(set(names)):
                    raise ValueError('Index group names must be unique. Rename duplicate groups in your index.')
                self.events.put(('groups', (names, cfg, str(target))))
            except Exception as exc:
                self.events.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    def run(self):
        cfg = self.snapshot()
        if not self.loaded or any(cfg[k] != self.loaded[k] for k in ('gmx', 'tpr', 'ndx')):
            messagebox.showerror('Load groups', 'Load groups again after changing the executable, TPR or index.')
            return
        if not Path(cfg['xtc']).is_file() or not Path(cfg['tpr']).is_file():
            messagebox.showerror('Missing input', 'Select valid XTC and TPR files.')
            return
        if not cfg['out'] or not cfg['stem'] or Path(cfg['stem']).name != cfg['stem'] or cfg['stem'] in ('.', '..'):
            messagebox.showerror('Output', 'Choose an output folder and a filename prefix without directory separators.')
            return
        selections = [box.get() for box in self.boxes]
        if not all(x in self.groups for x in selections):
            messagebox.showerror('Groups', 'Select centering, fitting and final output groups.')
            return
        self.cancel.clear()
        self.set_busy(True)
        self.progress['value'] = 0
        def worker():
            try:
                parent = Path(cfg['out']).expanduser().resolve()
                parent.mkdir(parents=True, exist_ok=True)
                directory = Path(tempfile.mkdtemp(prefix=cfg['stem'] + '_processed_', dir=parent))
                ndx = directory / 'processing_groups.ndx'
                ndx.write_bytes(Path(self.generated_index).read_bytes())
                steps = pipeline(cfg['gmx'], Path(cfg['tpr']).resolve(), Path(cfg['xtc']).resolve(), ndx,
                                 directory, cfg['stem'], *selections)
                manifest = dict(cfg, selected_groups=selections, steps=[{'command': a, 'stdin': b} for a,b,_ in steps])
                (directory / 'run.json').write_text(json.dumps(manifest, indent=2))
                with (directory / 'processing.log').open('w') as log:
                    for i, (args, stdin, target) in enumerate(steps, 1):
                        self.events.put(('status', f'Step {i}/4: {target.name}'))
                        log.write('\n$ ' + shlex.join(args) + '\nSelections: ' + repr(stdin) + '\n')
                        self.execute(args, stdin, log)
                        if not target.is_file() or not target.stat().st_size:
                            raise RuntimeError('GROMACS did not produce a nonempty output: ' + str(target))
                        self.events.put(('progress', i))
                    # Use the full pre-fit input and identical selections, ensuring subset atom order matches the final XTC.
                    args, stdin, final = steps[-1]
                    structure = final.with_suffix('.gro')
                    structure_args = list(args)
                    structure_args[structure_args.index('-o') + 1] = str(structure)
                    structure_args += ['-dump', '0']
                    self.events.put(('status', 'Preparing matching VMD structure…'))
                    log.write('\n$ ' + shlex.join(structure_args) + '\n')
                    self.execute(structure_args, stdin, log)
                    if not structure.is_file() or not structure.stat().st_size:
                        raise RuntimeError('Could not create matching VMD structure.')
                self.events.put(('done', {'trajectory': str(final), 'structure': str(structure)}))
            except Exception as exc:
                self.events.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    @staticmethod
    def terminate(process):
        try:
            if os.name == 'posix':
                os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
        except ProcessLookupError:
            pass

    def stop(self):
        self.cancel.set()
        if self.process:
            self.terminate(self.process)
        self.status.set('Stopping…')

    def poll(self):
        # Bound each drain so a verbose subprocess cannot starve the interface.
        for _ in range(250):
            try:
                kind, value = self.events.get_nowait()
            except queue.Empty:
                break
            if kind == 'log':
                self.log.insert('end', value)
                if int(self.log.index('end-1c').split('.')[0]) > 12000:
                    self.log.delete('1.0', '2000.0')
                self.log.see('end')
            elif kind == 'groups':
                self.groups, self.loaded, self.generated_index = value
                for box in self.boxes:
                    box['values'] = self.groups
                    box.set('')
                self.boxes[-1].set('System')
                self.status.set('Groups loaded. Choose centering and fitting groups.')
                self.set_busy(False)
            elif kind == 'progress':
                self.progress['value'] = value
            elif kind == 'status':
                self.status.set(value)
            elif kind == 'error':
                self.set_busy(False)
                self.status.set('Stopped / failed. See log.')
                messagebox.showerror('Processing stopped', value)
            elif kind == 'vmd_error':
                messagebox.showerror('VMD exited', value)
            elif kind == 'done':
                self.last_outputs = value
                self.use_final()
                self.set_busy(False)
                self.status.set('Completed all four steps.')
                messagebox.showinfo('Completed', 'Final fitted trajectory:\n' + value['trajectory'] + '\n\nMatching GRO structure ready. Click Open in VMD to visualize.')
        self.root.after(100, self.poll)

    def close(self):
        if self.busy:
            messagebox.showinfo('Processing active', 'Stop processing before closing the window.')
            return
        self.temp.cleanup()
        self.root.destroy()


if __name__ == '__main__':
    App(tk.Tk())
    tk.mainloop()
