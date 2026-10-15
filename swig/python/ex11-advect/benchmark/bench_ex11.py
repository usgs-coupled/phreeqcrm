r"""Time the ex11-advect transport loop on two PhreeqcRM builds and check they agree.

Each build is selected by putting its swig/python directory on PYTHONPATH, so any
two build trees (or installed packages) can be compared. Runs alternate in
A B B A order to spread thermal and background noise evenly over both builds.

Usage, building both sides from source (works on Windows, macOS and Linux):
    python bench_ex11.py --build-against origin/master

    This builds the git ref given (here origin/master) and the current working tree with
    pip, then compares them. It needs git, a C++ compiler and internet access.

Usage, with two builds you already have (replace the example paths with your own):
    python bench_ex11.py --base ../phreeqcrm-master/build/swig/python \
                         --branch build/swig/python \
                         --cells 1000 10000 --rounds 2
"""
import argparse
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX11 = os.path.dirname(HERE)
OUTPUTS = ["solution_total_molality_Cl", "solution_total_molality_Na",
           "solution_total_molality_K", "solution_total_molality_Ca",
           "exchange_X_species_log_molality_X-", "exchange_X_species_log_molality_NaX",
           "exchange_X_species_log_molality_KX", "exchange_X_species_log_molality_CaX2"]


def run_model(workdir, nxyz, nsteps, tag):
    """Same setup and loop as ex11-advect.ipynb, without printing or plotting."""
    import phreeqcrm
    try:
        from phreeqcrm import YAMLPhreeqcRM
    except ImportError:  # build trees ship yamlphreeqcrm.py as a loose module
        from yamlphreeqcrm import YAMLPhreeqcRM
    os.chdir(workdir)

    yrm = YAMLPhreeqcRM()
    yrm.YAMLSetGridCellCount(nxyz)
    yrm.YAMLSetErrorHandlerMode(1)
    yrm.YAMLSetComponentH2O(False)
    yrm.YAMLSetRebalanceFraction(0.5)
    yrm.YAMLSetRebalanceByCell(True)
    yrm.YAMLUseSolutionDensityVolume(False)
    yrm.YAMLSetPartitionUZSolids(False)
    yrm.YAMLSetUnitsSolution(2)
    for f in ("PPassemblage", "Exchange", "Surface", "GasPhase", "SSassemblage", "Kinetics"):
        getattr(yrm, f"YAMLSetUnits{f}")(1)
    yrm.YAMLSetTimeConversion(1.0 / 86400.0)
    yrm.YAMLSetRepresentativeVolume([1] * nxyz)
    yrm.YAMLSetDensityUser([1.0] * nxyz)
    yrm.YAMLSetPorosity([0.2] * nxyz)
    yrm.YAMLSetSaturationUser([1] * nxyz)
    yrm.YAMLLoadDatabase("phreeqc.dat")
    yrm.YAMLRunFile(True, True, True, "advect.pqi")
    yrm.YAMLRunString(True, False, True, "DELETE; -all")
    yrm.YAMLAddOutputVars("AddOutputVars", "true")
    yrm.YAMLFindComponents()
    yrm.YAMLInitialSolutions2Module([1] * nxyz)
    yrm.YAMLInitialExchanges2Module([1] * nxyz)
    yaml_file = f"ex11-{tag}-{nxyz}.yaml"
    yrm.WriteYAMLDoc(yaml_file)

    rm = phreeqcrm.BMIPhreeqcRM()
    rm.initialize(yaml_file)
    ncomps = rm.get_value_ptr("ComponentCount")[0]
    c = rm.get_value_ptr("Concentrations")
    bc = rm.InitialPhreeqc2Concentrations([0])
    bufs = [np.empty(nxyz) for _ in OUTPUTS]
    results = np.empty((nsteps, len(OUTPUTS), nxyz))

    t_update = 0.0
    t0 = time.perf_counter()
    for s in range(nsteps):
        a = time.perf_counter()
        rm.update()
        t_update += time.perf_counter() - a
        # advect one cell downstream; numpy shift is the notebook's loop, vectorized
        cc = c.reshape(ncomps, nxyz)
        cc[:, 1:] = cc[:, :-1].copy()
        cc[:, 0] = bc
        rm.set_value("Concentrations", c)
        for k, name in enumerate(OUTPUTS):
            bufs[k] = rm.get_value(name, bufs[k])
            results[s, k] = bufs[k]
    t_loop = time.perf_counter() - t0
    rm.finalize()
    np.save(f"out-{tag}-{nxyz}.npy", results)
    return dict(tag=tag, nxyz=nxyz, nsteps=nsteps, loop_s=t_loop, update_s=t_update,
                phreeqcrm_file=phreeqcrm.__file__)


def spawn(pythonpath, workdir, nxyz, nsteps, tag):
    env = dict(os.environ, PYTHONPATH=os.path.abspath(pythonpath))
    out = subprocess.run([sys.executable, os.path.abspath(__file__), "--child",
                          workdir, str(nxyz), str(nsteps), tag],
                         env=env, check=True, capture_output=True, text=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def pip_build(src, dest):
    """Build the phreeqcrm wheel from src into dest and return the folder to put on PYTHONPATH."""
    print(f"# building {src}", flush=True)
    # OpenMP off so both builds run on one thread, whatever the compiler supports
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--target", dest,
                    "--config-settings=cmake.define.PHREEQCRM_DISABLE_OPENMP=ON", src],
                   check=True)
    # the wheel installs phreeqcrm.py, _phreeqcrm and yamlphreeqcrm.py into this folder
    return os.path.join(dest, "phreeqcrm")


def build_both(ref, tmp):
    repo = subprocess.run(["git", "-C", HERE, "rev-parse", "--show-toplevel"],
                          check=True, capture_output=True, text=True).stdout.strip()
    src = os.path.join(tmp, "src-base")
    subprocess.run(["git", "-C", repo, "worktree", "add", "--detach", src, ref], check=True)
    try:
        base = pip_build(src, os.path.join(tmp, "base"))
    finally:
        subprocess.run(["git", "-C", repo, "worktree", "remove", "--force", src], check=True)
    branch = pip_build(repo, os.path.join(tmp, "branch"))
    return base, branch


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        workdir, nxyz, nsteps, tag = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        print(json.dumps(run_model(workdir, nxyz, nsteps, tag)))
        return

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-against", metavar="REF",
                    help="build this git ref and the current working tree, then compare them")
    ap.add_argument("--base", help="swig/python dir of the reference build")
    ap.add_argument("--branch", help="swig/python dir of the build under test")
    ap.add_argument("--cells", type=int, nargs="+", default=[1000, 10000])
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--rounds", type=int, default=2, help="A B B A rounds per grid size")
    ap.add_argument("--json", help="also write all timings to this file")
    a = ap.parse_args()
    if not a.build_against and not (a.base and a.branch):
        ap.error("give --build-against REF, or both --base and --branch")

    if not a.build_against:
        return compare(a)
    builddir = tempfile.mkdtemp(prefix="ex11-build-")
    try:
        a.base, a.branch = build_both(a.build_against, builddir)
        compare(a)
    finally:
        shutil.rmtree(builddir, ignore_errors=True)


def compare(a):
    workdir = tempfile.mkdtemp(prefix="ex11-bench-")
    for f in ("phreeqc.dat", "advect.pqi"):
        shutil.copy(os.path.join(EX11, f), workdir)

    print(f"# {platform.platform()}, {platform.processor() or platform.machine()}, "
          f"Python {platform.python_version()}, {a.steps} steps, {a.rounds} A-B-B-A rounds", flush=True)
    records, rows = [], []
    for n in a.cells:
        for _ in range(a.rounds):
            for tag in ("base", "branch", "branch", "base"):
                r = spawn(a.base if tag == "base" else a.branch, workdir, n, a.steps, tag)
                records.append(r)
                print(f"  {n:6d} cells  {tag:6s}  loop {r['loop_s']:8.2f} s", flush=True)
        same = np.array_equal(np.load(os.path.join(workdir, f"out-base-{n}.npy")),
                              np.load(os.path.join(workdir, f"out-branch-{n}.npy")))
        med = {t: statistics.median(r["loop_s"] for r in records if r["nxyz"] == n and r["tag"] == t)
               for t in ("base", "branch")}
        rows.append((n, med["base"], med["branch"], same))

    head = ["Cells", "Base (s)", "Branch (s)", "Speedup", "Time saved", "Outputs identical"]
    body = [[f"{n:,}", f"{b:.2f}", f"{p:.2f}", f"{b / p:.2f}x", f"{(1 - p / b) * 100:.0f}%",
             "yes" if same else "NO"] for n, b, p, same in rows]
    # padded so it lines up in a terminal and still pastes as a markdown table
    w = [max(len(r[i]) for r in [head] + body) for i in range(len(head))]
    print("\n| " + " | ".join(h.ljust(w[i]) for i, h in enumerate(head)) + " |")
    print("|" + "|".join("-" * (w[i] + 1) + ":" for i in range(5)) + "|:" + "-" * w[5] + ":|")
    for r in body:
        print("| " + " | ".join(c.rjust(w[i]) if i < 5 else c.center(w[i]) for i, c in enumerate(r)) + " |")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(records, fh, indent=1)
    shutil.rmtree(workdir)
    if not all(r[3] for r in rows):
        sys.exit("outputs differ between builds")


if __name__ == "__main__":
    main()
