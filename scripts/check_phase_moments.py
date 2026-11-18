from pathlib import Path
import subprocess,shutil,json
ROOT=Path(__file__).resolve().parents[1]
case=json.loads((ROOT/'examples/catalog.json').read_text())['07_immersion_3d']['directory']
compiler=shutil.which('gfortran');assert compiler,'Install GNU Fortran for the phase-moment kernel checks'
work=ROOT/'build/phase-moment-checks';work.mkdir(parents=True,exist_ok=True)
command=[compiler,'-O2','-fno-reciprocal-math','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(ROOT/case/'src/phase_moment_core.f90'),str(ROOT/'tests/phase_moment_checks.f90'),'-o',str(work/'check')]
subprocess.run(command,cwd=work,check=True,timeout=60)
subprocess.run([str(work/'check')],cwd=work,check=True,timeout=30)
