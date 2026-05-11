import sys
from pathlib import Path
sys.path.append(str(Path("/home/diegoa/dev/pysoc")))

from pysoc.io.SOC import Calculator

try:
    calc = Calculator("/home/diegoa/dev/pysoc/sample_orca/tddft_pysoc.out", num_singlets=5, num_triplets=5, QM_program="ORCA")
    table = calc.calculate()
    print(table)
except Exception as e:
    import subprocess
    if isinstance(e, subprocess.CalledProcessError):
        print("STDOUT:")
        print(e.stdout)
    else:
        raise
