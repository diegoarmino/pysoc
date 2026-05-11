import sys
from pathlib import Path
import tempfile
import subprocess
from pysoc.io.SOC import Calculator

calc = Calculator('/home/diegoa/dev/pysoc/sample_orca/tddft_pysoc.out', num_singlets=5, num_triplets=5, QM_program='ORCA')
calc.molsoc.parse()

# Now parse the .mkl file for basis set
mkl_path = '/home/diegoa/dev/pysoc/sample_orca/tddft_pysoc.mkl'

lines = open(mkl_path).readlines()
in_basis = False
basis_blocks = []
current_block = []

for line in lines:
    if line.startswith('$BASIS'):
        in_basis = True
        continue
    if in_basis:
        if line.startswith('$') and not line.startswith('$$'):
            basis_blocks.append(current_block)
            break
        if line.startswith('$$'):
            basis_blocks.append(current_block)
            current_block = []
        else:
            current_block.append(line)

# Now translate basis_blocks to Gaussian gfinput
SHELLS = {"S": 1, "P": 3, "D": 5, "F": 7, "G": 9, "H": 11, "I": 13}
basis_set = []

for i, atom in enumerate(calc.molsoc.geometry):
    element = atom[0]
    basis_set.append(f"{element}  0\n")
    
    block = basis_blocks[i]
    j = 0
    while j < len(block):
        header = block[j].split()
        if len(header) >= 3 and header[1] in SHELLS:
            shell_type = header[1]
            scale = header[2]
            j += 1
            primitives = []
            while j < len(block) and len(block[j].split()) == 2:
                primitives.append(block[j])
                j += 1
            basis_set.append(f"{shell_type}   {len(primitives)} {scale}\n")
            basis_set.extend(primitives)
        else:
            j += 1
    basis_set.append("****\n")

# Pop the last ****
if len(basis_set) > 0 and basis_set[-1] == "****\n":
    basis_set.pop()

calc.molsoc.basis_set = basis_set

# Recalculate ao_basis based on the new basis_set
ao_basis = []
for line in basis_set:
    parts = line.split()
    if len(parts) > 0 and parts[0] in SHELLS:
        ao_basis.append(SHELLS[parts[0]])

calc.molsoc.ao_basis = ao_basis
# Don't set ao_basis_sum, it's a property calculated from ao_basis
ao_basis_sum_val = sum(ao_basis)

# Now test molsoc
temp_dir = tempfile.mkdtemp()
calc.molsoc.output = Path(temp_dir)
calc.molsoc.write_molsoc_input({'d_basis_ncart': calc.molsoc.ao_ncart})

print("MOLSOC.INP CONTENT:\\n", Path(temp_dir, "molsoc.inp").read_text())
subprocess.run([calc.molsoc.MOLSOC_PATH, 'molsoc.inp'], cwd=temp_dir)
overlap = Path(temp_dir, 'molsoc_overlap.dat').read_text()

# Check diagonals
diagonals = []
elements = [float(x) for x in overlap.split() if "AO_overlap" not in x]
idx = 0
for k in range(1, ao_basis_sum_val+1):
    idx += k
    if idx - 1 < len(elements):
        diagonals.append(elements[idx - 1])

print("Min diagonal:", min(diagonals))
print("Max diagonal:", max(diagonals))

import shutil
shutil.rmtree(temp_dir)
