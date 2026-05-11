from pathlib import Path
import re
import cclib
import periodictable

from pysoc.io.molsoc import Molsoc

class Orca_parser(Molsoc):
    """
    Class for parsing the required output from ORCA.
    """
    
    # Molsoc ALWAYS uses Cartesian basis functions internally!
    # Therefore, we MUST use Cartesian shell sizes for ao_basis.dat.
    SHELLS = {"S": 1, "P": 3, "D": 6, "F": 10, "G": 15, "H": 21, "I": 28}

    def __init__(self, out_file_name, requested_singlets, requested_triplets):
        super().__init__(requested_singlets, requested_triplets)
        self.out_file_name = Path(out_file_name)
        self.inp_file_name = self.out_file_name
        self.num_frozen_orbitals = 0 
        self.basis_set = []

    @classmethod
    def from_output_files(cls, out_file_name, **kwargs):
        out_file_name = Path(out_file_name)
        return cls(out_file_name=out_file_name, **kwargs)

    def parse(self):
        self._parse_cclib()
        self._parse_mkl_basis()
        self._parse_text()

    def _parse_cclib(self):
        with open(self.out_file_name, 'r') as out_file:
            ccdata = cclib.io.ccread(out_file)
            
            self.num_atoms = ccdata.natom
            self.geometry = [[str(periodictable.elements[proton_num]), coord[0], coord[1], coord[2]] 
                             for proton_num, coord in zip(ccdata.atomnos, ccdata.atomcoords[-1])]
            
            self.num_orbitals = ccdata.nmo
            self.num_occupied_orbitals = ccdata.homos[0] + 1
            self.num_virtual_orbitals = self.num_orbitals - self.num_occupied_orbitals
            
            if hasattr(ccdata, 'etsyms'):
                for es_index, es_symm in enumerate(ccdata.etsyms):
                    if "Singlet" in es_symm:
                        es_list = self.singlet_states
                    elif "Triplet" in es_symm:
                        es_list = self.triplet_states
                    else:
                        continue
                    es_list.append([es_index+1, round(self.wavenumbers_to_energy(ccdata.etenergies[es_index]), 4)])

    def _parse_mkl_basis(self):
        # Determine paths
        mkl_file = self.out_file_name.with_suffix('.mkl')
        
        if not mkl_file.exists():
            import subprocess
            print(f"Warning: {mkl_file.name} not found. Attempting to generate it with orca_2mkl.")
            try:
                subprocess.run(['orca_2mkl', self.out_file_name.stem, '-mkl'], cwd=self.out_file_name.parent, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                print(f"Warning: orca_2mkl failed or not found. Make sure orca_2mkl is in your PATH and {self.out_file_name.with_suffix('.gbw').name} exists.")
        
        if not mkl_file.exists():
            raise FileNotFoundError(f"Required {mkl_file.name} not found and could not be generated.")
            
        lines = mkl_file.read_text().splitlines()
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
                    
        for i, atom in enumerate(self.geometry):
            element = atom[0]
            self.basis_set.append(f"{element}  0\n")
            
            block = basis_blocks[i]
            j = 0
            while j < len(block):
                header = block[j].split()
                if len(header) >= 3 and header[1] in self.SHELLS:
                    shell_type = header[1]
                    scale = header[2]
                    j += 1
                    primitives = []
                    while j < len(block) and len(block[j].split()) == 2:
                        primitives.append(block[j] + '\n')
                        j += 1
                    self.basis_set.append(f"{shell_type}   {len(primitives)} {scale}\n")
                    self.basis_set.extend(primitives)
                    self.ao_basis.append(self.SHELLS[shell_type])
                else:
                    j += 1
            self.basis_set.append("****\n")

        # Pop the last ****
        if len(self.basis_set) > 0 and self.basis_set[-1] == "****\n":
            self.basis_set.pop()

    def _parse_text(self):
        with open(self.out_file_name, 'r') as f:
            lines = f.readlines()
        
        overlap_matrix = []
        in_overlap = False
        
        mo_energies = [0.0] * self.num_orbitals
        mo_coeffs = [[0.0] * self.num_orbitals for _ in range(self.num_orbitals)]
        in_mos = False
        is_energy_line = False
        ao_idx = 0
        aonames = []
        cols = []
        
        in_singlets = False
        in_triplets = False
        
        ci_dict_singlets = {i: [0.0]*self.ndim for i in range(1, len(self.singlet_states)+1)}
        ci_dict_triplets = {i: [0.0]*self.ndim for i in range(1, len(self.triplet_states)+1)}
        current_state = None
        
        i = 0
        while i < len(lines):
            line = lines[i]
            
            # 2. OVERLAP MATRIX
            if "OVERLAP MATRIX" in line:
                in_overlap = True
                i += 2
                overlap_matrix = [[0.0]*self.num_orbitals for _ in range(self.num_orbitals)]
                continue
            if in_overlap:
                if "Time for" in line or "ORBITAL ENERGIES" in line or "MOLECULAR ORBITALS" in line or "DIPOLE" in line:
                    in_overlap = False
                    continue
                if line.strip() == "" or "----------------" in line:
                    continue
                elif len(line.split()) > 0 and all(x.isdigit() for x in line.split()):
                    cols = [int(x) for x in line.split()]
                elif len(line.split()) > 1 and line.split()[0].isdigit():
                    parts = line.split()
                    row_idx = int(parts[0])
                    for col_idx_offset, val in enumerate(parts[1:]):
                        overlap_matrix[row_idx][cols[col_idx_offset]] = float(val)
                            
            # 3. MOLECULAR ORBITALS
            if "MOLECULAR ORBITALS" in line:
                in_mos = True
                i += 2
                continue
            if in_mos:
                if len(line.strip()) == 0 and len(lines[i+1].strip()) == 0:
                    in_mos = False
                elif len(line.split()) > 0 and all(x.isdigit() for x in line.split()):
                    cols = [int(x) for x in line.split()]
                    ao_idx = 0
                    is_energy_line = True
                elif "--------" in line:
                    pass
                elif len(line.split()) == len(cols): 
                    # Energies or occupancies. Energies are usually < 0 or > 2, or the first row.
                    parts = line.split()
                    try:
                        if is_energy_line:
                            for c, e in zip(cols, parts):
                                mo_energies[c] = float(e)
                            is_energy_line = False
                    except ValueError:
                        pass
                else: 
                    parts = line.split()
                    if len(parts) >= len(cols) + 2: 
                        coeffs = parts[-len(cols):]
                        for c, val in zip(cols, coeffs):
                            mo_coeffs[c][ao_idx] = float(val)
                        if cols[0] == 0:
                            aonames.append(parts[1])
                        ao_idx += 1
                        
            # 4. CI COEFFICIENTS
            if "TD-DFT EXCITED STATES (SINGLETS)" in line:
                in_singlets = True
                in_triplets = False
            elif "TD-DFT EXCITED STATES (TRIPLETS)" in line:
                in_singlets = False
                in_triplets = True
            elif "Entering triplet calculation" in line or "RPA-DIAGONALIZATION" in line:
                in_singlets = False
            
            if in_singlets or in_triplets:
                if line.startswith("STATE "):
                    current_state = int(line.split()[1][:-1])
                elif "->" in line:
                    parts = line.split()
                    occ = int(parts[0][:-1]) 
                    virt = int(parts[2][:-1]) 
                    weight = float(parts[4])
                    
                    occ_idx = occ
                    virt_idx = virt - self.num_occupied_orbitals
                    
                    if 0 <= occ_idx < self.num_occupied_orbitals and 0 <= virt_idx < self.num_virtual_orbitals:
                        transition_idx = occ_idx * self.num_virtual_orbitals + virt_idx
                        if in_singlets and current_state in ci_dict_singlets:
                            ci_dict_singlets[current_state][transition_idx] = weight
                        elif in_triplets and current_state in ci_dict_triplets:
                            ci_dict_triplets[current_state][transition_idx] = weight
            
            i += 1
            
        # Post-process parsed data into Molsoc properties
            
        # Overlaps (1D list, full matrix since we hijack tddftb branch in Fortran for Spherical harmonics)
        self.AO_overlaps = []
        for i in range(self.num_orbitals):
            for j in range(self.num_orbitals):
                # Orca overlap matrix is symmetric, but we only read the lower half.
                # So we can mirror it to get the full matrix.
                if j <= i:
                    self.AO_overlaps.append(overlap_matrix[i][j])
                else:
                    self.AO_overlaps.append(overlap_matrix[j][i])
                
        # MOs (1D list of Alpha coeffs, energies)
        self.MO_energies = mo_energies
        
        # Permute D and F shells to match DFTB+ Spherical Harmonic ordering
        # ORCA D: 0(z2), 1(xz), 2(yz), 3(x2y2), 4(xy)
        # DFTB D: 0(xy), 1(yz), 2(z2), 3(xz), 4(x2y2)
        # Mapping ORCA -> DFTB+: ORCA[4, 2, 0, 1, 3]
        d_reorder = [4, 2, 0, 1, 3]
        # ORCA F: 0(z3), 1(xz2), 2(yz2), 3(zx2y2), 4(xyz), 5(xx23y2), 6(y3x2y2)
        # DFTB F: 0(y3x2y2), 1(xyz), 2(yz2), 3(z3), 4(xz2), 5(zx2y2), 6(xx23y2)
        # Mapping ORCA -> DFTB+: ORCA[6, 4, 2, 0, 1, 3, 5]
        f_reorder = [6, 4, 2, 0, 1, 3, 5]
        
        ao_idx = 0
        while ao_idx < len(aonames):
            name = aonames[ao_idx].lower()
            if 'dz2' in name:
                # Reorder the 5 D functions
                for mo_idx in range(self.num_orbitals):
                    d_coeffs = [mo_coeffs[mo_idx][ao_idx + k] for k in range(5)]
                    for k in range(5):
                        mo_coeffs[mo_idx][ao_idx + k] = d_coeffs[d_reorder[k]]
                ao_idx += 5
            elif 'fz3' in name:
                # Reorder the 7 F functions
                for mo_idx in range(self.num_orbitals):
                    f_coeffs = [mo_coeffs[mo_idx][ao_idx + k] for k in range(7)]
                    for k in range(7):
                        mo_coeffs[mo_idx][ao_idx + k] = f_coeffs[f_reorder[k]]
                ao_idx += 7
            else:
                ao_idx += 1
                
        self.MOA_coefficients = []
        for i in range(self.num_orbitals):
            self.MOA_coefficients.extend(mo_coeffs[i])
            
        # CI Coeffs
        # soc_td with 'gauss_tddft' expects 2 * ndim coefficients per state for X+Y, 
        # followed by 2 * ndim coefficients per state for X-Y.
        ci_xpy = []
        ci_xmy = []
        
        for i in self.requested_singlets:
            if i in ci_dict_singlets:
                for coeff in ci_dict_singlets[i]:
                    # Assume Alpha and Beta components for restricted singlet (coeff / sqrt(2) or just coeff)
                    ci_xpy.extend([coeff, coeff])
                    ci_xmy.extend([coeff, coeff])

        for i in self.requested_triplets:
            if i in ci_dict_triplets:
                for coeff in ci_dict_triplets[i]:
                    # Assume Alpha and Beta components for restricted triplet
                    ci_xpy.extend([coeff, -coeff])
                    ci_xmy.extend([coeff, -coeff])

        self.CI_coefficients = ci_xpy + ci_xmy

    @property
    def num_transitions(self):
        """
        The number of transitions (AO matrix size) for molsoc integrals.
        Since molsoc always uses Cartesian basis sets, this must use ao_ncart.
        """
        return self.ao_ncart[0] ** 2

    @property
    def ao_ncart(self):
        return [self.ao_basis_sum, self.num_occupied_orbitals, self.num_virtual_orbitals]
