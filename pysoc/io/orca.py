from pathlib import Path
import json
import re
import math
import subprocess
import cclib
import periodictable
import warnings

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
        self.active_occupied_orbitals = []
        self.active_virtual_orbitals = []
        self._active_occ_position = {}
        self._active_virt_position = {}
        self._sparse_ci_threshold = None
        self._td_root_counts = []

    @classmethod
    def from_output_files(cls, out_file_name, **kwargs):
        out_file_name = Path(out_file_name)
        return cls(out_file_name=out_file_name, **kwargs)

    def parse(self):
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

    def _parse_mkl_lines(self):
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
        return mkl_file.read_text().splitlines()
    
    def _parse_mkl_coord(self, lines):
        # lines is a list of strings (the whole MKL file)
        in_coord = False
        for line in lines:
            if line.startswith('$COORD'):
                in_coord = True
                continue
            if in_coord:
                if line.startswith('$END'):
                    break
                parts = line.split()
                if len(parts) >= 4:
                    atnum = int(parts[0])
                    x, y, z = map(float, parts[1:4])
                    element = str(periodictable.elements[atnum])
                    self.geometry.append([element, x, y, z])
        self.num_atoms = len(self.geometry)

    def _parse_mkl_occ(self, lines):
        in_occ = False
        occs = []
        for line in lines:
            if line.startswith('$OCC_ALPHA'):
                in_occ = True
                continue
            if in_occ:
                if line.startswith('$END'):
                    break
                # extract all floats from the line
                for token in line.split():
                    try:
                        occs.append(float(token))
                    except ValueError:
                        pass
        self.num_orbitals = len(occs)
        occupied = [idx for idx, occ in enumerate(occs) if occ > 0.5]
        virtual = [idx for idx, occ in enumerate(occs) if occ <= 0.5]
        self._set_active_orbitals(occupied, virtual)

    def _set_active_orbitals(self, occupied, virtual):
        self.active_occupied_orbitals = list(occupied)
        self.active_virtual_orbitals = list(virtual)
        self._active_occ_position = {
            orbital: idx for idx, orbital in enumerate(self.active_occupied_orbitals)
        }
        self._active_virt_position = {
            orbital: idx for idx, orbital in enumerate(self.active_virtual_orbitals)
        }
        self.num_occupied_orbitals = len(self.active_occupied_orbitals)
        self.num_virtual_orbitals = len(self.active_virtual_orbitals)
        self.num_frozen_orbitals = self.active_occupied_orbitals[0] if self.active_occupied_orbitals else 0

    def _parse_orbital_range(self, line):
        match = re.search(
            r"Operator\s+\d+:\s+Orbitals\s+(-?\d+)\s*\.\.\.\s*(-?\d+)\s+to\s+(-?\d+)\s*\.\.\.\s*(-?\d+)",
            line,
        )
        if match is None:
            return

        occ_start, occ_end, virt_start, virt_end = [int(value) for value in match.groups()]
        if min(occ_start, occ_end, virt_start, virt_end) < 0:
            return

        self._set_active_orbitals(
            range(occ_start, occ_end + 1),
            range(virt_start, virt_end + 1),
        )

    def _parse_root_count(self, line):
        match = re.search(r"Number of roots to be determined\s+\.\.\.\s+(\d+)", line)
        if match is not None:
            self._td_root_counts.append(int(match.group(1)))

    def _warn_if_selected_edge_roots(self):
        selected_edges = []
        root_sets = (
            ("singlet", self.requested_singlets, self._td_root_counts[0] if len(self._td_root_counts) > 0 else None),
            ("triplet", self.requested_triplets, self._td_root_counts[1] if len(self._td_root_counts) > 1 else None),
        )

        for multiplicity, requested_roots, available_roots in root_sets:
            if available_roots is not None and requested_roots and max(requested_roots) >= available_roots:
                selected_edges.append(f"{multiplicity} root {available_roots}")

        if selected_edges:
            warnings.warn(
                "The selected ORCA TD roots include the highest root printed "
                f"({', '.join(selected_edges)}). Davidson edge roots can be "
                "unstable or miss nearby states; request extra ORCA roots and "
                "select the lower subset with -s/-t for production comparisons.",
                RuntimeWarning,
            )

    @staticmethod
    def _permute_mo_ao_block(mo_coeffs, start, order):
        for mo_idx in range(len(mo_coeffs)):
            block = [mo_coeffs[mo_idx][start + offset] for offset in range(len(order))]
            for new_offset, old_offset in enumerate(order):
                mo_coeffs[mo_idx][start + new_offset] = block[old_offset]

    @staticmethod
    def _permute_square_block(matrix, start, order):
        size = len(order)
        for row in range(len(matrix)):
            block = [matrix[row][start + offset] for offset in range(size)]
            for new_offset, old_offset in enumerate(order):
                matrix[row][start + new_offset] = block[old_offset]

        rows = [matrix[start + offset][:] for offset in range(size)]
        for new_offset, old_offset in enumerate(order):
            matrix[start + new_offset] = rows[old_offset]

    @staticmethod
    def _extend_gaussian_spin_blocks(output, coefficients, beta_sign):
        scale = 1.0 / math.sqrt(2.0)
        output.extend(coeff * scale for coeff in coefficients)
        output.extend(beta_sign * coeff * scale for coeff in coefficients)

    def _orca_json_file_needs_update(self, json_file):
        if not json_file.exists():
            return True

        source_files = [
            self.out_file_name.with_suffix(suffix)
            for suffix in ('.gbw', '.cis', '.property.txt', '.json.conf')
        ]
        source_files = [path for path in source_files if path.exists()]
        if not source_files:
            return False

        json_mtime = json_file.stat().st_mtime
        return any(path.stat().st_mtime > json_mtime for path in source_files)

    def _write_orca_json_config(self):
        config_file = self.out_file_name.with_suffix('.json.conf')
        config = {}

        if config_file.exists():
            try:
                with config_file.open('r') as handle:
                    config = json.load(handle)
            except json.JSONDecodeError:
                warnings.warn(
                    f"Could not parse {config_file.name}; replacing it with "
                    "a minimal orca_2json config for TD vectors.",
                    RuntimeWarning,
                )

        required_config = {
            "CIS": True,
            "CISNRoots": True,
            "JSONFormats": ["json"],
            "MOCoefficients": False,
            "Basisset": False,
        }
        updated_config = dict(config)
        updated_config.update(required_config)
        if config_file.exists() and updated_config == config:
            return

        config = updated_config

        with config_file.open('w') as handle:
            json.dump(config, handle, indent=2)
            handle.write("\n")

    def _ensure_orca_json(self):
        json_file = self.out_file_name.with_suffix('.json')
        self._write_orca_json_config()
        if not self._orca_json_file_needs_update(json_file):
            return json_file

        gbw_file = self.out_file_name.with_suffix('.gbw')
        if not gbw_file.exists():
            warnings.warn(
                f"ORCA GBW file {gbw_file.name} is missing; falling back to "
                "sparse CI amplitudes printed in the .out file.",
                RuntimeWarning,
            )
            return None

        try:
            subprocess.run(
                ['orca_2json', gbw_file.name, '-json'],
                cwd=self.out_file_name.parent,
                check=True,
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
        except FileNotFoundError:
            warnings.warn(
                "orca_2json was not found in PATH; falling back to sparse CI "
                "amplitudes printed in the .out file.",
                RuntimeWarning,
            )
            return None
        except subprocess.CalledProcessError as error:
            warnings.warn(
                "orca_2json failed while exporting TD vectors; falling back to "
                f"sparse CI amplitudes printed in the .out file.\n{error.stdout}",
                RuntimeWarning,
            )
            return None

        if not json_file.exists():
            warnings.warn(
                f"orca_2json completed but did not create {json_file.name}; "
                "falling back to sparse CI amplitudes printed in the .out file.",
                RuntimeWarning,
            )
            return None

        return json_file

    def _load_orca_json(self):
        json_file = self._ensure_orca_json()
        if json_file is None:
            return None

        try:
            with json_file.open('r') as handle:
                return json.load(handle)
        except (OSError, json.JSONDecodeError) as error:
            warnings.warn(
                f"Could not read ORCA JSON TD-vector data from {json_file.name}: "
                f"{error}. Falling back to sparse CI amplitudes printed in the .out file.",
                RuntimeWarning,
            )
            return None

    def _json_ci_vector(self, root, key):
        if key not in root:
            return None

        orbwin = root.get("OrbWin")
        if not isinstance(orbwin, list) or len(orbwin) != 4:
            raise Exception(f"ORCA JSON TD root {root.get('IRoot')} is missing a valid OrbWin field")

        occ_start, occ_end, virt_start, virt_end = [int(value) for value in orbwin]
        expected_occ = occ_end - occ_start + 1
        expected_virt = virt_end - virt_start + 1

        blocks = [block for block in root[key] if isinstance(block, list)]
        if len(blocks) != expected_occ:
            raise Exception(
                f"ORCA JSON TD root {root.get('IRoot')} has {len(blocks)} {key} "
                f"occupied blocks; expected {expected_occ}"
            )

        coefficients = [0.0] * self.ndim
        for occ_offset, block in enumerate(blocks):
            if len(block) != expected_virt:
                raise Exception(
                    f"ORCA JSON TD root {root.get('IRoot')} has a {key} virtual "
                    f"block of length {len(block)}; expected {expected_virt}"
                )

            occ_idx = self._active_occ_position.get(occ_start + occ_offset)
            if occ_idx is None:
                continue

            for virt_offset, amplitude in enumerate(block):
                virt_idx = self._active_virt_position.get(virt_start + virt_offset)
                if virt_idx is None:
                    continue
                transition_idx = occ_idx * self.num_virtual_orbitals + virt_idx
                coefficients[transition_idx] = float(amplitude)

        return coefficients

    def _json_state_number(self, root, states):
        state_num = int(root["IRoot"])
        known_states = {state[0] for state in states}
        if state_num in known_states:
            return state_num

        # Some ORCA versions/tools report roots local to the multiplicity block.
        if 1 <= state_num <= len(states):
            return states[state_num - 1][0]

        return state_num

    def _parse_json_ci_coefficients(self):
        data = self._load_orca_json()
        if data is None:
            return None, None

        roots = data.get("Molecule", {}).get("TD-DFT", [])
        if not roots:
            warnings.warn(
                "ORCA JSON did not contain Molecule/TD-DFT vector data; "
                "falling back to sparse CI amplitudes printed in the .out file.",
                RuntimeWarning,
            )
            return None, None

        singlets = {}
        triplets = {}
        for root in roots:
            if "IRoot" not in root or "Multiplicity" not in root:
                continue

            x_coefficients = self._json_ci_vector(root, "X")
            if x_coefficients is None:
                continue

            y_coefficients = self._json_ci_vector(root, "Y")
            if y_coefficients is None:
                xpy = x_coefficients
                xmy = x_coefficients
                if root.get("TDA") != "ON":
                    warnings.warn(
                        f"ORCA JSON root {root['IRoot']} has no Y vector even "
                        "though it is not marked as TDA; using X for both X+Y and X-Y.",
                        RuntimeWarning,
                    )
            else:
                xpy = [x + y for x, y in zip(x_coefficients, y_coefficients)]
                xmy = [x - y for x, y in zip(x_coefficients, y_coefficients)]

            multiplicity = int(root["Multiplicity"])
            if multiplicity == 1:
                singlets[self._json_state_number(root, self.singlet_states)] = (xpy, xmy)
            elif multiplicity == 3:
                triplets[self._json_state_number(root, self.triplet_states)] = (xpy, xmy)

        return singlets, triplets

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
        self._parse_mkl_coord(lines)         
        self._parse_mkl_occ(lines)

        if self._parse_output_basis():
            return

        warnings.warn(
            "Could not find ORCA's raw 'BASIS SET IN INPUT FORMAT' block; "
            "falling back to the .mkl basis coefficients. These coefficients "
            "may not match molsoc's expected GTO contraction normalization.",
            RuntimeWarning,
        )

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

    def _parse_output_basis(self):
        try:
            lines = self.out_file_name.read_text().splitlines()
        except UnicodeDecodeError:
            return False

        try:
            start_idx = next(
                idx for idx, line in enumerate(lines)
                if line.strip() == "BASIS SET IN INPUT FORMAT"
            )
        except StopIteration:
            return False

        basis_by_element = {}
        current_element = None
        idx = start_idx + 1
        while idx < len(lines):
            line = lines[idx]
            stripped = line.strip()

            if basis_by_element and current_element is None and stripped.startswith("---"):
                break

            if stripped.startswith("NewGTO"):
                parts = stripped.split()
                if len(parts) >= 2:
                    current_element = parts[1]
                    basis_by_element[current_element] = []
                idx += 1
                continue

            if current_element is not None and stripped.lower().startswith("end;"):
                current_element = None
                idx += 1
                continue

            parts = stripped.split()
            if (
                current_element is not None
                and len(parts) >= 2
                and parts[0] in self.SHELLS
                and parts[1].isdigit()
            ):
                shell_type = parts[0]
                num_primitives = int(parts[1])
                primitives = []
                for primitive_line in lines[idx + 1: idx + 1 + num_primitives]:
                    primitive_parts = primitive_line.split()
                    if len(primitive_parts) < 3:
                        return False
                    exponent, coefficient = primitive_parts[-2:]
                    primitives.append(f"   {exponent:>18}  {coefficient:>16}\n")
                basis_by_element[current_element].append((shell_type, primitives))
                idx += num_primitives + 1
                continue

            idx += 1

        if not basis_by_element:
            return False

        basis_set = []
        ao_basis = []
        for atom in self.geometry:
            element = atom[0]
            if element not in basis_by_element:
                return False

            basis_set.append(f"{element}  0\n")
            for shell_type, primitives in basis_by_element[element]:
                basis_set.append(f"{shell_type}   {len(primitives)} 1.0\n")
                basis_set.extend(primitives)
                ao_basis.append(self.SHELLS[shell_type])
            basis_set.append("****\n")

        if basis_set and basis_set[-1] == "****\n":
            basis_set.pop()

        self.basis_set = basis_set
        self.ao_basis = ao_basis
        return True

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
        
        #ci_dict_singlets = {i: [0.0]*self.ndim for i in range(1, len(self.singlet_states)+1)}
        #ci_dict_triplets = {i: [0.0]*self.ndim for i in range(1, len(self.triplet_states)+1)}
        ci_dict_singlets = {}
        ci_dict_triplets = {}
        current_state = None
        
        i = 0
        while i < len(lines):
            line = lines[i]
            i += 1

            self._parse_orbital_range(line)
            self._parse_root_count(line)
            
            # 2. OVERLAP MATRIX
            if "OVERLAP MATRIX" in line:
                in_overlap = True
                overlap_matrix = [[0.0]*self.num_orbitals for _ in range(self.num_orbitals)]
                continue
            if in_overlap:
                if "Time for" in line or "ORBITAL ENERGIES" in line or "MOLECULAR ORBITALS" in line or "DIPOLE" in line or "D-I-I-S" in line:
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
            if "EXCITED STATES (SINGLETS)" in line:
                in_singlets = True
                in_triplets = False
            elif "EXCITED STATES (TRIPLETS)" in line:
                in_singlets = False
                in_triplets = True
            elif "Entering triplet calculation" in line or "RPA-DIAGONALIZATION" in line:
                in_singlets = False
            
            if in_singlets or in_triplets:
                if "the weight of the individual excitations are printed if larger than" in line:
                    match = re.search(r"larger than\s+([0-9.+\-Ee]+)", line)
                    if match is not None:
                        try:
                            self._sparse_ci_threshold = float(match.group(1))
                        except ValueError:
                            pass

                if "Storing amplitudes" in line:
                    in_singlets = False
                    in_triplets = False
                    continue
                if line.startswith("STATE "):
                    current_state = int(line.split()[1][:-1])
                        # example: "STATE  1:  E=   0.149148 au      4.059 eV    32734.2 cm**-1 ..."
                    parts = line.split()
                    # find the eV value: after 'au' there is something and then 'eV'
                    try:
                        # locate the index of 'eV' and take the previous element
                        idx_ev = parts.index('eV')
                        if idx_ev > 0:
                            energy_ev = float(parts[idx_ev-1])
                        else:
                            continue
                    except (ValueError, IndexError):
                        continue

                    if in_singlets:
                        self.singlet_states.append([current_state, energy_ev])
                        if current_state not in ci_dict_singlets:
                            ci_dict_singlets[current_state] = [0.0] * self.ndim
                    elif in_triplets:
                        self.triplet_states.append([current_state, energy_ev])
                        if current_state not in ci_dict_triplets:
                            ci_dict_triplets[current_state] = [0.0] * self.ndim

                elif "->" in line and ":" in line:
                    parts = line.split()
                    occ = int(parts[0][:-1]) 
                    virt = int(parts[2][:-1])
                    
                    # Extract the CI amplitude (c=...) – the number is the token
                    # immediately before the closing parenthesis.
                    amplitude = None
                    for idx, token in enumerate(parts):
                        if token.startswith('(c='):
                            # The coefficient may be in the same token or the next one
                            if token.endswith(')'):
                                amp_str = token[3:-1]   # remove 'c=' and ')'
                            else:
                                # the number is the next token, which ends with ')'
                                amp_str = parts[idx+1].rstrip(')')
                            try:
                                amplitude = float(amp_str)
                            except ValueError:
                                pass
                            break
                    if amplitude is None:
                        try:
                            weight = float(parts[4])
                            amplitude = math.sqrt(weight) if weight >= 0.0 else 0.0
                        except (IndexError, ValueError):
                            continue
                    
                    occ_idx = self._active_occ_position.get(occ)
                    virt_idx = self._active_virt_position.get(virt)
                    
                    if occ_idx is not None and virt_idx is not None:
                        transition_idx = occ_idx * self.num_virtual_orbitals + virt_idx
                        if in_singlets and current_state in ci_dict_singlets:
                            ci_dict_singlets[current_state][transition_idx] = amplitude
                        elif in_triplets and current_state in ci_dict_triplets:
                            ci_dict_triplets[current_state][transition_idx] = amplitude

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
                
        # Permute ORCA spherical harmonic shells to match the order produced
        # by the Cartesian -> spherical transform in soc_td/basis_match.f90.
        p_reorder = [2, 0, 1]  # ORCA: pz, px, py -> PySOC: py, pz, px

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
            if 'pz' in name:
                self._permute_mo_ao_block(mo_coeffs, ao_idx, p_reorder)
                ao_idx += 3
            elif 'dz2' in name:
                self._permute_mo_ao_block(mo_coeffs, ao_idx, d_reorder)
                ao_idx += 5
            elif 'fz3' in name:
                self._permute_mo_ao_block(mo_coeffs, ao_idx, f_reorder)
                ao_idx += 7
            else:
                ao_idx += 1
          
        # Also permute the spherical AO overlaps to match MO coefficient order
        n = self.num_orbitals
        overlap2d = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                overlap2d[i][j] = self.AO_overlaps[i*n + j]

        ao_idx = 0
        while ao_idx < len(aonames):
            name = aonames[ao_idx].lower()
            if 'pz' in name:
                self._permute_square_block(overlap2d, ao_idx, p_reorder)
                ao_idx += 3
            elif 'dz2' in name:
                self._permute_square_block(overlap2d, ao_idx, d_reorder)
                ao_idx += 5
            elif 'fz3' in name:
                self._permute_square_block(overlap2d, ao_idx, f_reorder)
                ao_idx += 7
            else:
                ao_idx += 1

        self.AO_overlaps = [overlap2d[i][j] for i in range(n) for j in range(n)]

        active_mo_indices = self.active_occupied_orbitals + self.active_virtual_orbitals
        self.MO_energies = [mo_energies[mo_idx] for mo_idx in active_mo_indices]

        self.MOA_coefficients = []
        for mo_idx in active_mo_indices:
            self.MOA_coefficients.extend(mo_coeffs[mo_idx])
            
        # CI Coeffs
        # soc_td with 'gauss_tddft' expects 2 * ndim coefficients per state for X+Y, 
        # followed by 2 * ndim coefficients per state for X-Y.
        self._warn_if_selected_edge_roots()

        json_ci_singlets, json_ci_triplets = self._parse_json_ci_coefficients()
        using_json_ci = json_ci_singlets is not None and json_ci_triplets is not None

        ci_xpy = []
        ci_xmy = []
        
        # Singlets (state numbers usually 1‑n_singl, but map generically)
        for i in self.requested_singlets:
            state_num = self.singlet_states[i-1][0]   # actual state number
            if using_json_ci:
                if state_num not in json_ci_singlets:
                    raise Exception(f"Missing ORCA JSON CI coefficients for singlet state {state_num}")
                xpy_coefficients, xmy_coefficients = json_ci_singlets[state_num]
            else:
                if state_num not in ci_dict_singlets:
                    raise Exception(f"Missing ORCA CI coefficients for singlet state {state_num}")
                coefficients = ci_dict_singlets[state_num]
                if not any(abs(coeff) > 0.0 for coeff in coefficients):
                    raise Exception(f"ORCA printed no CI coefficients for singlet state {state_num}")
                xpy_coefficients = coefficients
                xmy_coefficients = coefficients
            self._extend_gaussian_spin_blocks(ci_xpy, xpy_coefficients, beta_sign=1.0)
            self._extend_gaussian_spin_blocks(ci_xmy, xmy_coefficients, beta_sign=1.0)

        # Triplets (state numbers can be > n_singl)
        for i in self.requested_triplets:
            state_num = self.triplet_states[i-1][0]
            if using_json_ci:
                if state_num not in json_ci_triplets:
                    raise Exception(f"Missing ORCA JSON CI coefficients for triplet state {state_num}")
                xpy_coefficients, xmy_coefficients = json_ci_triplets[state_num]
            else:
                if state_num not in ci_dict_triplets:
                    raise Exception(f"Missing ORCA CI coefficients for triplet state {state_num}")
                coefficients = ci_dict_triplets[state_num]
                if not any(abs(coeff) > 0.0 for coeff in coefficients):
                    raise Exception(f"ORCA printed no CI coefficients for triplet state {state_num}")
                xpy_coefficients = coefficients
                xmy_coefficients = coefficients
            self._extend_gaussian_spin_blocks(ci_xpy, xpy_coefficients, beta_sign=-1.0)
            self._extend_gaussian_spin_blocks(ci_xmy, xmy_coefficients, beta_sign=-1.0)

        if not using_json_ci and self._sparse_ci_threshold is not None and self._sparse_ci_threshold > 0.0:
            warnings.warn(
                "ORCA printed only CI amplitudes above "
                f"{self._sparse_ci_threshold:g}; the parsed TD vectors are incomplete. "
                "Use a full-vector ORCA source for quantitative SOC values.",
                RuntimeWarning,
            )

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
