#!/usr/bin/env python3
"""
Boresch-style restraint helper
======================================================

It will:
  • Parse your PDB or GRO file to locate atom indices and coordinates.
  • Handle proper chain IDs for PDB files.
  • Compute distance (L1–P1), angles (L2–L1–P1, L1–P1–P2), and dihedrals
    (L3–L2–L1–P1, L2–L1–P1–P2, L1–P1–P2–P3).
  • Emit a GROMACS-format snippet under [ intermolecular_interactions ].

"""

from __future__ import annotations
import math
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

# ================================ CONFIG =================================== #
# Ligand name as it appears in the resname column (e.g., LDP, SKF, LY3, BMS)
LIGAND_NAME = "SKF"
# Optional: if multiple ligand copies exist, fix which one to use (integer resid or None)
LIGAND_RESID: Optional[int] = None

# Ligand atoms (names must match the atomname for the chosen ligand residue)
L1 = "N1"   
L2 = "O1"   
L3 = "O2"

# Protein atoms with (atom_name, residue_number, chain)
# For PDB: chain is the actual chain ID (A, B, R, etc.)
P1 = ("CA", "103", "R")  # e.g., CA of residue 103, chain R
P2 = ("CA", "198", "R")  # e.g., CB of residue 198, chain R
P3 = ("CA", "199", "R")  # e.g., CB of residue 199, chain R

# Input file
INPUT_FILE = "starting_structure.pdb"
# Optional output file. If None, snippet is printed to stdout.
OUTPUT_FILE: Optional[str] = None  # e.g., "boresch_snippet.itp"

# Force constants and types (tweak as desired)
FORCE_CONST = {
    "bond_kB": 418.4,      # bonds fcB
    "angle_kB": 209.2,     # angles fcB
    "dihedral_kB": 209.2,  # dihedrals fcB
}
TYPES = {
    "bond_type": 6,
    "angle_type": 1,
    "dihedral_type": 2,
}
# Formatting for printed values
DIST_FMT = "{:.3f}"  # nm
ANGL_FMT = "{:.1f}"  # degrees
# ============================== END CONFIG ================================= #

@dataclass
class Atom:
    index: int           # 1-based index
    resid: int
    resname: str
    atomname: str
    chain: str           # actual chain ID for PDB, heuristic for GRO
    x: float
    y: float
    z: float
    @property
    def pos(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

# ------------------------ file parsing ------------------------ #

def parse_pdb(path: str) -> List[Atom]:
    """Parse PDB file with proper chain handling and robust coordinate parsing."""
    atoms: List[Atom] = []
    index = 0
    
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    index += 1  # Sequential numbering since PDB atom numbers can be non-sequential
                    
                    atomname = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() if len(line) > 21 else ""
                    resid = int(line[22:26].strip())
                    
                    # Robust coordinate parsing to handle malformed PDB files
                    x_str = line[30:38].strip()
                    y_str = line[38:46].strip()
                    z_str = line[46:54].strip()
                    
                    # Fix common formatting issues like "126.540-" -> "126.540 -"
                    x_str = fix_coordinate_format(x_str)
                    y_str = fix_coordinate_format(y_str)
                    z_str = fix_coordinate_format(z_str)
                    
                    x = float(x_str) / 10.0  # Convert Å to nm
                    y = float(y_str) / 10.0
                    z = float(z_str) / 10.0
                    
                    atoms.append(Atom(index, resid, resname, atomname, chain, x, y, z))
                    
                except (ValueError, IndexError) as e:
                    print(f"[warn] Skipping malformed PDB line {line_num}: {line.strip()}", file=sys.stderr)
                    print(f"       Error: {e}", file=sys.stderr)
                    index -= 1  # Don't increment index for failed parsing
                    continue
    
    return atoms

def fix_coordinate_format(coord_str: str) -> str:
    """Fix common PDB coordinate formatting issues."""
    if not coord_str:
        return "0.0"
    
    # Handle cases like "126.540-" -> "126.540"
    if coord_str.endswith('-') and coord_str[:-1].replace('.', '').isdigit():
        return '-' + coord_str[:-1]
    
    # Handle cases like "126.540+" -> "126.540"
    if coord_str.endswith('+') and coord_str[:-1].replace('.', '').isdigit():
        return coord_str[:-1]
    
    # Handle missing decimal points or other issues
    if coord_str.replace('-', '').replace('.', '').isdigit():
        return coord_str
    
    # Try to extract valid number from string
    import re
    match = re.search(r'-?\d+\.?\d*', coord_str)
    if match:
        return match.group()
    
    # If all else fails, return 0.0 and let the caller handle it
    print(f"[warn] Could not parse coordinate '{coord_str}', using 0.0", file=sys.stderr)
    return "0.0"

def parse_gro(path: str) -> List[Atom]:
    """Parse GRO file with heuristic chain handling (legacy)."""
    atoms: List[Atom] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        _title = f.readline()
        n_line = f.readline()
        if not n_line:
            raise ValueError("Invalid .gro file: missing atom count line.")
        try:
            n_atoms = int(n_line.strip())
        except ValueError as e:
            raise ValueError("Invalid .gro file: atom count is not an integer.") from e
        
        for _ in range(n_atoms):
            line = f.readline()
            if not line:
                raise ValueError("Unexpected end of .gro while reading atoms.")
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            index = int(line[15:20])
            x = float(line[20:28])  # Already in nm
            y = float(line[28:36])
            z = float(line[36:44])
            
            # Heuristic chain extraction for GRO
            chain = extract_chain_from_resname(resname)
            
            atoms.append(Atom(index, resid, resname, atomname, chain, x, y, z))
    return atoms

def extract_chain_from_resname(resname: str) -> str:
    """Extract chain hint from GRO resname (heuristic)."""
    rn = resname.upper()
    # Look for common patterns: ALAr -> r, ALA_R -> R, ALA-R -> R
    if len(rn) > 3 and rn[-1].isalpha() and rn[-1] not in ['A', 'T', 'G', 'C']:
        return rn[-1]
    if "_" in rn:
        parts = rn.split("_")
        if len(parts) > 1 and len(parts[-1]) == 1:
            return parts[-1]
    if "-" in rn:
        parts = rn.split("-")
        if len(parts) > 1 and len(parts[-1]) == 1:
            return parts[-1]
    return ""  # No chain hint found

def parse_file(path: str) -> List[Atom]:
    """Auto-detect file type and parse accordingly."""
    if path.lower().endswith('.pdb'):
        return parse_pdb(path)
    elif path.lower().endswith('.gro'):
        return parse_gro(path)
    else:
        # Try to detect by content
        with open(path, "r") as f:
            first_line = f.readline()
            if first_line.startswith(("ATOM", "HETATM")):
                return parse_pdb(path)
            else:
                return parse_gro(path)

# ------------------------ selection helpers ------------------------ #

def match_chain_pdb(atom_chain: str, target_chain: Optional[str]) -> bool:
    """Exact chain matching for PDB files."""
    if not target_chain:
        return True
    return atom_chain.strip().upper() == target_chain.strip().upper()

def match_chain_gro(resname: str, target_chain: Optional[str]) -> bool:
    """Heuristic chain matching for GRO files (legacy)."""
    if not target_chain:
        return True
    c = target_chain.strip().upper()
    rn = resname.strip().upper()
    return rn.endswith(c) or ("_" + c) in rn or ("-" + c) in rn

def find_ligand_atom(atoms: List[Atom], ligand_name: str, atom_name: str, ligand_resid: Optional[int] = None) -> Atom:
    cands = [a for a in atoms if a.resname == ligand_name and a.atomname == atom_name]
    if ligand_resid is not None:
        cands = [a for a in cands if a.resid == ligand_resid]
    if not cands:
        hint = f" (resid {ligand_resid})" if ligand_resid is not None else ""
        raise ValueError(f"Could not find ligand atom {ligand_name}:{atom_name}{hint}")
    if len(cands) > 1 and ligand_resid is None:
        print(f"[warn] Multiple {ligand_name}:{atom_name} found; using the first occurrence. Consider setting LIGAND_RESID.", file=sys.stderr)
    return cands[0]

def find_protein_atom(atoms: List[Atom], resid: int, atom_name: str, chain: Optional[str], is_pdb: bool = True) -> Atom:
    """Find protein atom with proper chain handling."""
    if is_pdb:
        # PDB: use exact chain matching
        cands = [a for a in atoms if a.resid == resid and a.atomname == atom_name and match_chain_pdb(a.chain, chain)]
    else:
        # GRO: use heuristic chain matching, with fallback
        cands = [a for a in atoms if a.resid == resid and a.atomname == atom_name and match_chain_gro(a.resname, chain)]
        if not cands:
            # Fallback: ignore chain altogether (legacy GRO behavior)
            cands = [a for a in atoms if a.resid == resid and a.atomname == atom_name]
    
    if not cands:
        chain_txt = f" chain {chain}" if chain else ""
        raise ValueError(f"Could not find protein atom {atom_name} in residue {resid}{chain_txt}.")
    
    if len(cands) > 1:
        if is_pdb:
            # For PDB, this shouldn't happen with proper chain IDs
            print(f"[warn] Multiple matches for residue {resid} atom {atom_name} chain {chain}; using first.", file=sys.stderr)
        else:
            print(f"[warn] Multiple matches for residue {resid} atom {atom_name}; using first.", file=sys.stderr)
    
    return cands[0]

def residue_name(atoms: List[Atom], resid: int, chain: Optional[str], is_pdb: bool = True) -> Optional[str]:
    """Get residue name with proper chain handling."""
    if is_pdb:
        cands = [a for a in atoms if a.resid == resid and match_chain_pdb(a.chain, chain)]
    else:
        cands = [a for a in atoms if a.resid == resid and match_chain_gro(a.resname, chain)]
        if not cands:
            cands = [a for a in atoms if a.resid == resid]
    return cands[0].resname if cands else None

# ------------------------ geometry ------------------------ #

def vec(a: Atom, b: Atom) -> Tuple[float, float, float]:
    return (b.x - a.x, b.y - a.y, b.z - a.z)

def vdot(u: Tuple[float, float, float], v: Tuple[float, float, float]) -> float:
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def vnorm(u: Tuple[float, float, float]) -> float:
    return math.sqrt(vdot(u, u))

def vcross(u: Tuple[float, float, float], v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])

def distance(a: Atom, b: Atom) -> float:
    return vnorm(vec(a, b))

def angle(a: Atom, b: Atom, c: Atom) -> float:
    # angle at b formed by a-b-c
    ba = vec(b, a)
    bc = vec(b, c)
    nba = vnorm(ba)
    nbc = vnorm(bc)
    if nba == 0 or nbc == 0:
        return float('nan')
    cosang = max(-1.0, min(1.0, vdot(ba, bc) / (nba * nbc)))
    return math.degrees(math.acos(cosang))

def dihedral(a: Atom, b: Atom, c: Atom, d: Atom, debug: bool = False) -> float:
    """
    Calculate dihedral angle for atoms a-b-c-d.
    Uses the standard definition: angle between planes (a,b,c) and (b,c,d).
    """
    # Vectors along the bonds
    b1 = vec(a, b)  # a -> b
    b2 = vec(b, c)  # b -> c  
    b3 = vec(c, d)  # c -> d
    
    # Normal vectors to the two planes
    n1 = vcross(b1, b2)  # Normal to plane (a,b,c)
    n2 = vcross(b2, b3)  # Normal to plane (b,c,d)
    
    # Normalize the normal vectors
    n1_norm = vnorm(n1)
    n2_norm = vnorm(n2)
    
    if n1_norm == 0 or n2_norm == 0:
        return float('nan')
    
    n1_unit = (n1[0]/n1_norm, n1[1]/n1_norm, n1[2]/n1_norm)
    n2_unit = (n2[0]/n2_norm, n2[1]/n2_norm, n2[2]/n2_norm)
    
    # Calculate the angle between normal vectors
    cos_theta = max(-1.0, min(1.0, vdot(n1_unit, n2_unit)))
    
    # Determine sign using the middle bond vector b2
    b2_norm = vnorm(b2)
    if b2_norm == 0:
        return float('nan')
    
    b2_unit = (b2[0]/b2_norm, b2[1]/b2_norm, b2[2]/b2_norm)
    sin_theta = vdot(vcross(n1_unit, n2_unit), b2_unit)
    
    angle_rad = math.atan2(sin_theta, cos_theta)
    angle_deg = math.degrees(angle_rad)
    
    if debug:
        print(f"Debug dihedral {a.atomname}-{b.atomname}-{c.atomname}-{d.atomname}:")
        print(f"  a: ({a.x:.3f}, {a.y:.3f}, {a.z:.3f})")
        print(f"  b: ({b.x:.3f}, {b.y:.3f}, {b.z:.3f})")
        print(f"  c: ({c.x:.3f}, {c.y:.3f}, {c.z:.3f})")
        print(f"  d: ({d.x:.3f}, {d.y:.3f}, {d.z:.3f})")
        print(f"  b1 (a->b): ({b1[0]:.3f}, {b1[1]:.3f}, {b1[2]:.3f})")
        print(f"  b2 (b->c): ({b2[0]:.3f}, {b2[1]:.3f}, {b2[2]:.3f})")
        print(f"  b3 (c->d): ({b3[0]:.3f}, {b3[1]:.3f}, {b3[2]:.3f})")
        print(f"  n1: ({n1[0]:.3f}, {n1[1]:.3f}, {n1[2]:.3f})")
        print(f"  n2: ({n2[0]:.3f}, {n2[1]:.3f}, {n2[2]:.3f})")
        print(f"  cos_theta: {cos_theta:.3f}")
        print(f"  sin_theta: {sin_theta:.3f}")
        print(f"  angle: {angle_deg:.1f}°")
    
    return angle_deg

# ------------------------ formatting ------------------------ #

def format_boresch_snippet(idx: Dict[str, int], geom: Dict[str, float]) -> str:
    lines: List[str] = []
    lines.append("[ intermolecular_interactions]")

    # Bonds (distance): L1 - P1
    lines.append("[ bonds ]")
    lines.append("; ai     aj    type   bA      kA     bB      kB")
    b = geom["d_L1_P1"]
    lines.append(
        f"  {idx['L1']:>6}  {idx['P1']:>6}  {TYPES['bond_type']:<2}   {DIST_FMT.format(b):>6}   0.0    {DIST_FMT.format(b):>6}   {FORCE_CONST['bond_kB']}"
    )
    lines.append("")

    # Angles: L2-L1-P1 and L1-P1-P2
    lines.append("[ angles ]")
    lines.append("; ai     aj    ak     type    thA      fcA        thB      fcB")
    a1 = geom["ang_L2_L1_P1"]
    a2 = geom["ang_L1_P1_P2"]
    lines.append(
        f"  {idx['L2']:>6}  {idx['L1']:>6}  {idx['P1']:>6}   {TYPES['angle_type']:<1}     {ANGL_FMT.format(a1):>5}      0.0       {ANGL_FMT.format(a1):>5}      {FORCE_CONST['angle_kB']}"
    )
    lines.append(
        f"  {idx['L1']:>6}  {idx['P1']:>6}  {idx['P2']:>6}   {TYPES['angle_type']:<1}     {ANGL_FMT.format(a2):>5}      0.0       {ANGL_FMT.format(a2):>5}      {FORCE_CONST['angle_kB']}"
    )
    lines.append("")

    # Dihedrals: L3-L2-L1-P1, L2-L1-P1-P2, L1-P1-P2-P3
    lines.append("[ dihedrals ]")
    lines.append("; ai     aj    ak    al    type     thA      fcA       thB      fcB")
    d1 = geom["dih_L3_L2_L1_P1"]
    d2 = geom["dih_L2_L1_P1_P2"]
    d3 = geom["dih_L1_P1_P2_P3"]
    lines.append(
        f"  {idx['L3']:>6}  {idx['L2']:>6}  {idx['L1']:>6}  {idx['P1']:>6}  {TYPES['dihedral_type']:<1}      {ANGL_FMT.format(d1):>5}     0.0    {ANGL_FMT.format(d1):>7}      {FORCE_CONST['dihedral_kB']}"
    )
    lines.append(
        f"  {idx['L2']:>6}  {idx['L1']:>6}  {idx['P1']:>6}  {idx['P2']:>6}  {TYPES['dihedral_type']:<1}      {ANGL_FMT.format(d2):>5}     0.0    {ANGL_FMT.format(d2):>7}      {FORCE_CONST['dihedral_kB']}"
    )
    lines.append(
        f"  {idx['L1']:>6}  {idx['P1']:>6}  {idx['P2']:>6}  {idx['P3']:>6}  {TYPES['dihedral_type']:<1}      {ANGL_FMT.format(d3):>5}     0.0    {ANGL_FMT.format(d3):>7}      {FORCE_CONST['dihedral_kB']}"
    )
    return "\n".join(lines) + "\n"

# ------------------------ main flow ------------------------ #

def main() -> None:
    atoms = parse_file(INPUT_FILE)
    is_pdb = INPUT_FILE.lower().endswith('.pdb')
    
    print(f"Parsed {len(atoms)} atoms from {INPUT_FILE} ({'PDB' if is_pdb else 'GRO'} format)")

    # Resolve ligand atoms
    try:
        L1_atom = find_ligand_atom(atoms, LIGAND_NAME, L1, LIGAND_RESID)
        L2_atom = find_ligand_atom(atoms, LIGAND_NAME, L2, LIGAND_RESID)
        L3_atom = find_ligand_atom(atoms, LIGAND_NAME, L3, LIGAND_RESID)
    except Exception as e:
        print(f"[error] Ligand selection failed: {e}", file=sys.stderr)
        raise

    # Resolve protein atoms
    def parse_P(t: Tuple[str, str, Optional[str]]) -> Tuple[str, int, Optional[str]]:
        atom, resid_str, chain = t
        return atom, int(resid_str), (chain if chain else None)

    P1_atomname, P1_resid, P1_chain = parse_P(P1)
    P2_atomname, P2_resid, P2_chain = parse_P(P2)
    P3_atomname, P3_resid, P3_chain = parse_P(P3)

    try:
        P1_atom = find_protein_atom(atoms, P1_resid, P1_atomname, P1_chain, is_pdb)
        P2_atom = find_protein_atom(atoms, P2_resid, P2_atomname, P2_chain, is_pdb)
        P3_atom = find_protein_atom(atoms, P3_resid, P3_atomname, P3_chain, is_pdb)
    except Exception as e:
        print(f"[error] Protein selection failed: {e}", file=sys.stderr)
        raise

    # Sanity: print residue names
    for resid, chain in [(P1_resid, P1_chain), (P2_resid, P2_chain), (P3_resid, P3_chain)]:
        rn = residue_name(atoms, resid, chain, is_pdb)
        rn_txt = rn if rn is not None else "<not found>"
        ch_txt = f", chain={chain}" if chain else ""
        print(f"Residue {resid}{ch_txt}: {rn_txt}")

    # Indices
    idx = {
        "L1": L1_atom.index, "L2": L2_atom.index, "L3": L3_atom.index,
        "P1": P1_atom.index, "P2": P2_atom.index, "P3": P3_atom.index,
    }

    # Geometry
    geom = {
        "d_L1_P1": distance(L1_atom, P1_atom),
        "ang_L2_L1_P1": angle(L2_atom, L1_atom, P1_atom),
        "ang_L1_P1_P2": angle(L1_atom, P1_atom, P2_atom),
        "dih_L3_L2_L1_P1": dihedral(L3_atom, L2_atom, L1_atom, P1_atom, debug=True),
        "dih_L2_L1_P1_P2": dihedral(L2_atom, L1_atom, P1_atom, P2_atom, debug=True),
        "dih_L1_P1_P2_P3": dihedral(L1_atom, P1_atom, P2_atom, P3_atom, debug=True),
    }

    snippet = format_boresch_snippet(idx, geom)

    # Output
    if OUTPUT_FILE:
        with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
            f.write(snippet)
        print(f"\nWrote GROMACS snippet to: {OUTPUT_FILE}")
    else:
        print("\n" + snippet)

    # Summary
    print("Selections used:")
    print(f"  Ligand {LIGAND_NAME} (resid={LIGAND_RESID if LIGAND_RESID is not None else 'first match'})")
    for tag, a in (('L1', L1_atom), ('L2', L2_atom), ('L3', L3_atom), ('P1', P1_atom), ('P2', P2_atom), ('P3', P3_atom)):
        chain_info = f" chain={a.chain}" if a.chain else ""
        print(f"   {tag}: index={a.index:>6}  resid={a.resid:<5} resname={a.resname:<6} atom={a.atomname}{chain_info}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)
