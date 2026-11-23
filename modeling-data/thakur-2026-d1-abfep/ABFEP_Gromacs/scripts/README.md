# ABFEP setup and restraint scripts

These Python scripts document the preparation and restraint calculations used
for the released GROMACS ABFEP systems. They use only the Python standard
library.

The deposited systems under `../systems/` are already prepared and ready for
GROMACS preprocessing. Do **not** run the charge-modification script on those
released coordinates again.

## Scripts

### [`modify_GRO_file_ABFEP_charge_setup.py`](modify_GRO_file_ABFEP_charge_setup.py)

Prepares the coordinate component of the charge-compensation scheme used for a
positively charged ligand. It edits `solvated.gro` **in place** by converting
the last two complete `TP3` waters:

1. the second-last water becomes an additional `Na+` ion; and
2. the last water becomes the alchemical molecule `ION` (`cl`, `du1`, and
   `du2`).

The two hydrogens removed when the second-last water becomes `Na+` reduce the
GRO atom count by two. The script assumes that the final two `TP3` residues are
complete three-atom waters.

For dopamine, the endpoint charge bookkeeping is:

| Component | State A | State B |
|---|---:|---:|
| Dopamine | +1 | 0 |
| Added `Na+` | +1 | +1 |
| Alchemical `ION` | -1 | 0 |
| These three components | +1 | +1 |

The rest of the initially neutralized system has charge -1, so the complete
simulation box is neutral at both endpoints. Equivalently, dopamine changes by
-1 while `ION` changes by +1, giving a net charge change of zero.

Run it only on a copy of the unmodified coordinate file:

```bash
python3 modify_GRO_file_ABFEP_charge_setup.py
```

This script modifies coordinates only. The matching topology must separately
include the additional `Na+`, the alchemical `ION`, `ions-cl.itp`, and the
required `IM` and `du` atom types. These topology components are already
included in each applicable released condition under
`../systems/protein_bound/<condition>/toppar/`.

After preparation, verify the residue identities, atom count, topology molecule
counts, and total charge at both alchemical endpoints.

### [`boresch_restraint.py`](boresch_restraint.py)

Calculates the Boresch restraint geometry from a PDB or GRO structure and
prints a GROMACS `[ intermolecular_interactions ]` block containing one
distance, two angles, and three dihedrals.

Before running it, edit the configuration section near the top of the script:

- ligand residue and ligand anchor atoms (`L1`, `L2`, and `L3`);
- protein anchor atoms (`P1`, `P2`, and `P3`);
- input structure and optional output filename; and
- restraint force constants, if required.

```bash
python3 boresch_restraint.py
```

Check every selected atom index and the resulting geometry before inserting the
printed block into a topology.

### [`restraint_free_energy.py`](restraint_free_energy.py)

Evaluates the analytical standard-state correction for the six Boresch
restraints. Edit the temperature, equilibrium distance and angles, and six
force constants in the input section, then run:

```bash
python3 restraint_free_energy.py
```

The script reports `dG_off` and `dG_on` in kcal/mol. Use the sign required by
the thermodynamic cycle rather than adding both values.

