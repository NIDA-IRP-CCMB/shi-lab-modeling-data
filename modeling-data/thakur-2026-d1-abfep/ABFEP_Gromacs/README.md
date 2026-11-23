# D1R absolute binding free-energy data

This package contains the molecular inputs, alchemical parameters, analysis
outputs, and summary tables used for the D1R absolute binding free-energy
(ABFEP) calculations reported in the associated manuscript.

## Scope

The release contains the selected calculations used in the manuscript. Each
starting structure is identified by an anonymous name such as `system_001`.
Simulation trajectory numbers and snapshot times are intentionally not encoded
in the public names. Raw molecular-dynamics trajectories and checkpoints are
not included.

The private mapping from the anonymous names to the laboratory source folders
is maintained separately and is not part of this public package.

## System nomenclature

The D1R systems were prepared with dopamine or SKF81297 as the orthosteric
agonist and with no, one, or three PAMs. Apomorphine was included only as a
PAM-free benchmark. System identifiers use the form
`D1_<ligand/PAM composition>_<G-protein state>`. The orthosteric ligand and
PAM codes are `d` = dopamine, `a` = apomorphine, `s` = SKF81297, `l` =
LY3154207, `u9` = UNC9815, `u10` = UNC10062, and `b` = BMS-A1. The suffix
`_Gs` denotes a receptor coupled to the heterotrimeric Gs protein, whereas
`_noGs` denotes a Gs-free receptor.

The PAM-free systems are designated `D1_d_Gs` and `D1_d_noGs` for
dopamine-bound D1R, `D1_a_Gs` for apomorphine-bound D1R, and `D1_s_Gs` and
`D1_s_noGs` for SKF81297-bound D1R. The corresponding single-PAM systems
containing LY3154207 are designated `D1_dl_Gs` and `D1_dl_noGs`, and
`D1_sl_Gs` and `D1_sl_noGs`, respectively. The triple-PAM systems are
designated `D1_dlu9b_Gs` and `D1_dlu9b_noGs` (dopamine, LY3154207, UNC9815,
and BMS-A1), `D1_dlu10b_Gs` (dopamine, LY3154207, UNC10062, and BMS-A1), and
`D1_slu10b_Gs` and `D1_slu10b_noGs` (SKF81297, LY3154207, UNC10062, and
BMS-A1).

The source workbook label `D1_d_meta_up_Gs` corresponds to the public system
name `D1_d_Gs` used in the manuscript and directory structure. The additional
metadynamics-derived condition `D1_d_meta_down_Gs` remains separate.

## Directory layout

```text
ABFEP_Gromacs/
├── README.md
├── MDP/
│   ├── protein_bound/
│   └── bulk_solvent/
├── scripts/
│   ├── boresch_restraint.py
│   ├── modify_GRO_file_ABFEP_charge_setup.py
│   └── restraint_free_energy.py
├── systems/
│   ├── protein_bound/
│   │   └── <condition>/
│   │       ├── toppar/
│   │       ├── index.ndx
│   │       └── system_<NNN>/
│   └── bulk_solvent/
│       ├── apomorphine/
│       ├── dopamine/
│       └── skf81297/
└── analysis/
    ├── abfep_results.xls
    ├── protein_bound/<condition>/system_<NNN>/
    ├── bulk_solvent/
```

## Molecular inputs

Each `systems/protein_bound/<condition>/system_<NNN>/` directory contains its
authoritative compressed coordinates (`solvated.gro.gz`) and system-specific
`solvated.top`.
Regenerable PDB copies are omitted. Identical topology includes and indexes are
stored once at the condition level under `toppar/` and `index.ndx`.
The `ligand` and `system_name` columns in `analysis/abfep_results.xls`
identify the ligand and anonymous calculation. See
[`systems/protein_bound/README.md`](systems/protein_bound/README.md) for
decompression and GROMACS usage.

All protein-bound systems use the shared 64-window templates and generation
scripts stored directly under `MDP/protein_bound/`. The ligand molecule type
is standardized as `LIG` in every released topology and in the shared MDP
files.

The `systems/bulk_solvent/` directory contains the corresponding
solvated-ligand inputs. All three ligands use the shared 31-window MDP set
stored directly under `MDP/bulk_solvent/`; their topologies use the same
alchemical molecule name, `LIG`.
The bulk input and provenance status are described in
[`systems/bulk_solvent/README.md`](systems/bulk_solvent/README.md).

## Charge-preserving alchemical setup

The protein-bound calculations use the supplied
[`modify_GRO_file_ABFEP_charge_setup.py`](scripts/modify_GRO_file_ABFEP_charge_setup.py)
workflow. In charged-ligand systems,
the last two solvent molecules are replaced by an additional sodium ion and
an alchemical `ION` molecule. `ions-cl.itp` defines the alchemical molecule as
chloride in state A and a TIP3P water molecule in state B. The additional
sodium ion remains present, preserving the total endpoint charge during the
ligand transformation. The source protocol includes `posres_ION.itp` as an
intentionally empty compatibility placeholder, so no additional position
restraint is active for `ION`. The required `IM` and `du` atom types are
included in the released force-field topology.

The preparation script and usage notes are in
[`scripts/`](scripts/README.md). The matching topology components are already
included in each applicable released system.

## Boresch restraint

The six-coordinate Boresch restraint was generated with
[`boresch_restraint.py`](scripts/boresch_restraint.py). The analytic standard-state correction is
reported for every calculation in its `boresch_correction.txt` file.

## Free-energy convention

For every protein-bound calculation, the released summary uses

```text
ΔG_bind = −ΔG_protein + ΔG_solvent + ΔG_Boresch
```

All values in [`abfep_results.xls`](analysis/abfep_results.xls) are in
kcal/mol. Rows whose `system_name` is `system_<NNN>` contain individual
GROMACS calculations. A `SUMMARY` row follows each condition and contains
Excel formulas for the mean, sample standard deviation, and sample count. The
columns are:

- `condition`: D1R–ligand condition and G-protein state;
- `g_protein`: whether the G protein was included;
- `system_name`: anonymous public calculation identifier;
- `ligand`: orthosteric ligand used in the calculation;
- `protein_leg`: 5 ns protein-leg forward estimate;
- `bulk_solvent_leg`: ligand bulk-solvent decoupling estimate;
- `boresch_correction`: analytic restraint correction;
- `binding_dG`: calculated binding free energy;
- `mean_dG`: condition mean on `SUMMARY` rows;
- `stdev`: condition sample standard deviation on `SUMMARY` rows; and
- `sample_n`: number of calculations in the condition.

The raw selected forward estimates and Boresch outputs are retained under
`analysis/protein_bound/`. This permits the workbook values to be checked
against the source analysis outputs without releasing the trajectories. No
analysis notebook is included; the released values and summary formulas are
in `analysis/abfep_results.xls`.

## Calculation-level replicates and statistics

Each anonymous `system_<NNN>` is one selected starting structure used for a
separate ABFEP calculation. These are calculation-level replicates; structures
sampled from the same parent MD trajectory should not be interpreted as
independent MD trajectories. Reported uncertainties across these systems are
sample standard deviations. Between-condition tests reported in the manuscript
are two-sided Welch tests. Values below 0.001 should be reported as `P < 0.001`.

