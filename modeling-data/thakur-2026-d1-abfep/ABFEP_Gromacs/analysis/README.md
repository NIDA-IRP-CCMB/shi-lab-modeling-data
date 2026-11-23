# Analysis files

`abfep_results.xls` is the single results table. Rows whose `system_name` is
`system_<NNN>` contain one anonymous GROMACS ABFEP starting structure with its
protein, bulk-solvent, Boresch, and final binding-free-energy values. A
`SUMMARY` row follows each condition and contains Excel formulas for the mean,
sample standard deviation, and sample count.

For every GROMACS row, the corresponding directory under
`protein_bound/<condition>/<system_name>/` contains:

- `protein_forward_5ns.txt`: the 5000-ps protein-bound estimate used in the
  released summary; and
- `boresch_correction.txt`: the analytic Boresch correction and parameters.

Each directory under `bulk_solvent/<ligand>/` contains
`bulk_solvent_forward_20ns.txt`, the 20000-ps bulk-solvent estimate used for
that ligand.

The `ligand` column records the bound ligand. All released protein-bound
inputs use the same MDP templates and the common topology molecule name
`LIG`.
System paths follow the fixed layout
`systems/protein_bound/<condition>/<system_name>/solvated.gro.gz` and
`solvated.top`, with condition-level `toppar/` and `index.ndx` files.

No analysis notebook is included. The calculation-level values and summary
formulas required to reproduce the released table are contained in
`abfep_results.xls`.
