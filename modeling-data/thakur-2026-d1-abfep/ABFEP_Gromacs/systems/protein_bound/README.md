# Protein-bound input systems

Each condition contains anonymous `system_<NNN>` directories. Public names do
not encode the parent trajectory or snapshot time.

Each condition stores identical files only once:

- `toppar/`: the shared topology includes; and
- `index.ndx`: the default shared GROMACS index.

Each system directory contains only its unique inputs:

- `solvated.gro.gz`: authoritative gzip-compressed coordinates; and
- `solvated.top`: the system-specific topology containing its Boresch
  reference distance, angles, and dihedrals.


```bash
gunzip -k solvated.gro.gz
gmx grompp -f ../../../../MDP/protein_bound/enmin.mdp -c solvated.gro -p solvated.top -n ../index.ndx
```

To generate a full-system PDB after decompression:

```bash
gmx editconf -f solvated.gro -o starting_structure.pdb
```

All released ligand topology molecule types are named `LIG`, so the same MDP
templates under `MDP/protein_bound/` apply to every condition. Anonymous
system names and their calculated values are recorded in
`analysis/abfep_results.xls`.
