# Bulk-solvent calculations

This directory contains one solvated-ligand input set for each ligand. All
three inputs use the shared 31-window alchemical MDP set stored directly under
`MDP/bulk_solvent/`. Each topology names the decoupled molecule `LIG`.

The extracted ligand-specific 20 ns bulk-solvent values are:

| Ligand | Bulk-solvent decoupling (kcal/mol) |
| --- | ---: |
| Apomorphine | 108.012 |
| Dopamine | 157.602 |
| SKF81297 | 145.491 |

The extracted values are stored under
`analysis/bulk_solvent/<ligand>/bulk_solvent_forward_20ns.txt`.
