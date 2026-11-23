# Alchemical MDP files

Protein-bound MDP templates and their generation scripts are stored once under
`protein_bound/`. This shared 64-window protocol is used for apomorphine,
dopamine, and SKF81297. Every released protein-bound topology names its
alchemical ligand molecule `LIG`, matching `couple-moltype = LIG` in the MDP
files. Per-lambda generated MDP files are omitted because they can be
regenerated from the shared templates.

Bulk-solvent MDP templates and generation scripts are stored once under
`bulk_solvent/`. The shared 31-window protocol also uses
`couple-moltype = LIG` for all three ligand inputs.

Do not use the bulk-solvent MDP files for the protein-bound leg; the lambda
schedules and preparation protocols are different.
