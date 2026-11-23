# D1R molecular dynamics simulation inputs

Each receptor condition is self-contained. For example, `D1_d_Gs/` contains its
starting input files, minimization MDP, equilibration MDPs and production MDPs. Raw trajectory files
are not included.

Production protocols are numbered continuously within each condition as
`prod_0.mdp`, `prod_1.mdp`, and so on. These public names do not expose the
original internal trajectory identifiers. If a branch used a separate
continuation protocol, it is provided as `prod_<N>_cont.mdp`.

The original `gen-vel` and `gen-seed` settings are preserved directly in the
equilibration and production MDP files. A seed of `-1` instructs GROMACS to
select a random seed at run time.

