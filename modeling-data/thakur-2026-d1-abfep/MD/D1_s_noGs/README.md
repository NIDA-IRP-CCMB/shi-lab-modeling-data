# D1_s_noGs

This directory contains the starting inputs and MDP protocols for 5 selected
production trajectories.

- `input/`: starting coordinates, topology, index, and force-field parameters.
- `MDP/minimization/`: energy-minimization parameters.
- `MDP/equilibration/`: staged equilibration and restraint-release parameters.
- `MDP/production/`: continuously numbered production protocols from
  `prod_0.mdp` through `prod_4.mdp`, plus distinct continuation
  protocols where required.

The original MDP values, including explicit random seeds, are unchanged.
