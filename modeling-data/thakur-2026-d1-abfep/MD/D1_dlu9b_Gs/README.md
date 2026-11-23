# D1_dlu9b_Gs

This directory contains the starting inputs and MDP protocols for 9 selected
production trajectories.

- `input/`: the original starting coordinates, topology, index, and
  force-field parameters.
- `MDP/minimization/`: energy-minimization parameters.
- `MDP/equilibration/`: staged equilibration and restraint-release parameters.
- `MDP/production/`: continuously numbered production protocols from
  `prod_0.mdp` through `prod_8.mdp`. The single shared
  `prod_cont.mdp` contains the continuation settings used when a new
  continuation TPR was generated.

The original MDP values, including explicit random seeds, are unchanged.
