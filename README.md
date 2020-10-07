# cg_bond_calculator
Calculates bonded interactions from atomistic-to-CG mapped trajectories


```python3
import mdtraj as md
import unyt as u
from cg_bond_calculator.bond_calculator import BondCalculator

traj = md.load("trajectory.dcd", top="topology.hoomdxml")
calc = BondCalculator(traj, T=305)

# Bond and angle parameters will be stored here:
calc.bond_params
calc.angle_params
```
