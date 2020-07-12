# 2D_jamming_CPP
Pay attention to the energy minimization step.
Do not change the relative order of each steps!
For example, FIRE must come before the second step of Velocity Verlet algorithm.
Otherwise the system will behave strangely and will not converge.
