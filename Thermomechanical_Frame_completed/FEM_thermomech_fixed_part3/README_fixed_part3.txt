Fixed thermomechanical scripts for the third review batch.

Files:
1. Tri3_ThermoMech.m
   Added missing Tri3 thermomechanical element routine.
2. Quad4_ThermoMech.m
   Uses standard Q4 node order, adds Tref contribution, adds detJ and load checks.
3. post_Quad4_ThermoMech.m
   Uses the same standard Q4 node order and detJ check.
4. assembleSystem.m
   Skips empty element sets before resolving element metadata and adds basic checks.
5. main_mixed_thermomech_fixed.m
   A runnable validation script for mixed Tri3 + Quad4 free thermal expansion. Plotting calls are guarded by exist(...,'file').
