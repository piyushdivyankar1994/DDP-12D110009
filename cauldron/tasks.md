# Faster Search
  18 July 2016,

  Currently the searching is done using binary search. But since data is well ordered it isn't needed to do this. Making a hash function for the data is relativly simple. It is most likely to be linear. See the eam files for details.
  If form of the function is clear then a function may be written for calculating parameters of this function and hence generating hash function in runtime.

# Chemical Potential
  18 July 2016,

  Chemical potential is the amount by which energy changes when you replace one atom in a material. So  for first attempt we try it in a naive way.

  ## Approch

  1. Calculate energy for the entire matrix, at the begining.
  2. at each montecarlo step there is some change in the lattice, and we know the energy associated with this change so we can calculate the chemical potential at the point .. ?
  3. Update at each step.
