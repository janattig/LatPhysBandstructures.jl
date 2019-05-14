################################################################################
#
#   module LatPhysBandstructures
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LinearAlgebra
#
#   --> Definition of Bond Hamiltonians
#   --> Definition of Hamiltonians
#   --> Definition of Bandstructures
#
################################################################################










module LatPhysBandstructures


# Dependencies
using LatPhysBase
using LatPhysReciprocal
using LinearAlgebra




# BOND HAMILTONIANS

# Abstract type definition
include("bond_hamiltonian/abstract_bond_hamiltonian.jl")

# Concrete Hamiltonian: Heisenberg
include("bond_hamiltonian/concrete_bond_spin_hamiltonian_heisenberg.jl")
# Concrete Hamiltonian: Kitaev
include("bond_hamiltonian/concrete_bond_spin_hamiltonian_kitaev.jl")
# Concrete Hamiltonian: Sum of two Hamiltonians
include("bond_hamiltonian/concrete_bond_hamiltonian_sum.jl")



end # module
