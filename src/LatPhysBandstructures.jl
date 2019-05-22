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

# Concrete Hamiltonian: Sum of two Hamiltonians
include("bond_hamiltonian/concrete_bond_hamiltonian_sum.jl")

# Concrete Hamiltonian: Heisenberg
include("bond_hamiltonian/concrete_bond_spin_hamiltonian_heisenberg.jl")
# Concrete Hamiltonian: Kitaev
include("bond_hamiltonian/concrete_bond_spin_hamiltonian_kitaev.jl")

# Concrete Hamiltonian: Simple Hopping Hamiltonian (independent on label of bonds)
include("bond_hamiltonian/concrete_bond_hopping_hamiltonian_simple.jl")
# Concrete Hamiltonian: Simple Hopping Hamiltonian (NN)
include("bond_hamiltonian/concrete_bond_hopping_hamiltonian_simple_NN.jl")




# HAMILTONIANS

# Abstract type definition
include("hamiltonian/abstract_hamiltonian.jl")

# Concrete type definition
include("hamiltonian/concrete_hamiltonian.jl")




# BANDSTRUCTURE

# Abstract type definition
include("bandstructure/abstract_bandstructure.jl")

# Concrete type definition
include("bandstructure/concrete_bandstructure.jl")




# CONSTANT ENERGY MANIFOLD

# Abstract type definition
include("energy_manifold/abstract_energy_manifold.jl")

# Concrete type definition
include("energy_manifold/concrete_energy_manifold.jl")

end # module
