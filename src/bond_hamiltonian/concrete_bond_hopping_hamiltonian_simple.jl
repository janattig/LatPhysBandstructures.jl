################################################################################
#
#	ABSTRACT TYPE
#
#   BondHoppingHamiltonianSimple <: AbstractBondHamiltonian{L,1}
#   --> L is the label type of bonds
#   --> N is the dimension the bond term matrix (NxN matrix)
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - generator function
#
################################################################################






################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
mutable struct BondHoppingHamiltonianSimple{L} <: AbstractBondHamiltonian{L,1}

    # coupling
    coupling :: Float64

end

# export the type
export BondHoppingHamiltonianSimple



################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################



# get bond term of Hamiltonian
function bondterm(
            h :: BondHoppingHamiltonianSimple{L},
            b :: AbstractBond{L,NB}
        ) :: Float64 where {N,L,NB,NS}

    # every bond has the same coupling
    return h.coupling
end

# export the function
export bondterm



function getHoppingHamiltonianSimple(
            unitcell :: U
        ) :: BondHoppingHamiltonianSimple{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

    # just create a new hamiltonian
    return BondHoppingHamiltonianSimple{L}(1.0)
end

export getHoppingHamiltonianSimple
