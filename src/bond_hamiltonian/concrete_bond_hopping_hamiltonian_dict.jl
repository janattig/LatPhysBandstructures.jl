################################################################################
#
#	ABSTRACT TYPE
#
#   BondHoppingHamiltonianDict <: AbstractBondHamiltonian{L,1}
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
mutable struct BondHoppingHamiltonianDict{L} <: AbstractBondHamiltonian{L,1}

    # coupling
    couplings :: Dict{L,Float64}

end

# export the type
export BondHoppingHamiltonianDict



################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################



# get bond term of Hamiltonian
function bondterm(
            h :: BondHoppingHamiltonianDict{L},
            b :: AbstractBond{L,NB}
        ) :: Float64 where {N,L,NB,NS}

    # return the matching coupling and 0.0 if not found
    return get(h.couplings, label(b), 0.0)
end

# export the function
export bondterm


# set a couplint within the dictonary of couplings and add it if not present yet
function setCoupling(
            h :: BondHoppingHamiltonianDict{L},
            coupling :: L,
            coupling_strength :: Real
        ) where {L}

    # add the coupling to the dictonary
    h.couplings[coupling] = coupling_strength
end
function setCoupling(
            h :: BondHoppingHamiltonianDict{L},
            bond :: B,
            coupling_strength :: Real
        ) where {L,NB,B<:AbstractBond{L,NB}}

    # add the coupling to the dictonary
    h.couplings[label(bond)] = coupling_strength
end

# export the function
export setCoupling


function getHoppingHamiltonianDict(
            unitcell :: U
        ) :: BondHoppingHamiltonianDict{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

    # just create a new hamiltonian
    return BondHoppingHamiltonianDict{L}(Dict{L,Float64}())
end
function getHoppingHamiltonianDict(
            unitcell  :: U,
            couplings :: Dict{L,<:Real}
        ) :: BondHoppingHamiltonianDict{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

    # just create a new hamiltonian
    hb = BondHoppingHamiltonianDict{L}(Dict{L,Float64}())
    # set all copulings
    for c in keys(couplings)
        setCoupling(hb, c, couplings[c])
    end
    # return the bond hamiltonian
    return hb
end

export getHoppingHamiltonianDict
