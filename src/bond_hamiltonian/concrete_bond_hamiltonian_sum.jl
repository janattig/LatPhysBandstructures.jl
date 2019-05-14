################################################################################
#
#	ABSTRACT TYPE
#
#   BondHamiltonianSum <: AbstractBondHamiltonian{L,N}
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
mutable struct BondHamiltonianSum{L,N,H1<:AbstractBondHamiltonian{L,N},H2<:AbstractBondHamiltonian{L,N}} <: AbstractBondHamiltonian{L,N}

    # Hamiltonian 1
    h1 :: H1

    # Hamiltonian 2
    h2 :: H2

end

# export the type
export BondHamiltonianSum



################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################



# get bond term of Hamiltonian
function bondterm(
            h :: BondHamiltonianSum{L,N,H1,H2},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex} where {L,N,NB,H1,H2}

    # return the sum of the two sub-hamiltonians
    return bondterm(h.h1,b) .+ bondterm(h.h2,b)
end

# export the function
export bondterm






# Overwrite the base.+ method

# FALLBACK 1 - no agreeing (spin) dimension
function Base.:+(
            h1 :: H1,
            h2 :: H2
        ) :: BondHamiltonianSum{L1,N1,H1,H2} where {L1,L2,N1,N2,H1<:AbstractBondHamiltonian{L1,N1},H2<:AbstractBondHamiltonian{L2,N2}}

    # throw an error
    error("Hamiltonians 1 and 2 don't agree in (spin) dimension: N1="*string(N1)*", N2="*string(N2))
end

# FALLBACK 2 - no agreeing labels
function Base.:+(
            h1 :: H1,
            h2 :: H2
        ) :: BondHamiltonianSum{L1,N,H1,H2} where {L1,L2,N,H1<:AbstractBondHamiltonian{L1,N},H2<:AbstractBondHamiltonian{L2,N}}

    # throw an error
    error("Hamiltonians 1 and 2 don't work on the same label type: L1="*string(L1)*", L2="*string(L2))
end


# Fitting case of same L and N
function Base.:+(
            h1 :: H1,
            h2 :: H2
        ) :: BondHamiltonianSum{L,N,H1,H2} where {L,N,H1<:AbstractBondHamiltonian{L,N},H2<:AbstractBondHamiltonian{L,N}}

    # return a new bond hamiltonian
    return BondHamiltonianSum{L,N,H1,H2}(h1,h2)
end







# CUSTOM SHOWING

# all labels
function show(io::IO, h::H) where {L,N,H<:BondHamiltonianSum{L,N}}
    print(io, "Sum (Bond) Hamiltonian for two Hamiltonians which are:\n1) ")
    show(io, h.h1)
    print("\n2) ")
    show(io, h.h2)
end
