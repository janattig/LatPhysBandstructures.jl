################################################################################
#
#	ABSTRACT TYPE
#
#   BondHoppingHamiltonianSimpleNN <: AbstractBondHamiltonian{L,1}
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
mutable struct BondHoppingHamiltonianSimpleNN{L} <: AbstractBondHamiltonian{L,1}

    # coupling
    coupling :: Float64

    # bond label that identifies NN
    label :: L

end

# export the type
export BondHoppingHamiltonianSimpleNN



################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################



# get bond term of Hamiltonian
function bondterm(
            h :: BondHoppingHamiltonianSimpleNN{L},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex{Float64}} where {N,L,NB,NS}

    # get the bond label
    l = label(b)

    # check if it is the NN label
    if l == h.label
        return zeros(Complex{Float64},1,1) .+ h.coupling
    else
        return zeros(Complex{Float64},1,1)
    end
end

# export the function
export bondterm



function getHoppingHamiltonianSimpleNN(
            unitcell :: U
        ) :: BondHoppingHamiltonianSimpleNN{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

    # obtain all couplings
    couplings = unique!(label.(bonds(unitcell)))

    # check if the list contains a default label
    if getDefaultLabel(L) in couplings
        return BondHoppingHamiltonianSimpleNN{L}(1.0, getDefaultLabel(L))
    else
        return BondHoppingHamiltonianSimpleNN{L}(1.0, couplings[1])
    end
end

export getHoppingHamiltonianSimpleNN
