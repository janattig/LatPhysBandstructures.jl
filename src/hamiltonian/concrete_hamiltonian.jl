################################################################################
#
#   DEFINITION OF CONCRETE HAMILTONIAN <: AbstractHamiltonian
#   - label type L (of bond labels)
#   - unitcell type U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
#   - bond hamiltonian type HB<:AbstractBondHamiltonian{L,NS} where {NS}
#
################################################################################





################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
struct Hamiltonian{L,UC,HB} <: AbstractHamiltonian{L,UC,HB}

    # the unitcell
    unitcell :: UC

    # the bond hamiltonian
    h_bond   :: HB

end

# export the type
export Hamiltonian








################################################################################
#
#	INTERFACING / ACCESSING HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################


# Function to obtain the unitcell
function unitcell(
            h :: Hamiltonian{L,UC,HB}
        ) :: UC where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS}}

    # return the unitcell
    return h.unitcell
end



# Function to obtain the bond hamiltonian
function bondHamiltonian(
            h :: Hamiltonian{L,UC,HB}
        ) :: HB where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS}}

    # return the bond hamiltonian
    return h.h_bond
end




# easy constructor
function Hamiltonian(
            unitcell :: UC,
            h_bond   :: HB
        ) :: Hamiltonian{L,UC,HB} where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS}}

    # return a newly created hamiltonian
    return Hamiltonian{L,UC,HB}(unitcell, h_bond)
end
