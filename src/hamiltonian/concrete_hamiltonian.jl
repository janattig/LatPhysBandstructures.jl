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













function saveHamiltonian(
        h  :: H,
        fn :: AbstractString,
        group :: AbstractString = "hamiltonian"
        ;
        append :: Bool = false
    ) where {L,UC<:AbstractUnitcell,N,HB<:AbstractBondHamiltonian{L,N}, H<:Hamiltonian{L,UC,HB}}

    # determine the mode based on if one wants to append stuff
    if append
        mode = "r+"
    else
        mode = "w"
    end

    # group for bond hamiltonian
    group_hb = group*"/bond_hamiltonian"
    # group for unitcell
    group_uc = group*"/unitcell"

    # open the file in mode
    h5open(fn, mode) do file
        # create the group in which the bonds are saved
        group_h = g_create(file, group)
        # save the parameters
        attrs(group_h)["N"] = N
        attrs(group_h)["L"] = string(L)
        # save the groups
        attrs(group_h)["bond_hamiltonian"] = group_hb
        attrs(group_h)["unitcell"]         = group_uc
    end

    # save bond hamiltonian and unitcell in respective groups (append!)
    saveBondHamiltonian(bondHamiltonian(h), fn, group_hb, append=true)
    saveUnitcell(unitcell(h), fn, group_uc, append=true)

    # return nothing
    return nothing
end

function loadHamiltonian(
        ::Type{H},
        fn :: AbstractString,
        group :: AbstractString = "hamiltonian"
    ) where {L,UC<:AbstractUnitcell,N,HB<:AbstractBondHamiltonian{L,N}, H<:Union{Hamiltonian{L,UC,HB}, Hamiltonian}}

    # read attribute data
    attr_data = h5readattr(fn, group)

    # load bond hamiltonian
    hb = loadBondHamiltonian(fn, attr_data["bond_hamiltonian"])
    # load unitcell
    uc = loadUnitcell(fn, attr_data["unitcell"])

    # return the new hamiltonian
    return Hamiltonian(uc, hb)
end

function loadHamiltonian(
        fn :: AbstractString,
        group :: AbstractString = "hamiltonian"
    )

    return loadHamiltonian(Hamiltonian, fn, group)
end
