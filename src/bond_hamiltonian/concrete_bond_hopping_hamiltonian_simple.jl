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











function getBondHamiltonianType(
        ::Val{:BondHoppingHamiltonianSimple}
    )

    # return the type
    return BondHoppingHamiltonianSimple
end

function saveBondHamiltonian(
        hb :: HB,
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
        ;
        append :: Bool = false
    ) where {L,HB<:BondHoppingHamiltonianSimple{L}}

    # determine the mode based on if one wants to append stuff
    if append
        mode = "r+"
    else
        mode = "w"
    end

    # open the file in mode
    h5open(fn, mode) do file
        # create the group in which the bonds are saved
        group_hb = g_create(file, group)
        # save the type identifier
        attrs(group_hb)["type"] = "BondHoppingHamiltonianSimple"
        # save the parameters
        attrs(group_hb)["N"] = 1
        attrs(group_hb)["L"] = string(L)
        # save the Float64 coupling
        group_hb["coupling"] = hb.coupling
    end

    # return nothing
    return nothing
end

function loadBondHamiltonian(
        ::Type{HB},
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
    ) where {LI,HB<:Union{BondHoppingHamiltonianSimple{LI},BondHoppingHamiltonianSimple}}

    # read attribute data
    attr_data = h5readattr(fn, group)
    # determine D based on this
    L = Meta.eval(Meta.parse(attr_data["L"]))
    N = attr_data["N"]

    # load all remaining data
    coupling = h5read(fn, group*"/coupling")

    # return the new bond hamiltonian
    return BondHoppingHamiltonianSimple{L}(coupling)
end
