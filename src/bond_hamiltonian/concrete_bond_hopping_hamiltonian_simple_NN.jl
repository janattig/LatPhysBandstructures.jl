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
        ) :: Float64 where {N,L,NB,NS}

    # check if it is the NN label
    if label(b) == h.label
        return h.coupling
    else
        return 0.0
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










function getBondHamiltonianType(
        ::Val{:BondHoppingHamiltonianSimpleNN}
    )

    # return the type
    return BondHoppingHamiltonianSimpleNN
end

function saveBondHamiltonian(
        hb :: HB,
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
        ;
        append :: Bool = false
    ) where {L,HB<:BondHoppingHamiltonianSimpleNN{L}}

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
        attrs(group_hb)["type"] = "BondHoppingHamiltonianSimpleNN"
        # save the parameters
        attrs(group_hb)["N"]  = 1
        attrs(group_hb)["L"]  = string(L)
        # save the Float64 coupling
        group_hb["coupling"]  = hb.coupling
        if L <: Number
            group_hb["label"] = hb.label
        else
            group_hb["label"] = string(hb.label)
        end
    end

    # return nothing
    return nothing
end

function loadBondHamiltonian(
        ::Type{HB},
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
    ) where {LI,HB<:Union{BondHoppingHamiltonianSimpleNN{LI},BondHoppingHamiltonianSimpleNN}}

    # read attribute data
    attr_data = h5readattr(fn, group)
    # determine D based on this
    L = Meta.eval(Meta.parse(attr_data["L"]))
    N = attr_data["N"]

    # load coupling
    coupling  = h5read(fn, group*"/coupling")

    # load label
    label     = L(h5read(fn, group*"/label"))

    # return the new bond hamiltonian
    return BondHoppingHamiltonianSimpleNN{L}(coupling, label)
end
