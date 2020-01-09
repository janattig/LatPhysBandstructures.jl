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












function getBondHamiltonianType(
        ::Val{:BondHoppingHamiltonianDict}
    )

    # return the type
    return BondHoppingHamiltonianDict
end

function saveBondHamiltonian(
        hb :: HB,
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
        ;
        append :: Bool = false
    ) where {L,HB<:BondHoppingHamiltonianDict{L}}

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
        attrs(group_hb)["type"] = "BondHoppingHamiltonianDict"
        # save the parameters
        attrs(group_hb)["N"]    = 1
        attrs(group_hb)["L"]    = string(L)
        # save the labels
        if L <: Number
            group_hb["labels"]  = [p[1] for p in hb.couplings]
        else
            group_hb["labels"]  = [string(p[1]) for p in hb.couplings]
        end
        # save the Float64 couplings
        group_hb["couplings"]   = [p[2] for p in hb.couplings]
    end

    # return nothing
    return nothing
end

function loadBondHamiltonian(
        ::Type{HB},
        fn :: AbstractString,
        group :: AbstractString = "bond_hamiltonian"
    ) where {LI,HB<:Union{BondHoppingHamiltonianDict{LI},BondHoppingHamiltonianDict}}

    # read attribute data
    attr_data = h5readattr(fn, group)
    # determine D based on this
    L = Meta.eval(Meta.parse(attr_data["L"]))
    N = attr_data["N"]

    # load coupling
    couplings  = h5read(fn, group*"/couplings")

    # load label
    labels     = L.(h5read(fn, group*"/labels"))

    # return the new bond hamiltonian
    return BondHoppingHamiltonianDict{L}(Dict([(labels[i], couplings[i]) for i in 1:length(labels)]))
end
