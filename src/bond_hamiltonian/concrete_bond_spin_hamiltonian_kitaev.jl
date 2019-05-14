################################################################################
#
#	ABSTRACT TYPE
#
#   BondSpinHamiltonianKitaev <: AbstractBondHamiltonian{L,3}
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
mutable struct BondSpinHamiltonianKitaev{L} <: AbstractBondHamiltonian{L,3}

    # couplings
    Jx :: Float64
    Jy :: Float64
    Jz :: Float64

    # bond labels of these couplings
    x_bonds :: Vector{L}
    y_bonds :: Vector{L}
    z_bonds :: Vector{L}

end

# export the type
export BondSpinHamiltonianKitaev



################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################



# get bond term of Hamiltonian
function bondterm(
            h :: BondSpinHamiltonianKitaev{L},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex} where {L,NB}

    # get the bond label
    l = label(b)

    # define an empty coupling matrix
    coupling_matrix = zeros(3,3)

    # check for all neighbors if the label is in one of the categories
    if l in h.x_bonds
        coupling_matrix[1,1] += h.Jx
        return coupling_matrix
    elseif l in h.y_bonds
        coupling_matrix[2,2] += h.Jy
        return coupling_matrix
    elseif l in h.z_bonds
        coupling_matrix[3,3] += h.Jz
        return coupling_matrix
    end

    # if no bond found, return 0
    return coupling_matrix
end

# export the function
export bondterm





# FUNCTIONS TO FIND BOND LABELS

# FALLBACK ANY LABEL TYPE
# if specific function for label type L not found, denote all bonds as nearest neighbors
function getBondLabelsKitaevX(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a x
    return filter(c->occursin('x',lowercase(string(c))), all_couplings)
end
function getBondLabelsKitaevY(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a y
    return filter(c->occursin('y',lowercase(string(c))), all_couplings)
end
function getBondLabelsKitaevZ(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a z
    return filter(c->occursin('z',lowercase(string(c))), all_couplings)
end


# CONSTRUCTION FUNCTION
function getSpinHamiltonianKitaev(
            unitcell :: U
        ) :: BondSpinHamiltonianKitaev{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

    # obtain all couplings
    couplings = unique!(label.(bonds(unitcell)))

    # create default couplings
    Jx = 1.0
    Jy = 1.0
    Jz = 1.0
    # create list of labels
    x_bonds = getBondLabelsKitaevX(couplings)
    y_bonds = getBondLabelsKitaevY(couplings)
    z_bonds = getBondLabelsKitaevZ(couplings)

    # create and return a new object
    return BondSpinHamiltonianKitaev{L}(Jx,Jy,Jz,x_bonds,y_bonds,z_bonds)
end

# export the construction function
export getSpinHamiltonianKitaev









# CUSTOM SHOWING

# all labels
function show(io::IO, h::H) where {L,H<:BondSpinHamiltonianKitaev{L}}
    print(io, "Kitaev (Bond) Hamiltonian for nearest neighbors. Couplings strength and accepted labels are:")
    print(io, "\n--> Jx="*string(h.Jx)*", labels="*string(h.x_bonds))
    print(io, "\n--> Jy="*string(h.Jy)*", labels="*string(h.y_bonds))
    print(io, "\n--> Jz="*string(h.Jz)*", labels="*string(h.z_bonds))
end
