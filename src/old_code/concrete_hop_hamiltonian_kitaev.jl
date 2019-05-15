################################################################################
#
#	ABSTRACT TYPE
#
#   BondHopHamiltonianKitaev <: AbstractBondHamiltonian{L,1}
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
#   ABSTRACT TYPE DEFINITION OF HOP-HAMILTONIAN
#
################################################################################

mutable struct BondHopHamiltonianKitaev{L} <: AbstractBondHamiltonian{L,1}  # L is the label


    # couplings
    Jx :: Float64
    Jy :: Float64
    Jz :: Float64


    # bond labels of these couplings
    x_bonds :: Vector{L}
    y_bonds :: Vector{L}
    z_bonds :: Vector{L}

end

#export the type
export BondHopHamiltonianKitaev


#get the bond term of the hamiltonian

function bondtermhopping(
            h :: BondHopHamiltonianKitaev{L},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex} where {L,NB}

    # get the bond label
    l = label(b)

    # define an empty coupling matrix
    coupling_matrix = zeros(1,1)

    # check for all neighbors if the label is in one of the categories
    if l in h.x_bonds
        coupling_matrix[1,1] = h.Jx
        return coupling_matrix
    elseif l in h.y_bonds
        coupling_matrix[1,1] = h.Jy
        return coupling_matrix
    elseif l in h.z_bonds
        coupling_matrix[1,1] = h.Jz
        return coupling_matrix
    end

    # if no bond found, return 0
    return coupling_matrix
end


# FUNCTIONS TO FIND BOND LABELS
function getBondLabelsKitaevX(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a x (x-bond corresponds to label=1)
   
    return filter(c->c==1, all_couplings)
end
function getBondLabelsKitaevY(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a y (y-bond corresponds to label=2)
  
    return filter(c->c==2, all_couplings)
end
function getBondLabelsKitaevZ(
        all_couplings :: Vector{L}
    ) :: Vector{L} where {L}

    # all couplings that somehow contain a z (z-bond corresponds to label=3)
    
    return filter(c->c==3, all_couplings)
end

#construction function
function getHopHamiltonianKitaev(
            unitcell :: U
        ) :: BondHopHamiltonianKitaev{L} where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B}}

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
    return BondHopHamiltonianKitaev{L}(Jx,Jy,Jz,x_bonds,y_bonds,z_bonds)
end

# export the construction function
export getHopHamiltonianKitaev

