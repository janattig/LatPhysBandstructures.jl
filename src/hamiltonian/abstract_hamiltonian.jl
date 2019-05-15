################################################################################
#
#   DEFINITION OF ABSTRACT HAMILTONIAN
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
abstract type AbstractHamiltonian{
        L,
        UC<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
        HB<:AbstractBondHamiltonian{L,NS} where {NS}
    } end

# export the type
export AbstractHamiltonian

# DOCSTRING
"""
    abstract type AbstractHamiltonian{L,UC,HB}

The abstract Hamiltonian type that describes all quadratic hamiltonian implementations in reciprocal space. It is parametric in three types,
- `L` the label type of the bonds it covers
- `UC` the unitcell that it describes
- `HB` the bond hamiltonian that describes interactions on the bonds


# Examples

... not implemented yet ...

"""
AbstractHamiltonian







################################################################################
#
#	INTERFACING / ACCESSING HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################


# Function to obtain the unitcell
import LatPhysBase.unitcell
function unitcell(
            h :: H
        ) :: UC where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'unitcell' for Hamiltonian type " * string(typeof(h)))
end

# export the function
export unitcell



# Function to obtain the bond hamiltonian
function bondHamiltonian(
            h :: H
        ) :: HB where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'bondHamiltonian' for Hamiltonian type " * string(typeof(h)))
end

# export the function
export bondHamiltonian






################################################################################
#
#	MATRIX CALCULATION
#
################################################################################

# function to compute the matrix
function matrixAtK(
            h :: H,
            k :: Vector{<:Real}
        ) :: Matrix{Complex{Float64}} where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # new matrix
    h_matrix = zeros(Complex{Float64}, dim(h),dim(h))

    # add all bonds to the matrix
    @inbounds @simd for b in bonds(unitcell(h))
        # get the bond indices
        i_from  = from(b) :: Int64
        i_to    = to(b)   :: Int64
        delta_r = vector(b, unitcell(h))
        # get the interaction matrix and add it to the general matrix
        h_matrix[(i_from-1)*(NS)+1:(i_from)*(NS), (i_to-1)*(NS)+1:(i_to)*(NS)] .+= 0.5 .* bondterm(bondHamiltonian(h), b) .* exp(im*dot(k,delta_r))
    end

    # return the matrix
    return h_matrix
end

# export matrix function
export matrixAtK




# dimension calculation
function dim(
            h :: H
        ) :: Int64 where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # return the size of the hamiltonian
    return NS*numSites(unitcell(h))
end

# export dimension function
export dim
