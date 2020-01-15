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









function saveHamiltonian(
        h  :: H,
        fn :: AbstractString,
        group :: AbstractString = "hamiltonian"
        ;
        append :: Bool = false
    ) where {L,UC<:AbstractUnitcell,N,HB<:AbstractBondHamiltonian{L,N}, H<:AbstractHamiltonian{L,UC,HB}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'saveHamiltonian' for Hamiltonian type " * string(H))
end

function loadHamiltonian(
        ::Type{H},
        fn :: AbstractString,
        group :: AbstractString = "hamiltonian"
    ) where {L,UC<:AbstractUnitcell,N,HB<:AbstractBondHamiltonian{L,N}, H<:AbstractHamiltonian{L,UC,HB}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'loadHamiltonian' for Hamiltonian type " * string(H))
end

export loadHamiltonian, saveHamiltonian











################################################################################
#
#	MATRIX CALCULATION
#
################################################################################

# function to compute the matrix for any bond term dimension
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

# function to compute the matrix for bond term dimension 1
function matrixAtK(
            h :: H,
            k :: Vector{<:Real}
        ) :: Matrix{Complex{Float64}} where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,1},H<:AbstractHamiltonian{L,UC,HB}}

    # new matrix
    h_matrix = zeros(Complex{Float64}, dim(h),dim(h))

    # add all bonds to the matrix
    @inbounds @simd for b in bonds(unitcell(h))
        # get the interaction matrix and add it to the general matrix
        h_matrix[from(b), to(b)] += 0.5 * bondterm(bondHamiltonian(h), b) * exp(im*dot(k, vector(b, unitcell(h))))
    end

    # return the matrix
    return h_matrix
end

# export matrix function
export matrixAtK



# function to find a k where a certain energy is hit
function findKAtEnergy(
            h :: H,
            energy_cut :: Real,
            k_start :: Vector{<:Real}
            ;
            max_variation ::Real = 20*pi,
            max_errors :: Integer = 10,
            max_steps  :: Integer = 100,
            precision_energy :: Real = 1e-6,
            diff_step_size_k :: Real = 1e-6,
        ) :: Vector{Float64} where {L,LS,D,S<:AbstractSite{LS,D},N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # search until there are enough points
    for __try__ in 1:max_errors
        # assume k to be k_start
        k = k_start .+ Float64[2*(rand()-0.5)*max_variation*(__try__-1)/(max_errors-1) for i in 1:length(k_start)]
        H_eps = zeros(D)
        # start with the initial energy
        evs = eigvals(matrixAtK(h,k))
        e0  = minimum(abs.(evs.-energy_cut))
        # iterate i over 100 newton steps (maximum)
        for i in 1:max_steps
            # check if the energy is already converged
            if e0 < precision_energy
                #println("$(__try__-1) failed attempts, $(i) steps needed in attempt $(__try__)")
                # just return the k vector
                return k
            end
            # the current energy
            H_0 = e0
            # the energy when going along primary directions
            @simd for j in 1:D
                k_grad = deepcopy(k)
                @inbounds k_grad[j] += diff_step_size_k
                evs = eigvals(matrixAtK(h,k_grad))
                @inbounds H_eps[j] = minimum(abs.(evs.-energy_cut))
            end
            # the gradient of the energy
            dH = (H_eps .- H_0) ./ diff_step_size_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH .* (H_0 / dHdH)
            k .-= dk
            # calculate a new energy
            evs = eigvals(matrixAtK(h,k))
            e0  = minimum(abs.(evs.-energy_cut))
        end
    end

    # raise error
    error("Too many errors occured consecutively without finding points --> maybe no points with energy = $(energy_cut)?")
end

# export finding function
export findKAtEnergy



# dimension calculation
function dim(
            h :: H
        ) :: Int64 where {L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS},H<:AbstractHamiltonian{L,UC,HB}}

    # return the size of the hamiltonian
    return NS*numSites(unitcell(h))
end

# export dimension function
export dim
