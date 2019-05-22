################################################################################
#
#	ABSTRACT TYPE
#
#   AbstractEnergyManifold{hamiltonian}
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#
#
################################################################################


################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractEnergyManifold{
    H<:AbstractHamiltonian{L,UC,HB} where {L,UC,HB}
} end

# export the type
export AbstractEnergyManifold





################################################################################
#
#	INTERFACING / ACCESSING ENERGY MANIFOLDS
#	(functions have to be overwritten by concrete types)
#
################################################################################

# obtaining the hamiltonian
function hamiltonian(
            aem :: AEM
        ) :: H where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'hamiltonian' for energy manifold type " * string(typeof(aem)))
end

# export the function
export hamiltonian




# obtaining the k points of the manifold
function kpoints(
            aem :: AEM
        ) :: Vector{Vector{Float64}} where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'kpoints' for energy manifold type " * string(typeof(aem)))
end

# export the function
export kpoints



# obtaining the energy at which the manifold cuts
function energy(
            aem :: AEM
        ) :: Float64 where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'energy' for energy manifold type " * string(typeof(aem)))
end

# export the function
export energy


# obtaining the energy at which the manifold cuts
function numPoints(
            aem :: AEM
        ) :: Float64 where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}

    # return length of the point list
    return length(kpoints(aem))
end

# export the function
export numPoints


# recalculate the energy manifold (k points)
function recalculate!(
            aem      :: AEM,
            n_points :: Integer
            ;
            resolution :: Union{Int64, Vector{Int64}} = 100
        ) where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'recalculate!' for energy manifold type " * string(typeof(aem)))
end

# export the function
export recalculate!
