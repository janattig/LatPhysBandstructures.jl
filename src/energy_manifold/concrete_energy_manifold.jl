################################################################################
#
#	CONCRETE TYPE
#
#   EnergyManifold{hamiltonian}
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
mutable struct EnergyManifold{H} <: AbstractEnergyManifold{H}

    # the hamiltonian
    h :: H

    # energy of the cut
    E_cut :: Float64

    # the k points
    k_points :: Vector{Vector{Float64}}

end

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
            em :: EnergyManifold{H}
        ) :: H where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the internal field
    return em.h
end




# obtaining the k points of the manifold
function kpoints(
            em :: EnergyManifold{H}
        ) :: Vector{Vector{Float64}} where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the internal field
    return em.k_points
end



# obtaining the energy at which the manifold cuts
function energy(
            em :: EnergyManifold{H}
        ) :: Float64 where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the internal field
    return em.E_cut
end




# recalculate the energy manifold (k points)
function recalculate!(
            em       :: EnergyManifold{H},
            n_points :: Integer
            ;
            resolution :: Union{Int64, Vector{Int64}} = 100
        ) where {L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # ...

end
