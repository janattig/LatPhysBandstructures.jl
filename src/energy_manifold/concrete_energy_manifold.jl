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
export EnergyManifold





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
            bounds_lower :: Vector{<:Real} = ones(10) .* -2 .*pi,
            bounds_upper :: Vector{<:Real} = ones(10) .*  2 .*pi,
            max_errors :: Integer = 1000,
            max_steps  :: Integer = 100,
            precision_energy :: Real = 1e-6,
            diff_step_size_k :: Real = 1e-6,
            refold_to_first_BZ :: Bool = true
        ) where {L,LS,D,S<:AbstractSite{LS,D},B,UC<:AbstractUnitcell{S,B},HB,H<:AbstractHamiltonian{L,UC,HB}}

    # list of k values that contribute to the BZ
    k_values = Vector{Vector{Float64}}(undef, n_points)

    # obtain the hamiltonian
    h  = hamiltonian(em)
    uc = unitcell(h)

    # obtain the energy for the cut
    energy_cut = energy(em)


    # FIND ALL K POINTS
    @simd for i in 1:n_points
        # build starting point
        k = rand(D)
        @simd for j in 1:D
            @inbounds k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end
        # find the k value
        @inbounds k_values[i] = findKAtEnergy(h, energy_cut, k)
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        ruc = getReciprocalUnitcell(uc)
        @simd for i in 1:n_points
            @inbounds shiftToFirstBZ!(ruc, k_values[i])
        end
    end

    # set the list
    em.k_points = k_values

    # return nothing
    return nothing
end




# convenience function for obtaining the manifold
function getEnergyManifold(
            h     :: H,
            e_cut :: Real,
            n     :: Int64 = 1000
            ;
            recalculate :: Bool = true
        ) :: EnergyManifold{H} where {L,LS,D,S<:AbstractSite{LS,D},B,UC<:AbstractUnitcell{S,B},HB,H<:AbstractHamiltonian{L,UC,HB}}

    # create a new object
    em = EnergyManifold{H}(h, e_cut, Vector{Float64}[])

    # recalculate the numerical energy values
    if recalculate
        recalculate!(em, n)
    end

    # return the new object
    return em
end

# Convinience function if only unitcell and the bond hamiltonian is known (uses concrete Hamiltonian)
function getEnergyManifold(
        uc    :: UC,
        hb    :: HB,
        e_cut :: Real,
        n     :: Int64 = 1000
        ;
        recalculate :: Bool = true
    ) :: EnergyManifold{Hamiltonian{L,UC,HB}} where {L,N,LS,D,S<:AbstractSite{LS,D},B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS}}

    # create a new object
    em = EnergyManifold{Hamiltonian{L,UC,HB}}(Hamiltonian(uc,hb), e_cut, Vector{Float64}[])

    # recalculate the numerical energy values
    if recalculate
        recalculate!(em, n)
    end

    # return the new object
    return em
end

# Convinience function if only path and unitcell but not even the bond hamiltonian is known (uses concrete Hamiltonian and simple hopping hamiltonian)
function getEnergyManifold(
            uc    :: UC,
            e_cut :: Real,
            n     :: Int64 = 1000
            ;
            recalculate :: Bool = true
        ) :: EnergyManifold{Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}} where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B}}

    # create a new object
    em = EnergyManifold{Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}}(Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}(uc, getHoppingHamiltonianSimple(uc)), e_cut, Vector{Float64}[])

    # recalculate the numerical energy values
    if recalculate
        recalculate!(em, n)
    end

    # return the new object
    return em
end


# export the convenience function
export getEnergyManifold
