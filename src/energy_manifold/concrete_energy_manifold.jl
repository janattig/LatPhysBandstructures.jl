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
            resolution :: Union{Int64, Vector{Int64}} = 100,
            bounds_lower :: Vector{<:Real} = ones(10) .* -2 .*pi,
            bounds_upper :: Vector{<:Real} = ones(10) .*  2 .*pi,
            max_errors :: Integer = 1000,
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

    # the current search index
    index = 1
    # the current amount of errors (no k point found)
    errors_total = 0

    # search until there are enough points
    while index <= n_points

        # throw an error if too many consecutive errors occured
        if errors_total > max_errors
            if index == 1
                error("Too many errors occured consecutively without finding points --> maybe no points with energy = $(energy_cut)?")
            else
                error("Too many errors occured consecutively, only $(index-1) points of $(n_points) found.")
            end
        end

        # find a suitable starting point for the Newton algorithm
        k = rand(D)
        @simd for j in 1:D
            k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end

        # start with the initial energy
        evs = eigvals(matrixAtK(h,k))
        e0  = minimum(abs.(evs.-energy_cut))
        # iterate i over 100 newton steps (maximum)
        for i in 1:100
            # check if the energy is already converged
            if e0 < precision_energy
                # save the k vector
                k_values[index] = k
                # increment the index
                index = index+1
                # reset the error count
                errors_total = -1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the energy when going along primary directions
            H_eps = zeros(D)
            for j in 1:D
                k_grad = deepcopy(k)
                k_grad[j] += diff_step_size_k
                evs = eigvals(matrixAtK(h,k_grad))
                H_eps[j] = minimum(abs.(evs.-energy_cut))
            end
            # the gradient of the energy
            dH = (H_eps .- H_0) ./ diff_step_size_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                errors_total += 1
                break
            end
            # increment k
            dk = dH .* (H_0 / dHdH)
            k .-= dk
            # calculate a new energy
            evs = eigvals(matrixAtK(h,k))
            e0  = minimum(abs.(evs.-energy_cut))
        end
        # rais errors by 1
        errors_total += 1
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        ruc = getReciprocalUnitcell(uc)
        for i in 1:n_points
            shiftToFirstBZ!(ruc, k_values[i])
        end
    end

    # set the list
    em.k_points = k_values

    # return nothing
    return nothing
end
