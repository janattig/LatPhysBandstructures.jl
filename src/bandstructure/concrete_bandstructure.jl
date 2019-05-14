################################################################################
#
#	CONCRETE TYPE
#
#   Bandstructure{path} <: AbstractBandstructure{path}
#   --> path is the path in k-space
#
#   FILE CONTAINS
#       - concrete struct definition
#        - interface implementation
#
################################################################################



struct Bandstructure{P,H} <: AbstractBandstructure{P,H}

    # the path along which the band structure is calcualted
    path :: P
    # the hamiltonian
    h :: H

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
    bands :: Array{Array{Array{Float64, 1}, 1}, 1}

end


# export the type
export Bandstructure



################################################################################
#
#	INTERFACING / ACCESSING BANDSTRUCTURES
#	(functions have to be overwritten by concrete types)
#
################################################################################

# obtaining the hamiltonian
function hamiltonian(
            bs :: Bandstructure{P,H}
        ) :: H where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the field
    return bs.h
end

# export the function
export hamiltonian



# obtaining the path
function path(
            bs :: Bandstructure{P,H}
        ) :: P where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the field
    return bs.path
end

# export the function
export path



# recalculate the band structure (energy values)
function recalculate!(
            bs :: Bandstructure{P,H}
        ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'recalculate!' yet for bandstructure type " * string(typeof(bs)))
end

# export the function
export recalculate!
