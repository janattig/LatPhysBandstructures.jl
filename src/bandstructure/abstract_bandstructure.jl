################################################################################
#
#	ABSTRACT TYPE
#
#   AbstractBandstructure{path}
#   --> path is the path in k-space
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
abstract type AbstractBandstructure{
    P<:AbstractReciprocalPath{RP} where {RP},
    H<:AbstractHamiltonian{L,UC,HB} where {L,UC,HB}
} end

# export the type
export AbstractBandstructure





################################################################################
#
#	INTERFACING / ACCESSING BANDSTRUCTURES
#	(functions have to be overwritten by concrete types)
#
################################################################################

# obtaining the hamiltonian
function hamiltonian(
            bs :: BS
        ) :: H where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'hamiltonian' for bandstructure type " * string(typeof(bs)))
end

# export the function
export hamiltonian



# obtaining the path
function path(
            bs :: BS
        ) :: P where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'path' for bandstructure type " * string(typeof(bs)))
end

# export the function
export path



# obtaining the energy values
# bands for each segment
# bands[i] gives all bands of segment i
# bands[i][j] gives all energy values for band j of segment i
# bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
function energies(
            bs :: BS
        ) :: Vector{<:Vector{<:Vector{<:Real}}} where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'energies' for bandstructure type " * string(typeof(bs)))
end

# export the function
export energies



# recalculate the band structure (energy values)
function recalculate!(
            bs :: BS
            ;
            resolution :: Union{Int64, Vector{Int64}} = 100
        ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'recalculate!' for bandstructure type " * string(typeof(bs)))
end

# export the function
export recalculate!
