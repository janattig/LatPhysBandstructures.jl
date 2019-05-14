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



# recalculate the band structure (energy values)
function recalculate!(
            bs :: BS
        ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'recalculate!' for bandstructure type " * string(typeof(bs)))
end

# export the function
export recalculate!
