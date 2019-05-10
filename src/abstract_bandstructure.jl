

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




abstract type AbstractBandstructure{path<:AbstractReciprocalPath{P} where {P <: AbstractReciprocalPoint{D} where {D}}}
end

# export the type
export AbstractBandstructure

################################################################################
#
#	INTERFACING / ACCESSING BANDSTRUCTURES
#	(functions have to be overwritten by concrete types)
#
################################################################################

function newBandStructure(
           p :: path,
           energy_bands :: Array{Array{Array{Float64, 1}, 1}, 1}
        ) :: Bandstructure{path} where {path<: AbstractReciprocalPath{P} where {P <: AbstractReciprocalPoint{D} where {D}}}

 
    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newBandStructure' for concrete bandstructure type for the specified path in " * string(D) * " spatial dimensions")
end

# export the interface
export newBandStructure

