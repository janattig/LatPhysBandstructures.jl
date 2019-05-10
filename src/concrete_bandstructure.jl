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



struct Bandstructure{path} <: AbstractBandstructure{path} 
        

    # the path along which the band structure is calcualted
    path  ::  path

    # bands for each segment
    # bands[i] gives all bands of segment i
    # bands[i][j] gives all energy values for band j of segment i
    # bands[i][j][k] gives the energy value at kpoint index k of band j in segment i
     bands  :: Array{Array{Array{Float64, 1}, 1}, 1}

    


end


# export the type
export Bandstructure


################################################################################
#
#	INTERFACING BANDSTRUCTURE/ACCESSING BANDSTRUCTURE
#	
#
################################################################################

#Default constructor interface
# used for creation of new Bandstructures
function newBandStructure(
           p :: path,
           energy_bands :: Array{Array{Array{Float64, 1}, 1}, 1}
        ) :: Bandstructure{path} where {path<: AbstractReciprocalPath{P} where {P <: AbstractReciprocalPoint{D} where {D}}}

     #return the newly created object
    return Bandstructure{path}(p,energy_bands)
end

#export the interface
export newBandStructure