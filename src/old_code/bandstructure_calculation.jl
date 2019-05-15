################################################################################
#
#   BAND STRUCTURE CALCULATION
#
################################################################################
"""
    getBandstructure(
                unitcell::U,
                path::path
             [; resolution::Int64=-1,
                enforce_hermitian::Bool=false ]
            )

Calculates the band struture of a `Unitcell` object along some path given by a `Path` object and returns a `Bandstructure` object.
"""


function getBands(
                unitcell :: U,
                path :: RPA,
                bondInteractionMatrix:: HB,
                segment_resolution::Array{Int64, 1},
                kappa::Float64;
                enforce_hermitian::Bool=false,
                epsilon_degenerate::Float64=1e-6
            ) :: Bandstructure{RPA} where {
                 L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {
                 NS}, RPA<: AbstractReciprocalPath{P} where {
                 P <: AbstractReciprocalPoint{D} where {D}}}

   #get complete information of unitcell and couplings in a single object

    uc_bondInteractionMatrix = newHamiltonian(unitcell,bondInteractionMatrix)


    # build up the list of all bands of all segments
    segments_total    = Array{Array{Float64,1},1}[]


    # get the spin dimension

    spin_dimension = size(bondtermhopping(bondInteractionMatrix, bonds(unitcell)[rand(1:length(bonds(unitcell)))]))[1]

    # iterate over all path segments and push empty lists into the segments list
    for i in 1:length(segment_resolution)

        # build an empty band structure for this segment
        segment     = Array{Float64, 1}[]

        for b in 1:length(unitcell.sites)*spin_dimension
            push!(segment,     zeros(Float64, segment_resolution[i]))
           end

        # push the segment band structure into the complete segment list
        push!(segments_total,    segment)

    end



    # iterate over all segments
    for s in 1:length(segment_resolution)
        # get the grid in between two points
        resolution = segment_resolution[s]
        ss = (s+1)
        # get the local start and end point of the segment
        k1 = path.points[s].point
        k2 = path.points[ss].point
        kpoints  = getPointsOnLine(k1, k2, resolution)

        # calculate all energies
        for i in 1:resolution

            # get the current k
            k = kpoints[i]

            # get the interaction matrix for this k
            matrix = getMatrixAtK(uc_bondInteractionMatrix,k,kappa)
            # diagonalize the matrix
            eigenfactorization = eigen(matrix)
            eigenvalues  = (eigenfactorization.values)
            eigenvectors = eigenfactorization.vectors
            # save all the eigenvalues to their lists
            for b in 1:length(eigenvalues)
                segments_total[s][b][i] = eigenvalues[b]
            end

            ###########

        end
    end

    # generate a new band structure object
    bandstructure = newBandStructure(path, segments_total)

    # return the band structure
    return bandstructure
end
export getBands
