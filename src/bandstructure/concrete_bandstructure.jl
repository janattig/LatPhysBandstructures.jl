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



mutable struct Bandstructure{P,H} <: AbstractBandstructure{P,H}

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


# obtaining the energy values
function energies(
            bs :: Bandstructure{P,H}
        ) :: Array{Array{Array{Float64, 1}, 1}, 1} where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # return the energy array
    return bs.bands
end








function saveBandstructure(
        bs :: BS,
        fn :: AbstractString,
        group :: AbstractString = "bandstructure"
        ;
        append :: Bool = false
    ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS<:Bandstructure{P,H}}

    # determine the mode based on if one wants to append stuff
    if append
        mode = "r+"
    else
        mode = "w"
    end

    # group for hamiltonian
    group_h = group*"/hamiltonian"
    # group for path
    group_p = group*"/path"
    # group for energy bands
    group_e = group*"/bands"

    # open the file in mode
    h5open(fn, mode) do file
        # create the group in which the bonds are saved
        group_bs = g_create(file, group)
        # write number of segments
        attrs(group_bs)["segments"]    = length(bs.bands)
        # save the groups
        attrs(group_bs)["hamiltonian"] = group_h
        attrs(group_bs)["path"]        = group_p
        attrs(group_bs)["bands"]       = group_e
    end

    # save hamiltonian
    saveHamiltonian(hamiltonian(bs), fn, group_h, append=true)
    # save reciprocal path
    saveReciprocalPath(path(bs), fn, group_p, append=true)

    # save the individual segments within the bands group
    h5open(fn, "r+") do file
        # create the group in which the bonds are saved
        group_bands = g_create(file, group_e)
        # save the individual segments
        for si in 1:length(bs.bands)
            s = bs.bands[si]
            # reformat segment and save as matrix
            #segmat = [s[j][i] for i in 1:length(s[1]) for j in 1:length(s)]
            segmat = zeros(length(s[1]), length(s))
            for i in 1:length(s[1])
                for j in 1:length(s)
                    segmat[i,j] = s[j][i]
                end
            end
            # save
            group_bands["segment_$(si)"] = segmat
        end
    end

    # return nothing
    return nothing
end

function loadBandstructure(
        ::Type{BS},
        fn :: AbstractString,
        group :: AbstractString = "bandstructure"
    ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS<:Union{Bandstructure, Bandstructure{P,H}}}

    # read attribute data
    attr_data = h5readattr(fn, group)

    # group for hamiltonian
    group_h = attr_data["hamiltonian"]
    # group for path
    group_p = attr_data["path"]
    # group for energy bands
    group_e = attr_data["bands"]

    # load number of segments
    segments = attr_data["segments"]

    # load hamiltonian
    h = loadHamiltonian(fn, group_h)
    # load path
    p = loadReciprocalPath(fn, group_p)

    # load all energy bands
    bds = Vector{Vector{Float64}}[]
    # iterate over all expected segments
    for s in 1:segments
        # load segment as matrix
        segmat = h5read(fn, group_e*"/segment_$(s)")
        # reformat segment and push to list
        segment = Vector{Float64}[
            Float64[segmat[i,b] for i in 1:size(segmat,1)]
            for b in 1:size(segmat,2)
        ]
        # push to the list
        push!(bds, segment)
    end

    # return the new Bandstructure
    bs = getBandstructure(h,p,recalculate=false)
    bs.bands = bds
    return bs
end

function loadBandstructure(
        fn :: AbstractString,
        group :: AbstractString = "bandstructure"
    )

    return loadBandstructure(Bandstructure, fn, group)
end








# recalculate the band structure (energy values)
function recalculate!(
            bs :: Bandstructure{P,H}
            ;
            resolution :: Union{Int64, Vector{Int64}} = 100,
            batch_size :: Integer = -1,
            batch_size_min :: Integer = 1,
            batch_size_max :: Integer = 100
        ) where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # calculate the number of segments
    N_segments = length(path(bs)) - 1
    # parse the segment resolution
    seg_resolution = (zeros(Int64,N_segments).+1).*resolution
    if typeof(resolution)==Vector{Int64}
        if length(resolution) == N_segments
            seg_resolution = resolution
        else
            error("given resolution \"$(resolution)\" has not enough components ($(N_segments) needed)")
        end
    end

    # calculate number of bands
    N_bands = dim(hamiltonian(bs))


    # make a list of all k values
    k_values = Vector{Vector{Vector{Float64}}}(undef, N_segments)
    # make a list for all energy values
    e_values = Vector{Vector{Vector{Float64}}}(undef, N_segments)

    # fill all k values
    @inbounds for i in 1:N_segments
        # fill the k values with points on the path segment
        k_values[i] = getPointsOnLine(path(bs)[i], path(bs)[i+1], seg_resolution[i])
        # fill energy values with zeros
        e_values[i] = Vector{Float64}[
            zeros(Float64, seg_resolution[i]) for j in 1:N_bands
        ]
    end


    # calculate the energy values
    @inbounds for i in 1:N_segments
        # the actual batch size for pmap
        batch_size_used = 1
        if batch_size > 0
            batch_size_used = batch_size
        else
            batch_size_used=min(batch_size_min, max(batch_size_max, round(Int64, length(k_values[i])/(10*length(workers())))))
        end
        # get the energy values
        e_values_band = pmap(
            k -> sort(real.(eigvals(matrixAtK(hamiltonian(bs), k)))),
            k_values[i],
            batch_size=batch_size_used
        )
        # format the energy values
        for j in 1:seg_resolution[i]
            # put the eigenvalues into the bands
            @simd for b in 1:N_bands
                e_values[i][b][j] = e_values_band[j][b]
            end
        end
    end

    # set the energy values within band structure
    bs.bands = e_values

    # return nothing
    return nothing
end

# export the function
export recalculate!









# Convinience function if path and hamiltonian are known
function getBandstructure(
            path :: P,
            h    :: H
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{P,H} where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # create a new object
    bs = Bandstructure{P,H}(path, h, Vector{Vector{Vector{Float64}}}(undef, length(path)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

function getBandstructure(
            h    :: H,
            path :: P
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{P,H} where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}}

    # create a new object
    bs = Bandstructure{P,H}(path, h, Vector{Vector{Vector{Float64}}}(undef, length(path)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

function getBandstructure(
            h      :: H,
            points :: Union{Symbol, Tuple{Symbol, Int64}} ...
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{ReciprocalPath{ReciprocalPoint{D}},H} where {LS,D,S<:AbstractSite{LS,D},B,L,UC<:AbstractUnitcell{S,B},HB,H<:AbstractHamiltonian{L,UC,HB}}

    # create new path
    path = getReciprocalPath(unitcell(h), points...)
    # create a new object
    bs = Bandstructure{ReciprocalPath{ReciprocalPoint{D}},H}(path, h, Vector{Vector{Vector{Float64}}}(undef, length(path)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

# Convinience function if only path, unitcell and the bond hamiltonian is known (uses concrete Hamiltonian)
function getBandstructure(
            uc   :: UC,
            path :: P,
            h    :: HB
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{P,Hamiltonian{L,UC,HB}} where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B},NS,HB<:AbstractBondHamiltonian{L,NS}}

    # create a new object
    bs = Bandstructure{P,Hamiltonian{L,UC,HB}}(path, Hamiltonian{L,UC,HB}(uc, h), Vector{Vector{Vector{Float64}}}(undef, length(path)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

# Convinience function if only path and unitcell but not even the bond hamiltonian is known (uses concrete Hamiltonian and simple hopping hamiltonian)
function getBandstructure(
            uc   :: UC,
            path :: P
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{P,Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}} where {D,RP<:AbstractReciprocalPoint{D}, P<:AbstractReciprocalPath{RP}, L,S,N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B}}

    # create a new object
    bs = Bandstructure{P,Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}}(path, Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}(uc, getHoppingHamiltonianSimple(uc)), Vector{Vector{Vector{Float64}}}(undef, length(path)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

# Convinience function if only unitcell and points are known
function getBandstructure(
            uc     :: UC,
            points :: Union{Symbol, Tuple{Symbol, Int64}} ...
            ;
            recalculate :: Bool = true
        ) :: Bandstructure{ReciprocalPath{ReciprocalPoint{D}},Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}} where {LS,D,L,S<:AbstractSite{LS,D},N,B<:AbstractBond{L,N},UC<:AbstractUnitcell{S,B}}

    # create a new object
    bs = Bandstructure{ReciprocalPath{ReciprocalPoint{D}},Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}}(getReciprocalPath(uc, points...), Hamiltonian{L,UC,BondHoppingHamiltonianSimple{L}}(uc, getHoppingHamiltonianSimple(uc)), Vector{Vector{Vector{Float64}}}(undef, length(points)-1))

    # recalculate the numerical energy values
    if recalculate
        recalculate!(bs)
    end

    # return the new object
    return bs
end

# export the convenience function
export getBandstructure
