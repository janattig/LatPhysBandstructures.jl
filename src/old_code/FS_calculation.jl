# Function for calculating Fermi surface in 2D 
function getFermiSurface2D(
            unitcell:: U,
            bondInteractionMatrix:: HB,
            N_points::Int64,
            kappa::Float64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(2),
            bounds_upper::Array{Float64,1}=2*pi.*ones(2),
            refold_to_first_BZ::Bool=true
        ) where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}}


    
    # list of k values that contribute to the BZ
    k_values = zeros(Float64, N_points, 2)
    
    LTbondInteractionMatrix = newHamiltonian(unitcell,bondInteractionMatrix)

    # funtion of energy
    function energyfunction(kvector)
        eigenvalues = eigvals(getMatrixAtK(LTbondInteractionMatrix,kvector,kappa)) .- fermi_energy
        eigenvalues = eigenvalues .* eigenvalues
        return minimum(eigenvalues)
    end

    # the current search index
    index = 1
    # search until there are enough points
    while index <= N_points

        # find a suitable starting point for the Newton algorithm
        k = Float64[rand(), rand()]
        for j in 1:length(k)
            k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end

        # start with the initial energy
        e0 = energyfunction(k)
        # iterate i over 100 newton steps (maximum)
        for i in 1:100
            # check if the energy is already converged
            if e0 < epsilon
                # save the k vector
                k_values[index,:] = k
                # increment the index
                index = index+1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the gradient of the energy
            H_eps =
            [
                energyfunction(k .+ [epsilon_k, 0]),
                energyfunction(k .+ [0, epsilon_k])
            ]
            dH = (H_eps - H_0) ./ epsilon_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH * (H_0 / dHdH)
            k -= dk
            # calculate a new energy
            e0 = energyfunction(k)
        end
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        # get the lattice vectors
        a1 = latticeVectors(unitcell)[1]
        a2 = latticeVectors(unitcell)[2]
        # get the reciprocal lattice vectors
        b1 = [a2[2], -a2[1]]
        b2 = [a1[2], -a1[1]]
        # normalize the vectors
        b1 .*= 2*pi/sum(b1.*a1)
        b2 .*= 2*pi/sum(b2.*a2)

        # iterate over all k points
        for index in 1:N_points
            # assumend to be refolded to begin with
            refolded = true
            # iterate while still refolding possible
            while refolded
                # new iteration, not refolded
                refolded = false
                # try all possible refoldings
                for r in [b1, b2, -b1, -b2, b1+b2, -b2-b1, b1-b2, -b2+b1]
                    if sum((k_values[index,:].-r).*(k_values[index,:].-r)) < sum((k_values[index,:]).*(k_values[index,:]))
                        k_values[index,:] .-= r
                        refolded=true
                    end
                end
            end
        end
    end

    # return the list
    return k_values

end





# Function for calculating Fermi surface in 3D
function getFermiSurface3D(
            unitcell:: U,
            bondInteractionMatrix:: HB,
            N_points::Int64,
            kappa::Float64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}= -2 .*pi .*ones(3),
            bounds_upper::Array{Float64,1}=2 .* pi .*ones(3),
            refold_to_first_BZ::Bool=true
        ) where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}}

    # list of k values that contribute to the BZ
    k_values = zeros(Float64, N_points, 3)
    

    LTbondInteractionMatrix = newHamiltonian(unitcell,bondInteractionMatrix)
    # funtion of energy
    function energyfunction(kvector)
        eigenvalues = eigvals(getMatrixAtK(LTbondInteractionMatrix,kvector,kappa)) .- fermi_energy
        eigenvalues = eigenvalues .* eigenvalues
        return minimum(eigenvalues)
        #return eigmin(getMatrixAtK(LTbondInteractionMatrix,kvector))
    end

    # the current search index
    index = 1
    # search until there are enough points
    while index <= N_points

        # find a suitable starting point for the Newton algorithm
        k = Float64[rand(), rand(), rand()] 
        for j in 1:length(k)
            k[j] = k[j] * bounds_lower[j] + (1-k[j]) * bounds_upper[j]
        end

        # start with the initial energy
        e0 = energyfunction(k)
        # iterate i over 100 newton steps (maximum)
        for i in 1:1000
            # check if the energy is already converged
            if e0 < epsilon
                # save the k vector
                k_values[index,:] = k 
                # increment the index
                index = index+1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the gradient of the energy
            H_eps =
            [
                energyfunction(k .+ [epsilon_k, 0, 0]),
                energyfunction(k .+ [0, epsilon_k, 0]),
                energyfunction(k .+ [0, 0, epsilon_k])
            ]
            dH = (H_eps .- H_0) ./ epsilon_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH * (H_0 / dHdH)
            k -= dk
            # calculate a new energy
            e0 = energyfunction(k)
        end
    end

    # if refolding is enabled, refold all points
    if refold_to_first_BZ
        # get the lattice vectors
        a1 = latticeVectors(unitcell)[1]
        a2 = latticeVectors(unitcell)[2]
        a3 = latticeVectors(unitcell)[3]
        # get the reciprocal lattice vectors
        b1 = cross(a2, a3)
        b2 = cross(a3, a1)
        b3 = cross(a1, a2)
        # normalize the vectors
        b1 .*= 2*pi/sum(b1.*a1) 
        b2 .*= 2*pi/sum(b2.*a2)  
        b3 .*= 2*pi/sum(b3.*a3)
        # put together refolding
        refoldings = Array{Float64,1}[]
        for i in [-1,0,1]
        for j in [-1,0,1]
        for l in [-1,0,1]
            refolding = i*b1 + j*b2 + l*b3
            if sum(abs.(refolding))>1e-8
                push!(refoldings, refolding)
            end
        end
        end
        end
        # iterate over all k points
        for index in 1:N_points
            # assumend to be refolded to begin with
            refolded = true
            # iterate while still refolding possible
            while refolded
                # new iteration, not refolded
                refolded = false
                # try all possible refoldings
                for r in refoldings
                    
                    if sum((k_values[index,:].-r).*(k_values[index,:].-r)) < sum((k_values[index,:]).*(k_values[index,:]))
                        k_values[index,:] .-= r
                        refolded=true
                    end
                end
            end
        end
    end

    # return the list
    return k_values

end





# obtain Fermi surface

"""
    getFermiSurface(
                unitcell::Unitcell,
                N_points::Int64
             [; fermi_energy::Float64=0.0,
                enforce_hermitian::Bool=false,
                epsilon::Float64=1e-10,
                epsilon_k::Float64=1e-10,
                bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
                bounds_upper::Array{Float64,1}=2*pi.*ones(4),
                refold_to_first_BZ::Bool=true ]
            )
calculates `N_points` points which belong to the Fermi surface of the tight-binding model given by
the `Unitcell` object `unitcell`. The Fermi energy can be adjusted by using `fermi_energy=...`.
The points are located using Newton's method with random starting positions. The procedure can be
modified by using some of the optional keywords:
- `enforce_hermitian` determines if the matrix used for energy calculation is made hermitian by construction
- `epsilon` energy threshold for values to be considered close to the Fermi energy
- `epsilon_k` small k increment that is used for calculation of derivatives
- `bounds_upper` / `bounds_lower` the bounds between which the random starting location for the Newton method is searched
After the search has been completed, the values will optionally be refolded into the first brillouin zone by choosing `refold_to_first_BZ=true`.
# Examples
```julia-repl
julia> fermi_surace = getFermiSurface(unitcell, 100)    # 2D unitcell
100x2 Array{Float64}
...
julia> fermi_surace = getFermiSurface(unitcell, 100)    # 3D unitcell
100x3 Array{Float64}
...
julia> fermi_surace = getFermiSurface(unitcell, 100, fermi_energy=1.0)
100x2 Array{Float64}
...
julia> fermi_surace = getFermiSurface(unitcell, 100, refold_to_first_BZ=false)
100x2 Array{Float64}
...
```
"""


function getFermiSurface(
            unitcell:: U,
            bondInteractionMatrix:: HB,
            N_points::Int64,
            kappa::Float64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(3),
            bounds_upper::Array{Float64,1}=2*pi.*ones(3),
            refold_to_first_BZ::Bool=true
        ) where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}}
    
    
    # check if the unitcell can be put in either method
    if length(sites(unitcell)[1].point) != length(latticeVectors(unitcell))
       println("Unitcell has not the same number of lattice vectors as dimensions")
        return zeros(1,1)
    end
    
    
    # check which function to pass to
    if length(latticeVectors(unitcell)) == 2
        # 2D case
        return getFermiSurface2D(
                unitcell,
                bondInteractionMatrix,
                N_points,
                kappa,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ
            )
   
   elseif length(latticeVectors(unitcell)) == 3
        # 3D case
        return getFermiSurface3D(
                unitcell,
                bondInteractionMatrix,
                N_points,
                kappa,
                fermi_energy=fermi_energy,
                enforce_hermitian=enforce_hermitian,
                epsilon=epsilon,
                epsilon_k=epsilon_k,
                bounds_lower=bounds_lower,
                bounds_upper=bounds_upper,
                refold_to_first_BZ=refold_to_first_BZ
            )
    end
end


export getFermiSurface



