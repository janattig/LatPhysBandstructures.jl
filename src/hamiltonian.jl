
################################################################################
#
#   DEFINITION OF HAMILTONIAN
#   - label type L (of bond labels)
#   - unitcell type U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
#   - bond hamiltonian type HB<:AbstractBondHamiltonian{L,NS} where {NS}
#
#
#    NOTE: Hamiltonian defined here is in the level of Majoranas.
#
#
#    NOTE ABOUT GAUGE: Refer to sign_structure.jl
################################################################################




struct Hamiltonian{
        L,
        U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
        HB<:AbstractBondHamiltonian{L,NS} where {NS}
    }

    # the unitcell
    unitcell :: U

    # the bond Hamiltonian
    h_bond   :: HB

    # ...

end


# function to create such an object
function newHamiltonian(
            unitcell :: U,
            bond_hamiltonian :: HB
        ) :: Hamiltonian{L,U,HB} where {
             L,U<:AbstractUnitcell{S,B} where {
             S, B<:AbstractBond{L,N} where {N}},HB<:AbstractBondHamiltonian{L,NS} where {NS}}

    # return such an object
    return Hamiltonian{L,U,HB}(unitcell, bond_hamiltonian)
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           


# delta function for unitcell and bond to compute the displacement vector of the bond
function delta(
        uc :: U,
        b  :: B
    ) :: Vector{Float64} where {L,N,S,B<:AbstractBond{L,N},U<:AbstractUnitcell{S,B}}

    # build up the difference vector within the unitcell
    d = point(site(uc,to(b))) .- point(site(uc,from(b))) 
    # add all the bravais lattice vectors
    for i in 1:N
        d .+= wrap(b)[i] .* latticeVectors(uc)[i]
    end
    # return the vector
    return d
end

# function to obtain the Hamiltonian matrix at some k (nearest+next-nearest neighbor contribution)
function getMatrixAtK(
        hamiltonian :: Hamiltonian{L,U,HB},
        k :: Vector{<:Real},
        kappa :: Float64
    ) :: Matrix{Complex} where {L,NS,U,HB<:AbstractBondHamiltonian{L,NS}} 

    # new matrix
  
    h = zeros(Complex, NS*numSites(hamiltonian.unitcell),NS*numSites(hamiltonian.unitcell))
    adj_delta_r = 0+0im
    all_bonds = bonds(hamiltonian.unitcell)
    # add all bonds to the matrix
    for b in bonds(hamiltonian.unitcell)
        # get the bond indices
        i_from  = from(b)
        i_to    = to(b)
        delta_r = delta(hamiltonian.unitcell, b)
        
        # get the interaction matrix and add it to the general matrix
    
    h[(i_from-1)*(NS)+1:(i_from)*(NS), (i_to-1)*(NS)+1:(i_to)*(NS)] .+= 0.5 .*im .* plaquette(hamiltonian.unitcell, b) .* bondtermhopping(hamiltonian.h_bond, b).* exp(im*dot(k,delta_r))
      
     
  end
    h_NN = getMatrixAtKNN(hamiltonian,k,kappa) #hamiltonian for next-nearest neighbor contribution
    
    # return the matrix
    return (h .+ h_NN)
end



# function determine the energy of a k_vector for a given hamiltonian
function energy(
            hamiltonian :: Hamiltonian{L,U,HB},
            k           :: Vector{Float64},
           kappa :: Float64
        ) :: Float64 where {L,NS,U,HB<:AbstractBondHamiltonian{L,NS}}

    # first compute the matrix at the given k
    M_k = getMatrixAtK(hamiltonian, k,kappa)

    # return the minimal eigenvalue
    return eigmin(M_k)
end

# function to obtain the Hamiltonian matrix at some k (only next-nearest neighbor contribution)
function getMatrixAtKNN(
        hamiltonian :: Hamiltonian{L,U,HB},
        k :: Vector{<:Real},
        kappa :: Float64
    ) :: Matrix{Complex} where {L,NS,U,HB<:AbstractBondHamiltonian{L,NS}} 

    # new matrix
  
    h = zeros(Complex, NS*numSites(hamiltonian.unitcell),NS*numSites(hamiltonian.unitcell))
    
    all_bonds = getBondsNNNByBondDistance(hamiltonian.unitcell)
    # add all bonds to the matrix
    for b in all_bonds
        # get the bond indices
        i_from  = from(b)
        i_to    = to(b)
        delta_r = delta(hamiltonian.unitcell, b)
        
        # get the inteaction matrix and add it to the general matrix
          h[(i_from-1)*(NS)+1:(i_from)*(NS), (i_to-1)*(NS)+1:(i_to)*(NS)] .+= 0.5.* im.* plaquetteNNN(hamiltonian.unitcell, b) .* kappa .* exp(im*dot(k,delta_r))
      
           end
    # return the matrix
    return h
end
