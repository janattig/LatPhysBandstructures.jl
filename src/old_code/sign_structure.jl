###################################################################################################
#  Extremely manually keyed-in plaquette sign structure required for Hamiltonian
#
# (1) plaquette(uc,b)
#
#       WORKS ONLY FOR (10,3)a !!!
#
#
#    -> carries the signs for bonds between nearest neighbors
#    -> VERY manual/ Brute-force
#    -> The way the signs have been keyed in: Print the list of all bonds in a unitcell of (10,3)a
#    -> Refer to relevant paper for the signs a bond between each pair of nearest-sites carries
#    -> See which pair-of-sites it corresponds to in the bond-list printed  
#    -> Make matrix of signs "plaquette_convention" accordingly
#
#
#
#  (2) plaquetteNNN
#
#      -> Sign structure for bonds between next nearest neighbors
#      -> Works for sign-structure of any unitcell once the plaquette(uc,b) is set-up properly
#      -> Not as brute-force as the plaquette(uc,b) function
#
#######################################################################################################

function plaquette(
        uc :: U,
        b  :: B
    ) where {L,N,S,B<:AbstractBond{L,N},U<:AbstractUnitcell{S,B}}
    
   bond_index =  findall(bd -> bd == b,bonds(uc))
  
  # sign-structure ONLY FOR (10,3)a!!!

  plaquette_convention =   -[-1 1 1 1 -1 -1 1 -1 1 -1 -1 1]
   
    return( plaquette_convention[bond_index])
end


#sign-structure for next-nearest neighbors
function plaquetteNNN(
        uc :: U,
        b  :: B
    ) where {L,N,S,B<:AbstractBond{L,N},U<:AbstractUnitcell{S,B}}
    
    #b is the bond between pair of next-nearest neighbors
    i_from = from(b)
    i_to = to(b)
    z_eff = 0

    # to find the site lying in between the pair of NNN neighbors i_from and i_to

    # list of all bonds in NN-list in which i_from is starting site and i_to is ending site
    from_bonds = findall(bd -> from(bd)==i_from, bonds(uc))
    to_bonds = findall(bd -> to(bd)==i_to, bonds(uc))
    
    for i in 1:length(from_bonds)
        current_frombond = bonds(uc)[from_bonds[i]]
        for j in 1:length(to_bonds)
            current_tobond = bonds(uc)[to_bonds[j]]

            #adding up position vectors
            uc_add = wrap(current_tobond) .+ wrap(current_frombond)
            uc_sub = wrap(current_tobond) .- wrap(current_frombond)

            if(to(current_frombond)==from(current_tobond))
              
              #check if uc_add matches position vector of NNN-pair that we have
              if(uc_add==wrap(b))
                  
                  #extract signs of the two bonds involved
                  z1 = plaquette(uc,current_frombond)
                  z2 = plaquette(uc,current_tobond)

                  #to extract the type of bond (x,y or z)
                  dir_index = findall(bd -> from(bd)== i_from&&to(bd)==i_to, bonds(uc))
                   sign_bond = bonds(uc)[dir_index][1]

                   #enter the info of bond-labels to calculate handedness
                   perm = [label(current_frombond),label(current_tobond),label(sign_bond)]

                   #effective/resultant sign of NNN bond
                   z_eff = levicivita(perm).*z1.*z2

            end
        end
    end
 end
    return z_eff[1]
end