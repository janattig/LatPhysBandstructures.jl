################################################################################
#
#	ABSTRACT TYPE
#
#   BondHamiltonian{L,N}
#   --> L is the label type of bonds
#   --> N is the dimension the bond term matrix (NxN matrix)
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#
################################################################################






################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractBondHamiltonian{L,N} end

# export the type
export AbstractBondHamiltonian


# DOCSTRING
"""
    abstract type AbstractBondHamiltonian{L,N}

The abstract bond Hamiltonian type that describes all hamiltonian implementations on bonds. It is parametric in two types,
- `L` the label type of the bonds it covers
- `N` an integer number specifying the number of on-site degrees of freedom (e.g. number of spin components)

The abstract type is equipped with a single interface function `bondterm` which
returns the Hamiltonian matrix describing the Hamiltonian on that bond.
The function throws an adequate error when called but not implemented for concrete types.
See specific documentation for the function itself.


# Examples

```julia-REPL
julia> AbstractBondHamiltonian{Int64,3} <: AbstractBondHamiltonian{L,N} where{L,N}
true
```

"""
AbstractBondHamiltonian




################################################################################
#
#	INTERFACING / ACCESSING BOND HAMILTONIANS
#	(functions have to be overwritten by concrete types)
#
################################################################################




# Fallback if wrong label type is passed (DOES NOT HAVE TO BE OVERWRITTEN)
function bondterm(
            h :: AbstractBondHamiltonian{LH,NH},
            b :: AbstractBond{LB,NB}
        ) :: Matrix{Complex} where {LB,NB,LH,NH}

    # print an error because implementation for concrete type is missing
    error("Passed a bond with different label type to a bond Hamiltonian\nLabel type of bond: "*string(LB)*", label type of bond Hamiltonian: "*string(LH))
end

# get bond term of Hamiltonian
function bondterm(
            h :: AbstractBondHamiltonian{L,NH},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex} where {L,NB,NH}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'bondterm' for bond Hamiltonian type " * string(typeof(h)))
end

export bondterm



# DOCSTRING
"""
    function bondterm(
            h :: AbstractBondHamiltonian{L,NH},
            b :: AbstractBond{L,NB}
        ) :: Matrix{Complex} where {L,NB,NH}

The abstract bond Hamiltonian interface function. Returns a matrix with complex
entries for a passed `AbstractBond` object `b`.

E.g. for Hamiltonians of the form \$\\sum_{\\langle ij\\rangle} S_{i} J_{ij} S_{j} \\dots\$, this bond hamiltonian matrix would correspond
to the matrix \$J_{ij}\$.


# Examples

```julia-REPL
julia> bondterm(hamiltonian, b)
3×3 Array{Complex,2}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im
```

"""
bondterm




# fancy syntax of passing the bond to the object directly
function (hb :: Type{H})(b :: AbstractBond{LB,NB}) :: Matrix{Complex} where {LB,NB,LH,NH,H<:AbstractBondHamiltonian{LH,NH}}

    # return the bond term function
    return bondterm(hb, b)
end
