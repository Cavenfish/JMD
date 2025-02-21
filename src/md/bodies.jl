

#Atoms in simulation
#Needs mutable to swap masses on the fly
mutable struct Atom <: MyAtoms
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  s::Char
end

function getMols(bdys::Vector{MyAtoms}, rmax; D=3)
  r   = [i.r for i in bdys]

  pts = hcat(r...)

  ret = dbscan(pts[1:D, :], rmax)

  [i.core_indices for i in ret.clusters]
end

function getPairs(bdys::Vector{MyAtoms})

  # Get mols and N
  mols = if length(bdys) <= 3
    getMols(bdys, 1.5, D=length(bdys)-1) 
  else
    getMols(bdys, 1.5)
  end
  N    = size(mols)[1]

  # Make all pairs
  pars = Pair[]
  for i in 1:N
    for j in i+1:N
      push!(pars, Pair(mols[i],mols[j]))
    end
  end

  pars, mols
end