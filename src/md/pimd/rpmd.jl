
#Inputs for RPMD
struct RPMD
  n::UInt8
  
end

#Beads in simulation
mutable struct Bead
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  l::Int64 #link to real particle
end

function initBeads(bdys, nBeads, tj)
  beads  = Bead[]
  frames = 2:length(tj.t) |> collect
  
  for i in 1:length(bdys)
    for j = 1:nBeads
      atom = getFrame(tj, rand(frames))[i]
      push!(beads, Bead(atom.r, atom.v, atom.m, i))
    end
  end

  beads
end
