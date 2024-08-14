
#Inputs for RPMD
struct RPMD
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  nBeads::UInt16
  m::Vector
  potVars::PotVars
  EoM::Function
end

#Beads in simulation
mutable struct Bead
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  s::Char
  f::Int64 #family
end

function initBeads(bdys, nBeads, tj)
  beads  = Bead[]
  frames = 2:length(tj.t) |> collect
  
  for i in 1:nBeads
    frame = getFrame(tj, rand(frames))
    for j = 1:length(bdys)
      atom = frame[j]
      push!(beads, Bead(atom.r, atom.v, atom.m, atom.s, i))
    end
  end

  beads
end

function runRPMD_NVE(EoM, tspan, dt, bdys, nBeads, tj; kwargs...)
  beads = initBeads(bdys, nBeads, tj)
  
  pos   = [SVector{3}(i.r) for i in beads]
  vel   = [SVector{3}(i.v) for i in beads]
  mas   = [i.m for i in beads] 

  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  simu       = RPMD(beads, pars, mols, [], nBeads, mas, potVars, EoM)

  prob  = SecondOrderODEProblem(RPMD, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  solu
end

function RPMD(dv, v, u, p, t)

  E  = 0.0
  F  = zero(u)
  y0 = [j for i in u for j in i]
  G  = zero(y0)

  b = length(y0)
  a = div(b, p.nBeads)
  for i = 1:a:b
    @views p.EoM(true, G[i:i+(a-1)], y0[i:i+(a-1)], p)
  end

  for i in 1:3:b
    j     = div(i+2, 3)
    F[j] += -G[i:i+2]
  end

  #Sping forces go here
  
  dv .= F ./ p.m

end