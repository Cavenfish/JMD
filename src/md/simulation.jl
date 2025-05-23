struct Simulation
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  temp::Vector
  save::String
  PBC::Vector{Bool}
  NC::Vector{Int32}
  lattice::Matrix
  m::Vector
  potVars::PotVars
  NVT::Bool
  thermostat!
  thermoInps
end

function runMD(EoM, bdys::Vector{MyAtoms}, tspan::Tuple{Float64, Float64},
               dt::Float64; save="full", thermostat=nothing, thermoinps=nothing, kwargs...)
             
  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  lattice    = zeros(3,3)
  mas        = [i.m for i in bdys]
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  if thermoinps != nothing && thermostat != nothing
    NVT = true
  else
    NVT = false
  end

  simu = Simulation(
    bdys, pars, mols, [], [], [],
    save, PBC, NC, lattice, mas, potVars,
    NVT, thermostat, thermoinps
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end

function runMD(EoM, cell::MyCell, tspan::Tuple{Float64, Float64}, dt::Float64; 
               save="full", thermostat=nothing, thermoinps=nothing, kwargs...)

  bdys       = makeBdys(cell)
  potVars    = EoM(cell)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  if thermoinps != nothing && thermostat != nothing
    NVT = true
  else
    NVT = false
  end

  simu = Simulation(
    bdys, pars, mols, [], [], [],
    save, cell.PBC, cell.NC, cell.lattice, cell.masses, potVars,
    NVT, thermostat, thermoinps
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end