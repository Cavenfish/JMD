
struct NVEsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  save::String
  m::Vector
  potVars::PotVars
end

struct NVTsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  temp::Vector
  save::String
  m::Vector
  potVars::PotVars
  thermostat!::Function
  thermoInps
end

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

function run(EoM, bdys::Vector{MyAtoms}, NVT::Bool, tspan::Tuple{Float64, Float64},
             dt::Float64; save="full", kwargs...)
             
  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  lattice    = zeros(3,3)
  mas        = [i.m for i in bdys]
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  if NVT
    tFunc = kwargs[:thermostat]
    tInps = kwargs[:thermoinps]
  else
    tFunc = nothing
    tInps = nothing
  end

  simu = Simulation(
    bdys, pars, mols, [], [], [],
    save, PBC, NC, lattice, mas, potVars,
    NVT, tFunc, tInps
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end

function run(EoM, cell::MyCell, NVT::Bool, tspan::Tuple{Float64, Float64},
             dt::Float64; save="full", kwargs...)

  bdys       = makeBdys(cell)
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  if NVT
    tFunc = kwargs[:thermostat]
    tInps = kwargs[:thermoinps]
  else
    tFunc = nothing
    tInps = nothing
  end

  simu = Simulation(
    bdys, pars, mols, [], [], [],
    save, cell.PBC, cell.NC, cell.lattice, cell.masses, potVars,
    NVT, tFunc, tInps
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end

function runNVE(EoM, tspan, dt, bdys; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]
  
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], save, mas, potVars)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end

function runNVE(EoM, t1, t2, dt, bdys; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]
  tspan = (t1,t2)

  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], save, mas, potVars)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false, save_start=false)

  return solu
end

function runNVT(EoM, tspan, dt, bdys, thermostat, thermoInps; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]

  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  simu       = NVTsimu(bdys, pars, mols, [], [], [], save, mas, potVars, thermostat, thermoInps)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end