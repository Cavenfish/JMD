"""
SCME/f potential


SCME Papers:
https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp52097h

https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00598
"""

struct _SCME_PotVars <: PotVars
  molnum
  cm
  dp0
  qp0
  dp01
  qp01
  NC
  cell
  d1vH
  d2vH
  dd
  dq
  hpol
  qq
  rc_Elec
  iSlab
  convcrit
  d1vQM
  d2vQM
  te
  dpQM
  qpQM
  atoms_pbc
  tags

end

function SCME(bdys)

  n_atoms = length(bdys)
  coords  = [i.r for i in bdys]
  molnum  = floor(Int64, n_atoms/3)

  # Moments - irreducible
  dp = zeros(molnum, 3)
  qp = zeros(molnum, 3,3)
  op = zeros(molnum, 3,3,3)
  hp = zeros(molnum, 3,3,3,3)

  # Polarizabilities - irreducible
  dd   = zeros(molnum, 3,3)
  dq   = zeros(molnum, 3,3,3)
  qq   = zeros(molnum, 3,3,3,3)
  hpol = zeros(molnum, 3,3,3)

  # Higher order fields - irreducible
  d1vH = zeros(molnum, 3)
  d2vH = zeros(molnum, 3,3)

  NC   = [0,0,0]

  d1     = zeros(molnum, 3)
  d2     = zeros(molnum, 3)
  dd1    = zeros(molnum, 3,3,3)
  dd2    = zeros(molnum, 3,3,3)
  qatoms = zeros(n_atoms)
  dqdms  = zeros(molnum, 3,3,3)

  # MM fields and gradients
  d1v = zeros(molnum,3)
  d2v = zeros(molnum,3,3)
  d3v = zeros(molnum,3,3,3)
  d4v = zeros(molnum,3,3,3,3)
  d5v = zeros(molnum,3,3,3,3,3)

  # QM fields and gradients - filled in by GPAW - irreducible
  d1vQM = zeros(molnum,3)
  d2vQM = zeros(molnum,3,3)
  d3vQM = zeros(molnum,3,3,3)
  d4vQM = zeros(molnum,3,3,3,3)
  d5vQM = zeros(molnum,3,3,3,3,3)

  # QM induced dpole and qpole - irreducible
  dpQM  = zeros(molnum,3)
  qpQM  = zeros(molnum,3,3)

  # Energy and force terms
  # EOJ :: This needs a cleanup. We are opening up too many unneccessary arrays(?)
  u_ES = 0.0
  u_DS = 0.0
  u_RP = 0.0

  energy_ES = 0.0
  energy_DS = 0.0
  energy_RP = 0.0

  fa_ES = zeros(molnum*9)
  fa_DS = zeros(molnum*9)
  fa_RP = zeros(molnum*9)

  forces_ES = zeros(molnum*9)
  forces_DS = zeros(molnum*9)
  forces_RP = zeros(molnum*9)

  tags = repeat([false], n_atoms)

  # Load SCME init function as symbol
  sym = dlsym(libscme, :scme_initilization)
  
  # Call SCME init function 
  @ccall $sym(
    n_atoms::Cint,
    coords::Ptr{Cdouble},
    cell::Ptr{Cdouble},
    system::Ptr{Cint},
    nsys::Cint,
    cm::Ptr{Cdouble},
    dp::Ptr{Cdouble}, qp::Ptr{Cdouble}, op::Ptr{Cdouble}, hp::Ptr{Cdouble},
    dd::Ptr{Cdouble}, dq::Ptr{Cdouble}, qq::Ptr{Cdouble},
    hpol::Ptr{Cdouble},
    d1vH::Ptr{Cdouble}, d2vH::Ptr{Cdouble},
    zeros(3)::Ptr{Cint},
    te::Cdouble,
    rc_Elec::Cdouble,
    d1::Ptr{Cdouble},
    d2::Ptr{Cdouble},
    dd1::Ptr{Cdouble},
    dd2::Ptr{Cdouble},
    qatoms::Ptr{Cdouble},
    dqdms::Ptr{Cdouble},
    useDMS::Cint,
    atoms_pbc::Ptr{Cint},
    tags::Ptr{Cint}
  )::Cvoid

end

function SCME(dv, v, u, p, t)


  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function SCME(F, G, y0, p)

  if G != nothing
  
  end

  if F != nothing

  end

end