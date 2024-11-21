"""
SCME/f potential


SCME Papers:
https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp52097h

https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00598
"""

struct _SCME_PotVars <: PotVars
  molnum
  cm
  dp
  qp
  op
  hp
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
  irigidmolecules
  convcrit
  dqdms
  d1
  d2
  dd1
  dd2
  d1vQM
  d2vQM
  d3vQM
  d4vQM
  d5vQM
  te
  td
  dpQM
  qpQM
  u_ES
  u_DS
  u_RP
  fa_ES
  fa_DS
  fa_RP
  Ar
  Br
  Cr
  rc_Core
  atoms_pbc
  tags
end

function SCME(bdys)

  isdefined(JMD, :userTOML) ? tomlFile = userTOML : tomlFile = SCMEtoml
  tom       = TOML.parsefile(tomlFile)
  te        = tom["te"]
  td        = tom["td"]
  Ar        = tom["Ar"]
  Br        = tom["Br"]
  Cr        = tom["Cr"]
  cell      = tom["cell"]
  nsys      = tom["nsys"]
  system    = tom["system"]
  useDMS    = tom["useDMS"]
  rc_Elec   = tom["rc_Elec"]
  rc_Disp   = tom["rc_Disp"]
  rc_Core   = tom["rc_Core"]
  atoms_pbc = tom["atoms_pbc"] 


  n_atoms = length(bdys)
  coords  = [j for i in bdys for j in i.r]
  molnum  = floor(Int64, n_atoms/3)
  cm      = zeros(molnum, 3)

  # Moments - irreducible
  dp   = zeros(molnum, 3)
  dp0  = zeros(molnum, 3)
  dp01 = zeros(molnum, 3)
  qp   = zeros(molnum, 3,3)
  qp0  = zeros(molnum, 3,3)
  qp01 = zeros(molnum, 3,3)
  op   = zeros(molnum, 3,3,3)
  hp   = zeros(molnum, 3,3,3,3)

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

  # Maybe move these things into TOML file
  iSlab           = false
  irigidmolecules = false
  convcrit        = 1e-11

  tags = repeat([false], n_atoms)

  # Load SCME init function as symbol
  sym = dlsym(libscme, :scme_initilization)
  

  """
  int n_atoms,
  double coords[],
  double lattice[],
  int system[],
  int Nsys,
  double cm_init[][3],
  // assigned multipoles
  double dp_init[][3],
  double qp_init[][3][3],
  double op_init[][3][3][3],
  double hp_init[][3][3][3][3],
  // assigned polarizabilities
  double dd_init[][3][3],
  double dq_init[][3][3][3],
  double qq_init[][3][3][3][3],
  // assigned hyperpolarizabilities
  double hpol_init[][3][3][3],
  // calculated field and its first derivative
  double d1vH_init[][3],
  double d2vH_init[][3][3],
  // number of cells in x,y,z directions
  int NC[],
  // EOJ: ?
  double te,
  double rc_Elec,
  // // useDMS variables
  double d1_init[][3],
  double d2_init[][3],
  double dd1_init[][3][3][3],
  double dd2_init[][3][3][3],
  double qatoms[],
  double dqdms[][3][3][3],
  bool useDMS,
  bool pbc[3],
  bool tags[]   // tagging molecules as QM or MM
  """

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

  vars = _SCME_PotVars(
    molnum,
    cm,
    dp,
    qp,
    op,
    hp,
    dp0,
    qp0,
    dp01,
    qp01,
    NC,
    cell,
    d1vH,
    d2vH,
    dd,
    dq,
    hpol,
    qq,
    rc_Elec,
    iSlab,
    irigidmolecules,
    convcrit,
    dqdms,
    d1,
    d2,
    dd1,
    dd2,
    d1vQM,
    d2vQM,
    d3vQM,
    d4vQM,
    d5vQM,
    te,
    td,
    dpQM,
    qpQM,
    u_ES,
    u_DS,
    u_RP,
    fa_ES,
    fa_DS,
    fa_RP,
    Ar,
    Br,
    Cr,
    rc_Core,
    atoms_pbc,
    tags
  )

  vars
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