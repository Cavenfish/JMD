const SCMEtoml   = joinpath(@__DIR__, "scme.toml")
const libscme    = dlopen(joinpath(@__DIR__, "libscme.so"), Libdl.RTLD_GLOBAL)


function scme_initialize_system()
  sym = dlsym(libscme, :scme_initilization)

  z3 = zeros(3)

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
    z3::Ptr{Cint},
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

function scme_scf_step()
  sym = dlsym(libscme, :scme_scf_step)

  @ccall $sym(
    molnum::Cint,
    cm::Ptr{Cdouble},
    dp0::Ptr{Cdouble},
    qp0::Ptr{Cdouble},
    dp01::Ptr{Cdouble},
    qp01::Ptr{Cdouble},
    NC::Ptr{Cint},
    cell::Ptr{Cdouble},
    d1vH::Ptr{Cdouble},
    d2vH::Ptr{Cdouble},
    dd::Ptr{Cdouble},
    dq::Ptr{Cdouble},
    hpol::Ptr{Cdouble},
    qq::Ptr{Cdouble},
    rc_Elec::Cdouble,
    iSlab::Cint,
    convcrit::Cdouble,
    d1vQM::Ptr{Cdouble},
    d2vQM::Ptr{Cdouble},
    te::Cdouble,
    dpQM::Ptr{Cdouble},
    qpQM::Ptr{Cdouble},
    atoms_pbc::Ptr{Cint},
    tags::Ptr{Cint}
  )::Cvoid

end


function scme_ES()
  sym = dlsym(libscme, :scme_ES)

  dpoleQM = dp - 0.5 * dpQM
  qpoleQM = qp - 0.5 * qpQM

  @ccall $sym(
    molnum::Cint,
    coords::Ptr{Cdouble},
    cm::Ptr{Cdouble},
    rc_Elec::Cdouble,
    system::Ptr{Cint},
    iSlab::Cint,
    irigidmolecules::Cint,
    addCore::Cint,
    NC::Ptr{Cint},
    cell::Ptr{Cdouble},
    dp::Ptr{Cdouble},
    qp::Ptr{Cdouble},
    op::Ptr{Cdouble},
    hp::Ptr{Cdouble},
    dp0::Ptr{Cdouble},
    qp0::Ptr{Cdouble},
    u_ES::Cdouble,
    fa_ES::Ptr{Cdouble},
    te::Cdouble,
    useDMS::Cint,
    qatoms::Ptr{Cdouble},
    dqdms::Ptr{Cdouble},
    d1::Ptr{Cdouble},
    d2::Ptr{Cdouble},
    dd1::Ptr{Cdouble},
    dd2::Ptr{Cdouble},
    d1vQM::Ptr{Cdouble},
    d2vQM::Ptr{Cdouble},
    d3vQM::Ptr{Cdouble},
    d4vQM::Ptr{Cdouble},
    d5vQM::Ptr{Cdouble},
    dpoleQM::Ptr{Cdouble},
    qpoleQM::Ptr{Cdouble},
    atoms_pbc::Ptr{Cint},
    tags::Ptr{Cint}
  )::Cvoid

end

function scme_Disp()
  sym = dlsym(libscme, :scme_Disp)

  @ccall $sym(
    molnum::Cint,
    coords::Ptr{Cdouble},
    cell::Ptr{Cdouble},
    atoms_pbc::Ptr{Cint},
    u_DS::Cdouble,
    fa_DS::Ptr{Cdouble},
    td::Ptr{Cdouble},
    tags::Ptr{Cint},
    NC::Ptr{Cint},
    rc_Disp::Cdouble
  )::Cvoid

end

function scme_Core()
  sym = dlsym(libscme, :scme_Core)

  @ccall $sym(
    molnum::Cint,
    coords::Ptr{Cdouble},
    tags::Ptr{Cint},
    cell::Ptr{Cdouble},
    atoms_pbc::Ptr{Cint},
    addCore::Cint,
    u_RP::Cdouble,
    fa_RP::Ptr{Cdouble},
    Ar::Ptr{Cdouble},
    Br::Ptr{Cdouble},
    Cr::Ptr{Cdouble},
    rc_Core::Cdouble
  )::Cvoid

end