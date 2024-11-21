const SCMEtoml   = joinpath(@__DIR__, "scme.toml")
const libscme    = dlopen(joinpath(@__DIR__, "libscme.so"), Libdl.RTLD_GLOBAL)


function scme_initialize_system!()
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

function scme_scf_step!(vars)
  sym = dlsym(libscme, :scf_step)

  @ccall $sym(
    vars.molnum::Cint,
    vars.cm::Ptr{Cdouble},
    vars.dp0::Ptr{Cdouble},
    vars.qp0::Ptr{Cdouble},
    vars.dp01::Ptr{Cdouble},
    vars.qp01::Ptr{Cdouble},
    vars.NC::Ptr{Cint},
    vars.cell::Ptr{Cdouble},
    vars.d1vH::Ptr{Cdouble},
    vars.d2vH::Ptr{Cdouble},
    vars.dd::Ptr{Cdouble},
    vars.dq::Ptr{Cdouble},
    vars.hpol::Ptr{Cdouble},
    vars.qq::Ptr{Cdouble},
    vars.rc_Elec::Cdouble,
    vars.iSlab::Cint,
    vars.convcrit::Cdouble,
    vars.d1vQM::Ptr{Cdouble},
    vars.d2vQM::Ptr{Cdouble},
    vars.te::Cdouble,
    vars.dpQM::Ptr{Cdouble},
    vars.qpQM::Ptr{Cdouble},
    vars.atoms_pbc::Ptr{Cint},
    vars.tags::Ptr{Cint}
  )::Cvoid

end


function scme_ES!(vars, coords)
  sym = dlsym(libscme, :scme_ES)

  dpoleQM = vars.dp - 0.5 * vars.dpQM
  qpoleQM = vars.qp - 0.5 * vars.qpQM
  E       = Ref{Cdouble}(0)

  @ccall $sym(
    vars.molnum::Cint,
    coords::Ptr{Cdouble},
    vars.cm::Ptr{Cdouble},
    vars.rc_Elec::Cdouble,
    vars.system::Ptr{Cint},
    vars.iSlab::Cint,
    vars.irigidmolecules::Cint,
    vars.addCore::Cint,
    vars.NC::Ptr{Cint},
    vars.cell::Ptr{Cdouble},
    vars.dp::Ptr{Cdouble},
    vars.qp::Ptr{Cdouble},
    vars.op::Ptr{Cdouble},
    vars.hp::Ptr{Cdouble},
    vars.dp0::Ptr{Cdouble},
    vars.qp0::Ptr{Cdouble},
    E::Ref{Cdouble},
    vars.fa_ES::Ptr{Cdouble},
    vars.te::Cdouble,
    vars.useDMS::Cint,
    vars.qatoms::Ptr{Cdouble},
    vars.dqdms::Ptr{Cdouble},
    vars.d1::Ptr{Cdouble},
    vars.d2::Ptr{Cdouble},
    vars.dd1::Ptr{Cdouble},
    vars.dd2::Ptr{Cdouble},
    vars.d1vQM::Ptr{Cdouble},
    vars.d2vQM::Ptr{Cdouble},
    vars.d3vQM::Ptr{Cdouble},
    vars.d4vQM::Ptr{Cdouble},
    vars.d5vQM::Ptr{Cdouble},
    dpoleQM::Ptr{Cdouble},
    qpoleQM::Ptr{Cdouble},
    vars.atoms_pbc::Ptr{Cint},
    vars.tags::Ptr{Cint}
  )::Cvoid

  E[]
end

function scme_Disp!(vars, coords)
  sym = dlsym(libscme, :scme_Disp)

  E = Ref{Cdouble}(0)

  @ccall $sym(
    vars.molnum::Cint,
    coords::Ptr{Cdouble},
    vars.cell::Ptr{Cdouble},
    vars.atoms_pbc::Ptr{Cint},
    E::Ref{Cdouble},
    vars.fa_DS::Ptr{Cdouble},
    vars.td::Ptr{Cdouble},
    vars.tags::Ptr{Cint},
    vars.NC::Ptr{Cint},
    vars.rc_Disp::Cdouble
  )::Cvoid

  E[]
end

function scme_Core!(vars, coords)
  sym = dlsym(libscme, :scme_Core)

  E = Ref{Cdouble}(0)

  @ccall $sym(
    vars.molnum::Cint,
    coords::Ptr{Cdouble},
    vars.tags::Ptr{Cint},
    vars.atoms_pbc::Ptr{Cint},
    vars.addCore::Cint,
    vars.cell::Ptr{Cdouble},
    E::Ref{Cdouble},
    vars.fa_RP::Ptr{Cdouble},
    vars.Ar::Ptr{Cdouble},
    vars.Br::Ptr{Cdouble},
    vars.Cr::Ptr{Cdouble},
    vars.rc_Core::Cdouble
  )::Cvoid

  E[]
end