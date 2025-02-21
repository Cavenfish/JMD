const libmbx       = dlopen(joinpath(@__DIR__, "libmbx.so"), Libdl.RTLD_GLOBAL) 
const E_MBX2JMD    = 0.043 # kcal/mol --> eV
const MBX_GAS_JSON = joinpath(@__DIR__, "gas.json")
const MBX_PBC_JSON = joinpath(@__DIR__, "pbc.json")

function mbx_get_monomer_info(at_nams)
  monomers = Dict()
  monomers["nh4+"] = ["N","H","H","H","H"]
  monomers["nh3"] = ["N","H","H","H"]
  monomers["ch4"] = ["C","H","H","H","H"]
  monomers["pf6-"] = ["P","F","F","F","F","F","F"]
  monomers["co2"] = ["C","O","O"]
  monomers["li"] = ["Li"]
  monomers["na"] = ["Na"]
  monomers["k"] = ["K"]
  monomers["rb"] = ["Rb"]
  monomers["cs"] = ["Cs"]
  monomers["f"] = ["F"]
  monomers["cl"] = ["Cl"]
  monomers["br"] = ["Br"]
  monomers["i"] = ["I"]
  monomers["ar"] = ["Ar"]
  monomers["he"] = ["He"]
  monomers["h2"] = ["H","H"]
  monomers["h2o"] = ["O","H","H"]
  monomers["n2o5"] = ["O","N","N","O","O","O","O"]

  mon_nams = String[]
  num_ats  = Int32[]

  tmp = values(monomers) |> (x -> length.(x)) |> unique |> sort |> reverse

  a = 1
  b = 0
  while a <= length(at_nams)
    for i in tmp
      a+i-1 > length(at_nams) && continue
      cur_ats = at_nams[a:a+i-1]
      for j in keys(monomers)
        if monomers[j] == cur_ats
          push!(mon_nams, j)
          push!(num_ats, length(monomers[j]))
          a += length(monomers[j])
          break
        end
      end
    end
    b += 1
    if b > length(at_nams)
      @error "Check xyz file"
    end
  end

  mon_nams, num_ats
end

function mbx_initialize_system(vars)
  sym = dlsym(libmbx, :initialize_system_py_)

  @ccall $sym(
    vars.xyz::Ptr{Cdouble},
    vars.num_ats::Ptr{Cint},
    vars.at_nams::Ptr{Ptr{Cchar}},
    vars.mon_nams::Ptr{Ptr{Cchar}},
    vars.num_mon::Ref{Cint},
    vars.json::Ptr{Cchar}
  )::Cvoid
end

function mbx_initialize_system(xyz, num_ats, at_nams, mon_nams, num_mon, json)
  sym = dlsym(libmbx, :initialize_system_py_)

  @ccall $sym(
    xyz::Ptr{Cdouble},
    num_ats::Ptr{Cint},
    at_nams::Ptr{Ptr{Cchar}},
    mon_nams::Ptr{Ptr{Cchar}},
    num_mon::Ref{Cint},
    json::Ptr{Cchar}
  )::Cvoid
end

function mbx_get_energy(xyz, nats)
  sym = dlsym(libmbx, :get_energy_)

  E   = Ref{Cdouble}(0)

  @ccall $sym(
    xyz::Ptr{Cdouble},
    nats::Ref{Cint},
    E::Ref{Cdouble}
  )::Cvoid

  E[] * E_MBX2JMD
end

function mbx_get_energy_grad!(G, xyz, nats)
  sym = dlsym(libmbx, :get_energy_g_)

  E   = Ref{Cdouble}(0)

  @ccall $sym(
    xyz::Ptr{Cdouble}, 
    nats::Ref{Cint}, 
    E::Ref{Cdouble},
    G::Ptr{Cdouble}
  )::Cvoid

  G .*= E_MBX2JMD

  E[] * E_MBX2JMD
end

function mbx_get_total_dipole(μ)
  sym = dlsym(libmbx, :get_total_dipole_)

  @ccall $sym(μ::Ptr{Cdouble})::Cvoid
end

function mbx_finalize_system()
  sym = dlsym(libmbx, :finalize_system_)

  @ccall $sym()::Cvoid
end

function mbx_set_box(box)
  sym = dlsym(libmbx, :set_box_)

  l = length(box)

  @ccall $sym(
    l::Ref{Cint},
    box::Ptr{Cdouble}
  )::Cvoid

end