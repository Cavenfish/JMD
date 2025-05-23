"""
Anneals a cluster, based on given toml input card,
to get crystal structure.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
jldfile = "/home/brian/COjl/10K-MvH.jld2"
cluster = "250co"
cycles = 10

[NVT]
thermostat = "BDP!"
thermoInps = "BDP"
temp = 25
kB = "kB"
tau.num = 100
tau.uni = "fs"
time.num = 5
time.uni = "ps"

[OPT]
algo  = "LBFGS"
kwargs.x_tol = 1e-4
kwargs.f_tol = 1e-12
kwargs.g_tol = 1e-4
kwargs.iterations = 100000

[Saving]
xyz = "crystal.xyz"
"""

function anneal(inpFile::String)
  inp = TOML.parsefile(inpFile)

  #Split dict for easier usage
  cnfg = inp["Settings"]
  nvtd = inp["NVT"]

  #Load up EoM and opt algo
  EoM  = mkvar(cnfg["EoM"])
  algo = mkvar(inp["OPT"]["algo"])()

  #Prep nvt inputs
  T      = nvtd["temp"]
  tau    = nvtd["tau"]["num"] * mkvar(nvtd["tau"]["uni"])
  nvtInp = mkvar(nvtd["thermoInps"])
  inps   = nvtInp(T, mkvar(nvtd["kB"]), tau)
  thermo = mkvar(nvtd["thermostat"])

  #Build opt kwargs dict
  kwargs = Dict()
  for key in keys(inp["OPT"]["kwargs"])
    value = inp["OPT"]["kwargs"][key]
    kwargs[Symbol(key)] = value
  end

  #Get leftover vars
  jld = cnfg["jldfile"]
  clu = cnfg["cluster"]
  N   = cnfg["cycles"]
  t   = nvtd["time"]["num"] * mkvar(nvtd["time"]["uni"])

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = jd[clu]

  #Start annealing
  for i in 1:N

    #Run NVT
    nvt = runNVT(EoM, (0, t), fs, bdys, thermo, inps)

    #Update bdys
    getLastFrame!(bdys, nvt)

    #Free memory
    @free nvt 

    #Optimize geometry
    bdys = opt(EoM, algo, bdys; kwargs...)
    
  end

  writeXyz(inp["Saving"]["xyz"], bdys)

end