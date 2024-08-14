
struct Traj <: MyTraj
  t::Vector{Float64}
  E::Vector{Float64}
  T::Vector{Float64}
  F::Vector
  r::Vector
  v::Vector
  m::Vector{Float64}
  s::Vector{Char}
end

function processDynamics(solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  s = [i.s for i in solu.prob.p.bdys]
  T = [getTemp(m, i, kB, length(m)) for i in v]

  return Traj(t, E, T, F, r, v, m, s)
end

function processDynamics!(tj, solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  T = [getTemp(m, i, kB, length(m)) for i in v]

  push!(tj.t, t...)
  push!(tj.E, E...)
  push!(tj.T, T...)
  push!(tj.F, F...)
  push!(tj.r, r...)
  push!(tj.v, v...)
end

function processTmpFiles(files; kwargs...)
  
  f    = open(files[1], "r")
  solu = deserialize(f)
  tj   = processDynamics(solu; kwargs...) 
  close(f)

  # Free memory
  @free solu

  for file in files[2:end]
    f    = open(file, "r")
    solu = deserialize(f)
    processDynamics!(tj, solu; kwargs...)
    close(f)

    #Free memory
    @free solu
  end

  return tj
end

function trackVACF(files, safe)

  df = DataFrame()

  for file in files
    f    = open(file, "r")
    solu = deserialize(f)
    close(f)

    v,m  = getVelMas(solu)

    inp  = vacfInps(v[1:safe], m, 1e15, false, HannM, 8, true)
    out  = VDOS(inp)
    col  = replace(file, ".tmp" => "")

    if isempty(df)
      df[!, "v"] = out.v
    end

    df[!, col] = out.I

    @free solu
  end

  return df
end

function trackEnergyDissipation(traj, pot, mol; eignvec=nothing)

  m = traj.m
  N = length(traj.t)
  n = length(m)

  molVib, molRot, molTra = Float64[], Float64[], Float64[]
  avgVib, avgRot, avgTra = Float64[], Float64[], Float64[]

  others = [[i,i+1] for i in 1:2:n if !(i in mol)]

  zpe = zeros(N)

  for i in 1:N
    r = traj.r[i]
    v = traj.v[i]

    mvib = getVibEnergy(r[mol], v[mol], m[mol]; pot=pot)
    mrot = getRotEnergy(  r[mol], v[mol], m[mol])
    mtra = getTransEnergy(r[mol], v[mol], m[mol])

    push!(molVib, mvib)
    push!(molRot, mrot)
    push!(molTra, mtra)

    avib,arot,atra = 0.0,0.0,0.0
    for j in others
      avib += getVibEnergy(r[j], v[j], m[j]; pot=pot)
      arot += getRotEnergy(  r[j], v[j], m[j])
      atra += getTransEnergy(r[j], v[j], m[j])
    end

    if eignvec != nothing
      bdys = getFrame(traj, i)
      
      for e in eignvec
        zpe[i] += getVibEnergy(bdys, e)
      end
    end

    push!(avgVib, avib/n)
    push!(avgRot, arot/n)
    push!(avgTra, atra/n)
  end

  df = DataFrame(  time = traj.t,
                   temp = traj.T,
                 molVib = molVib,
                 molRot = molRot,
                 molTra = molTra,
                 avgVib = avgVib,
                 avgRot = avgRot,
                 avgTra = avgTra,
                   vzpe = zpe)

  return df
end

function trackAllVibEnergy(tj; pot=nothing)

  N   = length(tj.r[1])

  pos = tj.r[102]
  vel = tj.v[102]
  mas = tj.m

  tmp = [getCOVibEnergy(pos[i:i+1], vel[i:i+1], mas[i:i+1]) for i in 1:2:N]

  # This is the excited molecule
  ind = findall(e -> e == maximum(tmp), tmp)[1]

  dis = DataFrame()
  vib = DataFrame()

  dis[!, "time"] = tj.t
  vib[!, "time"] = tj.t

  for i in 1:div(N,2)
    str = "mol$i"
    
    dis[!, str] = zero(tj.t)
    vib[!, str] = zero(tj.t)
  end

  for i in 1:length(tj.t)
    for j in 1:div(N,2)
      #properties of excited molecule
      epos = tj.r[i][ind*2-1:ind*2]
      emas = tj.m[ind*2-1:ind*2]
      ecom = CoM(epos, emas)

      a   = j*2-1
      b   = j*2

      #properties of other molecule
      pos = tj.r[i][a:b]
      vel = tj.v[i][a:b]
      mas = tj.m[a:b]
      com = CoM(pos, mas)
      
      #distance and vib energy
      d = norm(com - ecom)
      v = getCOVibEnergy(pos, vel, mas; pot=pot)

      #store in DFs
      dis[i,j+1] = d
      vib[i,j+1] = v

    end
  end

  return dis, vib
end

function trackRadialEnergy(tj; pot=nothing, Rs=[[3,5],[5,9],[9,30]])

  #I'm working in let-end blocks to local scope some
  # variables so they don't get double referenced
  
  #Find indicies of excited molecule
  ind = let
    v = tj.v[102]

    i = findfirst(e -> e == maximum(v), v)

    isodd(i) ? [i, i+1] : [i-1, i]
  end 

  #Get all molecules (only works for CO)
  mols = let
    r  = tj.r[1]
    d  = [[norm(j - i) for j in r] for i in r]
    
    findall.(e -> e < 2, d) |> unique
  end

  #Remove excited molecule from mols list
  findfirst(e -> e == ind, mols) |> (x -> popat!(mols, x))

  #Pre-allocate distance array
  d = zeros(length(mols))

  #Pre-define mass since it is constant in time
  m = tj.m

  #Pre-build DataFrame
  df = DataFrame()

  df[!, "time"] = tj.t
  for i in 1:length(Rs)
    df[!, "traR$i"] = zero(tj.t)
    df[!, "rotR$i"] = zero(tj.t)
    df[!, "vibR$i"] = zero(tj.t)
  end

  #Very complicated nested loop
  for t in 1:length(tj.t)
    d .*= 0.0

    r = tj.r[t]
    v = tj.v[t]

    com  =  CoM(r[ind], m[ind])
    coms = [CoM(r[mol], m[mol]) for mol in mols]
    d  .+= [norm(i - com) for i in coms]

    for k in 1:length(Rs)
      low, hi = Rs[k]

      #Energy vars
      tra = 0.0
      rot = 0.0
      vib = 0.0

      #get indicies of mols in radial slice
      l = findall(e -> low < e < hi, d)

      for i in 1:length(l)
        j = mols[l[i]]

        tra += getTransEnergy(r[j], v[j], m[j])
        rot += getRotEnergy(  r[j], v[j], m[j])
        vib += getCOVibEnergy(r[j], v[j], m[j]; pot=pot)

      end#mols within low-hi loop

      df[t, "traR$k"] += tra / length(l)
      df[t, "rotR$k"] += rot / length(l)
      df[t, "vibR$k"] += vib / length(l)

    end#radial slices loop
  end#time loop

  df
end#function