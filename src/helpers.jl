
#Calculate the shift frequency for a molecule at vib energy E
freqShiftMorse(v0, D, E) = v0 * ( (D - E) / D )^0.5

#Calculate the energy for a given shifted frequency
engyShiftMorse(v0, D, v) = D - D * (v/v0)^2

function mkvar(x)
  fn = Symbol(x)
  X  = @eval $fn
  return X
end

function translateBdys!(bdys, v)
  for i in bdys
    i.r .+= v
  end
end

function diffDotSqrt(v2, v1)
  rvec = v2 - v1
  r    = sqrt(dot(rvec, rvec))
  (r, rvec)
end

function getAngle(r1, r2)
  x  = (dot(r1, r2)) / (norm(r1) * norm(r2))
  x  = round(x, digits=15)
  
  acos(x)
end

function getAngleCO(co1, co2)
  r1 = co1[1].r - co1[2].r
  r2 = co2[1].r - co2[2].r

  x  = (dot(r1, r2)) / (norm(r1) * norm(r2))
  x  = round(x, digits=15)
  
  acos(x)
end

function getFrame(tj, i::Int64)
  m = tj.m
  s = tj.s
  r = tj.r[i]
  v = tj.v[i]

  [Atom(r[j], v[j], m[j], s[j]) for j in 1:length(m)]
end

function getFrame!(bdys, tj, i::Int64)
  r = tj.r[i]
  v = tj.v[i]

  for j = 1:length(bdys)
    bdys[j].r = r[j]
    bdys[j].v = v[j]
  end

end

function getLastFrame(solu)
  n = length(solu.prob.p.bdys)
  
  new = Atom[]
  for i in 1:n
    r = solu.u[end].x[2][i]
    v = solu.u[end].x[1][i]
    m = solu.prob.p.bdys[i].m
    s = solu.prob.p.bdys[i].s
    push!(new, Atom(r, v, m, s))
  end
  return new
end

function getLastFrame!(bdys, solu)
  N = length(bdys)
  
  for i in 1:N
    bdys[i].r = solu.u[end].x[2][i]
    bdys[i].v = solu.u[end].x[1][i]
    bdys[i].m = solu.prob.p.bdys[i].m
    bdys[i].s = solu.prob.p.bdys[i].s
  end
end

function swapAtoms!(bdys, i, j)
  a = bdys[i].r
  b = bdys[j].r

  bdys[i].r = b
  bdys[j].r = a
end

function centerBdys!(bdys)
  com = CoM(bdys)
  for i in bdys
    i.r -= com
  end
end

function CoM(bdys)
  M = sum([i.m for i in bdys])
  r = sum([i.m*i.r for i in bdys])
  return r ./ M
end

function CoM(pos,mas)
  M = sum(mas)
  r = sum([mas[i]*pos[i] for i in 1:length(mas)])
  return r ./ M
end

function vCoM(bdys)
  M = sum([i.m for i in bdys])
  v = sum([i.m*i.v for i in bdys])
  return v ./ M
end

function vCoM(vel,mas)
  M = sum(mas)
  v = sum([mas[i]*vel[i] for i in 1:length(mas)])
  return v ./ M
end

function zeroVCoM!(bdys)
  N    = length(bdys)
  M    = sum([i.m for i in bdys])
  vcom = vCoM(bdys)
  x    = vcom / N * M
  
  for i in 1:N
    bdys[i].v -= (x / bdys[i].m)
  end
end

function reducedMass(bdys)
  μ = sum([i.m^-1 for i in bdys])^-1
  return μ
end

function reducedMass(mas::Vector{Float64})
  μ = sum([m^-1 for m in mas])^-1
  return μ
end

function swapIso!(bdys, swap, mas)
  for i in 1:length(swap)
    j         = swap[i]
    bdys[j].m = mas[i]
  end
end

function vibExcite!(mol, eignvec, E)
  M = [i.m for i in mol for j in 1:3]
  v = @. sqrt( 2E / M ) * eignvec

  vcom = vCoM(mol)
  for i in 1:3:length(v)
    j::UInt32 = (i+2)/3
    mol[j].v += v[i:i+2]
  end 

  if !isapprox(vCoM(mol), vcom; atol=1e-8)
    println("Uh Oh")
    zeroVCoM!(mol)
  end

end

function transExcite!(mol, KE)
  r  = randVector()
  ke = KE * 8.6173e-5 #convert Kelvin to eV

  for i in mol
    i.v += sqrt(2ke / i.m) .* r
  end
end

function getTransEnergy(mol)
  μ = reducedMass(mol)
  v = vCoM(mol)
  E = 0.5 * μ * dot(v,v)
  return E
end

function getTransEnergy(pos,vel,mas)
  μ = reducedMass(mas)
  v = vCoM(vel,mas)
  E = 0.5 * μ * dot(v,v)
  return E
end

function getRotEnergy(mol)
  vcom = vCoM(mol)
  com  =  CoM(mol)
  E    = 0.0

  for i in mol
    r  = i.r - com
    v  = i.v - vcom
    w  = cross(r,v) / dot(r,r)
    I  = i.m * dot(r,r) 
    E += 0.5 * I * dot(w,w)
  end

  return E
end

function getRotEnergy(pos,vel,mas)
  vcom = vCoM(vel,mas)
  com  =  CoM(pos,mas)
  E    = 0.0

  for i in 1:length(mas)
    r  = pos[i] - com
    v  = vel[i] - vcom
    w  = cross(r,v) / dot(r,r)
    I  = mas[i] * dot(r,r)
    E += 0.5 * I * dot(w,w)
  end
  return E
end

function getVibEnergy(mol, eignvec; pot=nothing)

  E = 0.0
  for i in 1:3:length(eignvec)
    j::UInt32 = (i+2)/3
    ehat      = eignvec[i:i+2] / norm(eignvec[i:i+2])
    E        += 0.5 * mol[j].m * dot(mol[j].v, ehat)^2
  end

  if pot != nothing
    E += getPotEnergy(pot, mol)
  end

  return E
end

#The general method above will be great for more complicated molecules
# but it is very slow. While below only works for CO but is fast.

function getCOVibEnergy(mol; pot=nothing)
  diff = mol[2].r - mol[1].r
  rhat = normalize(diff)
  v1   = dot(mol[1].v, rhat)
  v2   = dot(mol[2].v, rhat)
  E    = 0.5*mol[1].m*v1^2 + 0.5*mol[2].m*v2^2

  if pot != nothing
    E += getPotEnergy(pot, mol)
  end

  return E
end

function getCOVibEnergy(pos,vel,mas; pot=nothing)
  diff = pos[2] - pos[1]
  rhat = normalize(diff)
  v1   = dot(vel[1], rhat)
  v2   = dot(vel[2], rhat)
  E    = 0.5*mas[1]*v1^2 + 0.5*mas[2]*v2^2

  if pot != nothing
    x0 = [j for i in pos for j in i]
    E += pot(true, nothing, x0, optVars([[1,2]], []))
  end

  E
end

function getPotEnergy(EoM, bdys::Vector{MyAtoms})
  x0, vars = prep4pot(EoM, bdys)
  energy   = EoM(true, nothing, x0, vars)
  
  energy
end

function getPotEnergy(EoM, cell::MyCell)
  x0, vars = prep4pot(EoM, cell)
  energy   = EoM(true, nothing, x0, vars)
  
  energy
end

function getForces(EoM, bdys::Vector{MyAtoms})
  x0, vars = prep4pot(EoM, bdys)
  G        = zero(x0)
  EoM(nothing, G, x0, vars)

  [-G[i:i+2] for i = 1:3:length(G)]
end

function getForces(EoM, cell::MyCell)
  x0, vars = prep4pot(EoM, cell)
  G        = zero(x0)
  EoM(nothing, G, x0, vars)

  [-G[i:i+2] for i = 1:3:length(G)]
end

function getSurfaceMolecules(bdys; α=nothing)
  _, mols = getPairs(bdys)

  pts  = [i.r for i in bdys]

  A    = alphashape(pts; α=α)
  
  i = vcat(A.perimeter...) |> unique
  j = findall.(e -> e in i, mols) |> (x -> findall(e -> !isempty(e), x))

  surf = bdys[vcat(mols[j]...)]

  surf
end

function pickRandomMol(bdys, loc)

  surf = getSurfaceMolecules(bdys)
  bulk = [i for i in bdys if !(i in surf)]

  if loc == "surf"

    _,mols = getPairs(surf)
    return surf[rand(mols)]

  elseif loc == "bulk"

    _,mols = getPairs(bulk)
    return bulk[rand(mols)]

  end
  
end 