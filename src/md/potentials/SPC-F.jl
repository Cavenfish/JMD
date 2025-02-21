"""
SPC-F
"""

struct _SPCF_PotVars{F<:Float64} <: PotVars
  Kb::F
  req::F
  Kθ::F
  θeq::F
  σ::F
  ϵ::F
  Qo::F
  Qh::F
end

SPCF(bdys) = _SPCF_PotVars(
  48.05913,
  1.0,
  3.97,
  1.910611,
  3.145, 
  0.007,
  -2.959855,
  0.5 * 2.959855
)

function SPCF(dv, v, u, p, t)

  # initialize things
  E = 0.0
  F = zero(u)
  P = p.potVars

  for mol in p.mols
    o, h1, h2 = mol

    E += _harmonicBond!(F, u, o, h1, P.Kb, P.req)
    E += _harmonicBond!(F, u, o, h2, P.Kb, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.Kθ, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    E += _vdw!(F, u, o1, o2, P.ϵ, P.σ)

    E += _Coulomb!(F, u, o1, o2, P.Qo, P.Qo)

    if any(p.PBC)
      L  = eachrow(p.lattice)
      E +=_pbcInteractions!(F, u, o1, o2, _vdw, L, p.NC, (P.ϵ, P.σ))
      E +=_pbcInteractions!(F, u, o1, o2, _Coulomb, L, p.NC, (P.Qo, P.Qo))
    end

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(F, u, i, j, P.Qh, P.Qh)

        if any(p.PBC)
          L  = eachrow(p.lattice)
          E +=_pbcInteractions!(F, u, i, j, _Coulomb, L, p.NC, (P.Qh, P.Qh))
        end

      end
    end

    for i in [h1,h2]
      E += _Coulomb!(F, u, o2, i, P.Qo, P.Qh)

      if any(p.PBC)
        L  = eachrow(p.lattice)
        E +=_pbcInteractions!(F, u, o2, i, _Coulomb, L, p.NC, (P.Qo, P.Qh))
      end
    end

    for i in [h3,h4]
      E += _Coulomb!(F, u, o1, i, P.Qo, P.Qh)

      if any(p.PBC)
        L  = eachrow(p.lattice)
        E +=_pbcInteractions!(F, u, o1, i, _Coulomb, L, p.NC, (P.Qo, P.Qh))
      end
    end

  end

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function SPCF(F, G, y0, p)

  # initialize things
  E      = 0.0
  u      = Vector[]
  P      = p.potVars
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  for mol in p.mols
    o, h1, h2 = mol

    E += _harmonicBond!(forces, u, o, h1, P.Kb, P.req)
    E += _harmonicBond!(forces, u, o, h2, P.Kb, P.req)
    E += _harmonicBondAngle!(forces, u, h1, o, h2, P.Kθ, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    E += _vdw!(forces, u, o1, o2, P.ϵ, P.σ)

    E += _Coulomb!(forces, u, o1, o2, P.Qo, P.Qo)

    if any(p.PBC)
      L  = eachrow(p.lattice)
      E +=_pbcInteractions!(forces, u, o1, o2, _vdw, L, p.NC, (P.ϵ, P.σ))
      E +=_pbcInteractions!(forces, u, o1, o2, _Coulomb, L, p.NC, (P.Qo, P.Qo))
    end

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(forces, u, i, j, P.Qh, P.Qh)
        if any(p.PBC)
          L  = eachrow(p.lattice)
          E +=_pbcInteractions!(forces, u, i, j, _Coulomb, L, p.NC, (P.Qh, P.Qh))
        end
      end
    end

    for i in [h1,h2]
      E += _Coulomb!(forces, u, o2, i, P.Qo, P.Qh)
      if any(p.PBC)
        L  = eachrow(p.lattice)
        E +=_pbcInteractions!(forces, u, o2, i, _Coulomb, L, p.NC, (P.Qo, P.Qh))
      end
    end

    for i in [h3,h4]
      E += _Coulomb!(forces, u, o1, i, P.Qo, P.Qh)
      if any(p.PBC)
        L  = eachrow(p.lattice)
        E +=_pbcInteractions!(forces, u, o1, i, _Coulomb, L, p.NC, (P.Qo, P.Qh))
      end
    end

  end

  if G != nothing
    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

  if F != nothing
    return E
  end

end