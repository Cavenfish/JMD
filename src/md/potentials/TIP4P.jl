"""
TIP4P/2005f 
"""

struct _TIP4P_PotVars{F<:Float64} <: PotVars
  D::F
  a::F
  req::F
  K::F
  θeq::F
  ϵoo::F
  σoo::F
  drel::F
  Qh::F
  Qm::F
end

#PotVar building function 
TIP4P(bdys) = _TIP4P_PotVars(
  4.48339,    # eV
  2.287,      # \AA
  0.9419,     # \AA
  3.81209321, # eV rad^-2  1.16123e-3 * (2pi/360)^2
  1.87448,    # rad        107.4  * (2pi/360)
  8.03e-3,    # eV
  3.1644,     # \AA
  0.13194,    # \AA
  2.1113635,  # 
  -2 * 2.1113635
)

function TIP4P(dv, v, u, p, t)

  # initialize things
  E = 0.0
  F = zero(u)
  P = p.potVars


  for mol in p.mols
    o, h1, h2 = mol

    E += _Morse!(F, u, o, h1, P.D, P.a, P.req)
    E += _Morse!(F, u, o, h2, P.D, P.a, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.K, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    E += _vdw!(F, u, o1, o2, P.ϵoo, P.σoo)

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(F, u, i, j, P.Qh, P.Qh)
      end
    end

    E += _getMforces!(F, u, par[1], par[2], P.drel, P.Qh, P.Qm)

  end

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function TIP4P(F, G, y0, p)
  
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

    E += _Morse!(forces, u, o, h1, P.D, P.a, P.req)
    E += _Morse!(forces, u, o, h2, P.D, P.a, P.req)
    E += _harmonicBondAngle!(forces, u, h1, o, h2, P.K, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    E += _vdw!(forces, u, o1, o2, P.ϵoo, P.σoo)

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(forces, u, i, j, P.Qh, P.Qh)
      end
    end

    E += _getMforces!(forces, u, par[1], par[2], P.drel, P.Qh, P.Qm)

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

function _getMforces!(F, u, w1, w2, drel, Qh, Qm)
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  # Get Angle
  θ1  = getAngle(u[h1] - u[o1], u[h2] - u[o1])
  θ2  = getAngle(u[h3] - u[o2], u[h4] - u[o2])

  # Get r vectors
  r1o  = u[h1] - u[o1]
  r2o  = u[h2] - u[o1]
  r3o  = u[h3] - u[o2]
  r4o  = u[h4] - u[o2]

  # Get Norms
  d1o  = norm(r1o)
  d2o  = norm(r2o)
  d3o  = norm(r3o)
  d4o  = norm(r4o)

  # Get dm
  dm1 = drel * (d1o*cos(θ1/2) + d2o*cos(θ1/2))
  dm2 = drel * (d3o*cos(θ2/2) + d4o*cos(θ2/2))

  # Get bisectors
  rbi1 = r1o/d1o + r2o/d2o
  rbi2 = r3o/d3o + r4o/d4o

  # Get weights
  wh1 = dm1 / d1o / norm(rbi1)
  wh2 = dm1 / d2o / norm(rbi1)
  wh3 = dm2 / d3o / norm(rbi2)
  wh4 = dm2 / d4o / norm(rbi2)

  # Get m site vectors
  m1 = u[o1] + wh1*r1o + wh2*r2o
  m2 = u[o2] + wh3*r3o + wh4*r4o

  # H1 -- M2
  E,f    = _Coulomb(u[h1], m2, Qh, Qm)
  F[h1] -= f
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  # H2 -- M2
  e,f    = _Coulomb(u[h2], m2, Qh, Qm)
  E     += e
  F[h2] -= f
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  # H3 -- M1
  e,f    = _Coulomb(u[h3], m1, Qh, Qm)
  E     += e
  F[h3] -= f
  F[h1] += f * wh1
  F[h2] += f * wh2
  F[o1] += f * (1 - wh1 - wh2)

  # H4 -- M1
  e,f    = _Coulomb(u[h4], m1, Qh, Qm)
  E     += e
  F[h4] -= f
  F[h1] += f * wh1
  F[h2] += f * wh2
  F[o1] += f * (1 - wh1 - wh2)

  # M1 -- M2
  e,f    = _Coulomb(m1, m2, Qm, Qm)
  E     += e
  F[h1] -= f * wh1
  F[h2] -= f * wh2
  F[o1] -= f * (1 - wh1 - wh2)
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  E
end