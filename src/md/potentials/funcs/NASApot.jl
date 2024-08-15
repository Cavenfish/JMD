"""
NASA water monomer potential

H. Partridge and D. W. Schwenke, J. Chem. Phys. 106, 4618 (1997)
"""

function _Va!(F, u, i, j, D, a, r0)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = D * (exp(-2a*(r-r0)) - 2*exp(-a*(r-r0)))
  f     = D * 2a * (exp(-2a*(r-r0)) - exp(-a*(r-r0))) * rvec / r

  F[i] -= f
  F[j] += f

  E
end

function _Vb!(F, u, i, j, A, b)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = A * exp(-b*r)
  f     = b * E * rvec / r

  F[i] -= f
  F[j] += f

  E
end

function _Vc!(F, u, )

end