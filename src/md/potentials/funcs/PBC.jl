
function _pbcInteractions!(F, u, a, b, func, L, NC, p)
  E = 0.0
  for i = 1:3
    for j = -NC[i]:NC[i]
      j == 0 && continue
      
      r2    = u[a] + (L[i] * j)
      e,f   = func(u[b], r2, p...)
      E    += e
      F[b] -= f
      
      r2    = u[b] + (L[i] * j)
      e,f   = func(u[a], r2, p...)
      E    += e
      F[a] -= f
    end
  end
  E
end