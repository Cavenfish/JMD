
# This is copied from the ASE implementation 
function getNumericalStress(EoM, cell::MyCell; eps=1e-6)
  stress = zeros(3,3)
  vol    = getVolume(cell)
  tmp    = deepcopy(cell)

  for i = 1:3

    x = Matrix(1.0I, 3,3)

    x[i,i]       = 1.0 + eps
    tmp.lattice .= cell.lattice * x
    eplus        = getPotEnergy(EoM, tmp)

    x[i,i]       = 1.0 - eps
    tmp.lattice .= cell.lattice * x
    eminus       = getPotEnergy(EoM, tmp)

    stress[i,i]  = (eplus - eminus) / (2 * eps * vol)
    x[i,i]       = 1.0

    i+1 > 3 ? j =  (i+1)%3 : j = i+1

    x[i,j] = x[j,i] = +0.5 * eps
    tmp.lattice    .= cell.lattice * x
    eplus           = getPotEnergy(EoM, tmp)

    x[i,j] = x[j,i] = -0.5 * eps
    tmp.lattice    .= cell.lattice * x
    eminus          = getPotEnergy(EoM, tmp)

    stress[i,j] = stress[j,i] = (eplus - eminus) / (2 * eps * vol)
  end

  stress
end

# A version of the above that only operates on diagonal terms of matrix
function getNumericalStressOrthogonal(EoM, cell::MyCell; eps=1e-6)
  stress = zeros(3,3)
  vol    = getVolume(cell)
  tmp    = deepcopy(cell)

  for i = 1:3

    x = Matrix(1.0I, 3,3)

    x[i,i]       = 1.0 + eps
    tmp.lattice .= cell.lattice * x
    eplus        = getPotEnergy(EoM, tmp)

    x[i,i]       = 1.0 - eps
    tmp.lattice .= cell.lattice * x
    eminus       = getPotEnergy(EoM, tmp)

    stress[i,i]  = (eplus - eminus) / (2 * eps * vol)
    x[i,i]       = 1.0

  end

  stress
end