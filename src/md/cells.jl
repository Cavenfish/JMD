#TODO:
#   -Make function to get primitive cells
#   -Make function to clean duplicates

struct Cell <: MyCell
  lattice::Matrix{Float64}
  scaled_pos::Vector{Vector{Float64}}
  masses::Vector{Float64}
  symbols::Vector{Char}
  PBC::Vector{Bool}
  NC::Vector{Int32}
end

function makeCell(bdys::Vector{MyAtoms}, lattice; PBC=repeat([true], 3), NC=[1,1,1])

  Cell(
    lattice,
    getScaledPos(bdys, lattice),
    [i.m for i in bdys],
    [i.s for i in bdys],
    PBC, NC
  )
end

function makeBdys(cell)::Vector{MyAtoms}
  pos  = getPos(cell)
  vel  = [zeros(3) for i in pos]

  [Atom(pos[i], vel[i], cell.masses[i], cell.symbols[i]) for i = 1:length(pos)]
end

function getScaledPos(bdys::Vector{MyAtoms}, lattice)

  T    = inv(lattice)
  
  [T * i.r for i in bdys]
end

function getScaledPos(x0, lattice)

  T    = inv(lattice)
  
  [T * x0[i:i+2] for i = 1:3:length(x0)-1]
end

function getPos(cell)
  [cell.lattice * i for i in cell.scaled_pos]
end

function getVolume(cell)
  a,b,c  = eachrow(cell.lattice)

  cross(b,c) |> (x -> dot(a,x))
end 

function wrap!(cell)
  r = getPos(cell)
  
  f = [floor.(transpose(cell.lattice) \ i) for i in r]

  for i = 1:length(r)
    r[i] .-= transpose(cell.lattice) * f[i]
  end

  T = inv(cell.lattice)

  for i = 1:length(r)
    cell.scaled_pos[i] .= T * r[i]
  end

end

function wrap!(bdys, lattice)

  f = [floor.(transpose(lattice) \ i.r) for i in bdys]
  
  for i = 1:length(bdys)
    bdys[i].r .-= transpose(lattice) * f[i]
  end

end

# Only works for orthogonal cells 
#  TODO: make work for all cells
function center!(cell)
  #Maybe a way to consolodate the following segment would be nice?
  amax = [i[1] for i in cell.scaled_pos] |> maximum
  amin = [i[1] for i in cell.scaled_pos] |> minimum
  bmax = [i[2] for i in cell.scaled_pos] |> maximum
  bmin = [i[2] for i in cell.scaled_pos] |> minimum
  cmax = [i[3] for i in cell.scaled_pos] |> maximum
  cmin = [i[3] for i in cell.scaled_pos] |> minimum
  
  A = (1 - amax - amin) / 2
  B = (1 - bmax - bmin) / 2
  C = (1 - cmax - cmin) / 2

  for i in cell.scaled_pos
    i .+= [A, B, C]
  end

end

function replicate(cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  n      = prod(N)
  newS   = repeat(cell.symbols, n)
  newM   = repeat(cell.masses, n)
  newLat = cell.lattice * Diagonal(N)
  newPos = repeat(getPos(cell), n)

  # Is this style easier to read than inline?
  f = [i*a + j*b + k*c 
        for i = 0:N[1]-1 
          for j = 0:N[2]-1 
            for k = 0:N[3]-1 
              for q = 1:m]
  
  newPos .+= f

  newScaledPos = [inv(newLat) * r for r in newPos]

  Cell(newLat, newScaledPos, newM, newS, cell.PBC, cell.NC)
end

function getMIC(bdys, lattice)
  a,b,c  = eachrow(lattice)
  new    = MyAtoms[]
  s      = repeat([i.s for i in bdys], 27)
  m      = repeat([i.m for i in bdys], 27)

  # I think it is
  f = [i*a + j*b + k*c + bdys[q].r
        for i = -1:1 
          for j = -1:1 
            for k = -1:1 
              for q = 1:length(bdys)]

  for i = 1:length(f)
    push!(new, Atom(f[i], zeros(3), m[i], s[i]))
  end

  new
end

function makeSuperCell(cell, T)

  lattice = T * cell.lattice
  N       = [T[1,1], T[2,2], T[3,3]]
  
  super   = replicate(cell, N)
  bdys    = makeBdys(super)
  wrap!(bdys, lattice)

  makeCell(bdys, lattice, PBC=cell.PBC, NC=cell.NC)
end

function getPrimitiveCell(cell, symprec)
  amu  = TOML.parsefile(joinpath(@__DIR__, "../data/Atoms.toml"))["Mass"]
  num  = TOML.parsefile(joinpath(@__DIR__, "../data/Atoms.toml"))["Number"]
  rnum = Dict(values(num) .=> keys(num))
  atms = string.(cell.symbols) |> (x -> [num[i] for i in x])
  
  L    = transpose(cell.lattice)
  tmp  = Spglib.Cell(L, cell.scaled_pos, atms)
  stan = standardize_cell(tmp, symprec)
 
  syms = [rnum[i] for i in stan.atoms]
  mas  = [ amu[i] for i in syms]
  L    = transpose(stan.lattice)

  Cell(L, stan.positions, mas, only.(syms), cell.PBC, cell.NC)
end

function getMols(cell::MyCell, rmax; D=3)
  r   = getPos(cell)

  pts = hcat(r...)

  ret = dbscan(pts[1:D, :], rmax)

  [i.core_indices for i in ret.clusters]
end

function getPairs(cell::MyCell)

  n = length(cell.masses)

  # Get mols and N
  mols = if n <= 3
    getMols(cell, 1.5, D=n-1) 
  else
    getMols(cell, 1.5)
  end
  N    = size(mols)[1]

  # Make all pairs
  pars = Pair[]
  for i in 1:N
    for j in i+1:N
      push!(pars, Pair(mols[i],mols[j]))
    end
  end

  pars, mols
end