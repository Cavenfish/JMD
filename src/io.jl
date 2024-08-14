export readASExyz, readXyz, writeXyz, writeXyzTraj

function readASExyz(xyz)
  sys = readlines(xyz)
  N   = length(sys) - 2
  hed = sys[2]
  sys = split.(sys[3:end], " ")
  sys = deleteat!.(sys, findall.(e -> e == "", sys))
  amu = TOML.parsefile(joinpath(@__DIR__, "data/Atoms.toml"))["Mass"]
  set = Atom[]

  for i in range(1,N)
    props = parse.(Float64, sys[i][2:end])
    pos   = Vector(props[1:3])
    s     = sys[i][1]

    if occursin("masses", hed)
      mas = props[4]
      vel = Vector(props[5:7] ./ mas)
    else
      mas = amu[s]
      vel = Vector(props[4:6] ./ mas)
    end #if-else

    particle = Atom(pos, vel, mas, s[1])
    push!(set, particle)
  end #for loop

  return set
end #read_ase_xyz

function readXyz(xyz)
  stream = readlines(xyz)
  amu    = TOML.parsefile(joinpath(@__DIR__, "data/Atoms.toml"))["Mass"]
  set    = Atom[]

  #Skip header lines then parse file
  for line in stream[3:end]
    s = split(line, " ")
    s = deleteat!(s, findall(e -> e == "", s))

    pos  = parse.(Float64, s[2:4])
    pos  = Vector(pos)
    
    vel  = if length(s) >= 7
      tmp = parse.(Float64, s[5:end])
      Vector(tmp)
    else
      Vector([0.0,0.0,0.0]) 
    end
    
    mas  = amu[s[1]]
    sym  = s[1][1]

    atom = Atom(pos, vel, mas, sym)
    push!(set, atom)
  end

  return set
end

function writeXyzTraj(fileName::String, solu; dt=1)
  f    = open(fileName, "w")
  bdys = solu.prob.p.bdys
  N    = length(bdys)
  T    = length(solu.t)


  for i in 1:dt:T
    t = solu.t[i]
    u = solu.u[i].x[2] # x[1] -> vel || x[2] -> pos

    println(f, N)
    println(f, "i=$i, time=$t")

    for j in 1:N

      s     = bdys[j].s
      x,y,z = u[j]

      println(f, "$s   $x   $y   $z")
    end 

  end
  close(f)
end 

#Hacky and only works for CO right now
function writeXyzTraj(fileName::String, tj::MyTraj; dt=1)
  f = open(fileName, "w")
  S = map(x -> x > 15 ? 'O' : 'C', tj.m)
  T = length(tj.t)
  N = length(tj.m)

  for i in 1:dt:T
    t = tj.t[i]
    u = tj.r[i]

    println(f, N)
    println(f, "i=$i, time=$t")

    for j in 1:N

      s     = S[j]
      x,y,z = u[j]

      println(f, "$s   $x   $y   $z   $vx   $vy   $vz")
    end 

  end
  close(f)
end 


function writeXyz(fileName::String, bdys)
  f = open(fileName, "w")
  N = length(bdys)

  println(f, N)
  println(f, "Made by JMD")

  for j in 1:N

    s          = bdys[j].s
    x,y,z      = bdys[j].r
    vx, vy, vz = bdys[j].v 

    println(f, "$s   $x   $y   $z   $vx   $vy   $vz")
  end 

  close(f)
end 

