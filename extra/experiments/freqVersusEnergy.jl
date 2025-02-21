using TOML, JLD2

"""
Calculates the vibrational frequency at various
vibrational energies.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
clus = "/home/brian/myStuff/clusters.jld2"
json = "/home/brian/myStuff/jsonFiles/mol.json"
modes = "/home/brian/myStuff/modes.jld2"
Erange = [0.01, 0.01, 1.5]


[Saving]
data = "myData.jld2"
"""
function fve(inpFile::String)

  # Read input card
  inp = TOML.parsefile(inpFile)
  cfg = inp["Settings"]
  sav = inp["Saving"]["data"]

  # Load clusters and modes
  clus  = jldopen(cfg["clus"])
  modes = jldopen(cfg["modes"])

  # Load EoM 
  EoM = JMD.mkvar(cfg["EoM"])

  # Make energy range
  a = cfg["Erange"]
  E = collect(a[1]:a[2]:a[3])

  # Save energy range and init saving file
  jldsave(sav; E)

  for k1 in keys(clus)

    for k2 in keys(clus[k1])

      bdys = clus[k1][k2]
      
      # Open data file to write data in
      # I close between bdys incase one crashes you have some data still
      jldopen(sav, "a+") do file

        for k3 in keys(modes[k1][k2])

          m = modes[k1][k2][k3]
          ν = getFvE(MBX, bdys, E, m)

          file["$(k1)/$(k2)/$(k3)"] = ν

        end# k3

      end# close file

    end# k2

  end# k1


end