
# Communication between phonopy and JMD goes through the phonopy yaml file.
# Python should be used for all non-force calculations.

function phonopy_addForces(yamlFile, saveName, forces; symprec=1e-5)
  phonopy    = pyimport("phonopy")
  obj        = phonopy.load(yamlFile, symprec=symprec)
  obj.forces = forces
  obj.save(saveName)
end