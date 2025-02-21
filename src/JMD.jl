"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""
module JMD

  using TOML
  using JLD2
  using FFTW
  using Libdl
  using Optim
  using PyCall
  using LsqFit
  using Spglib
  using MiniQhull
  using DataFrames
  using Statistics#might not be used
  using Clustering
  using StaticArrays
  using Serialization
  using LinearAlgebra
  using Distributions#might not be used
  using KernelDensity
  using OrdinaryDiffEq
  using SpecialFunctions#might not be used
  using FiniteDifferences


  export
    #constants
    fs, ps, ns, kB,

    #io.jl
    readASExyz, readXyz, writeXyz, readCell, writeCell, writeXyzTraj,

    #helpers.jl
    CoM, vCoM, zeroVCoM!, swapIso!, vibExcite!, transExcite!, pickRandomMol, 
    getFrame, getFrame!, getLastFrame, getLastFrame!, getPotEnergy, getForces,
    centerBdys!,

    #bodies.jl
    getMols, getPairs,

    #cells.jl
    makeCell, makeBdys, getScaledPos, getPos, wrap!, replicate, makeSuperCell,
    getMIC, center!, getPrimitiveCell, getVolume,

    #potentials
    COCO, HGNN, MBX, SPCF, TIP4P, SCME,

    #simulation.jl
    runNVE, runNVT,

    #thermostats.jl
    Berendsen, Berendsen!, Langevin, Langevin!, BDP, BDP!, BDPnT, BDPnT!,

    #post-processing.jl
    processDynamics, processDynamics!, processTmpFiles, trackVACF,
    trackEnergyDissipation, trackAllVibEnergy, trackRadialEnergy,

    #vacf.jl
    vacfInps, VDOS, getVelMas,

    #desorb.jl

    #vibrations.jl
    getHarmonicFreqs, animateMode, getModePES, getModeInteractionPES,

    #optimizations.jl
    opt, optCell,

    #structural.jl
    rdf, adf, density,

    #stress.jl
    getNumericalStress, getNumericalStressOrthogonal,

    #freqShifts.jl
    getINM, getMolFreq, getAllFreqs, getFvE, getFreqCoupling,

    #participationRatio.jl
    getIPR, getPR,

    #neighbors.jl
    countNearestNeighbors,

    #vibCoup.jl
    getVibCoup,

    #alphashape.jl
    alphashape,

    #savitzkyGolay.jl
    savGol,

    #peakFinding.jl
    findPeaks, findTurningPoints,

    #anneal.jl
    anneal,

    #hitAndStick.jl
    hitAndStick, HnS,

    #phonopy.jl
    phonopy_addForces
  
  #end exports

  const fs   = 0.09822694788464063 # 1fs in ASE time
  const ps   = 1000 * fs
  const ns   = 1000 * ps
  const kB   = 8.617333262e-5 # eV / K

  include("./types.jl")

  include("./io.jl")
  include("./macros.jl")
  include("./helpers.jl")

  include("./lib/MBX/libmbx.jl")
  include("./lib/SCME/libscme.jl")
  include("./lib/Phonopy/phonopy.jl")

  include("./md/cells.jl")
  include("./md/bodies.jl")
  include("./md/potentials/MvHffCO.jl")
  include("./md/potentials/COCOff.jl")
  include("./md/potentials/HGNN.jl")
  include("./md/potentials/TIP4P.jl")
  include("./md/potentials/SPC-F.jl")
  include("./md/potentials/CH4.jl")
  include("./md/potentials/MBX.jl")
  include("./md/potentials/SCME.jl")
  include("./md/simulation.jl")
  include("./md/thermostats.jl")
  include("./md/post-processing.jl")
  include("./md/potentials/funcs/PBC.jl")
  include("./md/potentials/funcs/intra.jl")
  include("./md/potentials/funcs/inter.jl")
  include("./md/potentials/funcs/damping.jl")
  include("./md/potentials/funcs/TTMnrg.jl")

  include("./analysis/vacf.jl")
  include("./analysis/desorb.jl")
  include("./analysis/vibrations.jl")
  include("./analysis/optimizations.jl")
  include("./analysis/structral.jl")
  include("./analysis/decayRates.jl")
  include("./analysis/freqShifts.jl")
  include("./analysis/participationRatio.jl")
  include("./analysis/neighbors.jl")
  include("./analysis/vibCoup.jl")

  include("./mathtk/alphashape.jl")
  include("./mathtk/savitzkyGolay.jl")
  include("./mathtk/peakFinding.jl")
  include("./mathtk/stress.jl")

  include("./building/anneal.jl")
  include("./building/hitAndStick.jl")

  include("./QM/orcaConfig.jl")
end # module
