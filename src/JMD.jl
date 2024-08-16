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
using LsqFit
using CodecZlib
using MiniQhull
using DataFrames
using Statistics
using CairoMakie
using Clustering
using StaticArrays
using Serialization
using LinearAlgebra
using Distributions
using KernelDensity
using OrdinaryDiffEq
using SpecialFunctions
using FiniteDifferences

const fs   = 0.09822694788464063 # 1fs in ASE time
const ps   = 1000 * fs
const ns   = 1000 * ps
const kB   = 8.617333262e-5 # eV / K

include("./types.jl")

include("./io.jl")
include("./macros.jl")
include("./helpers.jl")

include("./lib/MBX/libmbx.jl")

include("./md/bodies.jl")
include("./md/potentials/MvHffCO.jl")
include("./md/potentials/COCOff.jl")
include("./md/potentials/HGNN.jl")
include("./md/potentials/TIP4P.jl")
include("./md/potentials/SPC-F.jl")
include("./md/potentials/CH4.jl")
include("./md/potentials/MBX.jl")
include("./md/simulation.jl")
include("./md/thermostats.jl")
include("./md/post-processing.jl")
include("./md/potentials/funcs/intra.jl")
include("./md/potentials/funcs/inter.jl")
include("./md/potentials/funcs/damping.jl")
include("./md/potentials/funcs/TTMnrg.jl")
include("./md/potentials/funcs/NASApot.jl")

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

include("./plotting/config.jl")
include("./plotting/vibDisp.jl")
include("./plotting/clusters.jl")

include("./mathtk/alphashape.jl")
include("./mathtk/savitzkyGolay.jl")
include("./mathtk/peakFinding.jl")

include("./building/anneal.jl")
include("./building/hitAndStick.jl")

include("./experiments/localHeating.jl")
include("./experiments/bindingEnergy.jl")
include("./experiments/energyDissipation.jl")

include("./QM/orcaConfig.jl")

include("./organization/dataManagement.jl")
end # module
