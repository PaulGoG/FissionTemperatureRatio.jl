using CSV, DataFrames, LsqFit, Plots, LaTeXStrings
cd(@__DIR__)
gr(size = (1000, 1000), dpi = 600)
if !isdir("outputData/")
    mkdir("outputData/")
end

#Input data parameters
const inputDir = "U233_nf_nuA"
const targetNucleusID = "U33"
const noZperA = 5
const noDigits = 6
const noFitLines = 5
const minimumLineLength = 5
const A₀, Z₀ = 234, 92
const AH_min = Int(ceil(A₀/2))
const maximumAH = 160

#Read mass excess and isobaric charge distribution of Wahl data from files
dataMassExcess = CSV.read("inputData/AME2020.ANA", DataFrame; 
    delim = ' ', ignorerepeated=true, header = ["Z", "A", "symbol", "D", "errD"]
)
dataMassExcess.D .*= 1e-3; dataMassExcess.errD .*= 1e-3
const Dᵖ = only(dataMassExcess.D[(dataMassExcess.A .== 1) .& (dataMassExcess.Z .== 1)])
const Dⁿ = only(dataMassExcess.D[(dataMassExcess.A .== 1) .& (dataMassExcess.Z .== 0)])

if isfile("inputData/DeltaZ_rms_A.$(targetNucleusID)")
    isobaricChgDistrib = CSV.read("inputData/DeltaZ_rms_A.$(targetNucleusID)", DataFrame; 
        delim = ' ', ignorerepeated=true, header = ["A", "DeltaZ", "rms"], skipto=2
    )
else
    isobaricChgDistrib = DataFrame(A = [])
end

#Define struct object for linear fits on data intervals
abstract type AbstractStruct end
struct FitObject <: AbstractStruct
    x_begin::Int
    x_end::Int
    y_begin::Number
    slope::Number
    σSlope::Number
    R²::Number
end
mutable struct FitResult <: AbstractStruct
    objects::Vector{FitObject}
    score::Number
end

if !isdir("outputData/$(inputDir[begin:end-4])/")
    mkdir("outputData/$(inputDir[begin:end-4])/")
    mkdir("outputData/$(inputDir[begin:end-4])/RTAH/")
end