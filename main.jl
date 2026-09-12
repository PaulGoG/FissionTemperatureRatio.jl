include("preamble.jl")
include("functionBodies.jl")

if isdir("inputData/$(inputDir)/")
    datafileList = readdir("inputData/$(inputDir)/")
else
    error("Data dir at PATH inputData/$(inputDir)/ does not exist!")
end
if noFitLines < 2
    error("Minimum number of fit lines is 2!")
end

pZA = fragmentationDomain(isobaricChgDistrib)
nuAPlot = scatter(
    title=inputDir, xlabel="A", ylabel="Prompt neutron multiplicity",
    minorgrid=true, framestyle=:box, ylims=(0,5)
)
rνPlot = scatter(
    title=inputDir, xlabel=L"A_H", ylabel=L"r_\nu =\nu_H/\nu_L + \nu_H",
    minorgrid=true, framestyle=:box, xlims=(AH_min,maximumAH),
    ylims=(0,1)
); hline!(rνPlot, [0.5], linestyle=:dashdot, linecolor=:black, label="")
rνFit = DataFrame(A = Int[], value = Float64[], σ = Float64[])
ratioTemperaturesPlot = scatter(
    title=inputDir, xlabel=L"A_H", ylabel=L"R_T = T_L/T_H",
    minorgrid=true, framestyle=:box, xlims=(AH_min,maximumAH)
)
#Parse nu(A) data in directory and calculate ratios rν & RT
for index in eachindex(datafileList)
    nuA = CSV.read("inputData/$(inputDir)/$(datafileList[index])", DataFrame; 
        header=["A", "nu", "errnu"], skipto=2, ignorerepeated=true, delim=' ', silencewarnings=true
    )
    legendLabel = replace(datafileList[index], inputDir => "")[begin+1:end-4]
    legendLabel = replace(legendLabel, '_' => ' ')
    if !isempty(skipmissing(nuA.errnu))
        Ind = Int[]
        for indices in first(axes(nuA))
            nuRelVariation = abs(nuA.errnu[indices] - nuA.nu[indices])/abs(nuA.nu[indices])
            if !ismissing(nuA.errnu[indices]) && (nuA.errnu[indices] >= nuA.nu[indices] || nuRelVariation <= 0.7 || nuRelVariation >= 0.996)
                push!(Ind, indices)
            elseif ismissing(nuA.errnu[indices])
                nuA.errnu[indices] = 0.0
            end
        end
        deleteat!(nuA, Ind)
    else
        nuA.errnu = [0.0 for i in eachindex(nuA.nu)]
    end
    scatter!(nuAPlot, nuA.A, nuA.nu, yerror=nuA.errnu, markershape=:xcross, label=legendLabel)
    rν = ratioNU(nuA)
    scatter!(rνPlot, rν.A, rν.value, yerror=rν.σ, markershape=:xcross, label=legendLabel)
    RT = ratioOfTemperatures(rν, pZA, dataMassExcess)
    scatter!(ratioTemperaturesPlot, RT.A, RT.value, yerror=RT.σ, markershape=:xcross, label=legendLabel)
    append!(rνFit, rν)
    RT.value .= round.(RT.value, digits=noDigits)
    RT.σ .= round.(RT.σ, digits=noDigits)
    CSV.write(
        "outputData/$(inputDir[begin:end-4])/RTAH/RT_$(datafileList[index])",
        RT, writeheader=true, newline="\r\n", delim=' ', header=["A", "value", "erro"]
    )
end
savefig(nuAPlot, "outputData/$(inputDir[begin:end-4])/nuA.png")

#Compute FitObject for rν and then plot
wt = 1 ./(rνFit.σ.^2)
if isempty(wt[.!isinf.(wt)])
    wt .= 1/length(wt)
else
    wt[isinf.(wt)] .= sqrt(sum(wt[.!isinf.(wt)].^2))/length(wt[.!isinf.(wt)]); wt ./= sum(wt)
end
if iseven(A₀)
    initial_rν = 0.5
else
    initial_rν = sum(rνFit.value[rνFit.A .== AH_min] .*wt[rνFit.A .== AH_min])
end
initialFitObject = FitObject[]
for i in 0:noFitLines-2
    push!(initialFitObject, FitObject(AH_min+i, AH_min+i+1, initial_rν, 0.0, 0.0, length(rνFit.A)))
end
push!(initialFitObject, FitObject(AH_min+noFitLines-1, last(rνFit.A), initial_rν, 0.0, 0.0, length(rνFit.A)))
fitResult = FitResult(initialFitObject, length(rνFit.A)*noFitLines)
fitResult = findLinFits(fitResult, rνFit, 1, maximum(rνFit.A)-noFitLines+1)
rνFitted = DataFrame(A = Int[], value = Float64[], σ = Float64[])
for i in eachindex(fitResult.objects)
    for A in fitResult.objects[i].x_begin:fitResult.objects[i].x_end
        Value = fitResult.objects[i].slope*(A - fitResult.objects[i].x_begin) + fitResult.objects[i].y_begin
        σ = fitResult.objects[i].σSlope
        if Value < 1
            push!(rνFitted, (A, Value, σ))
        end
    end
end
plot!(rνPlot, rνFitted.A, rνFitted.value, ribbon=rνFitted.σ, label="", linewidth=2)
rνAverage = round(sum(rνFit.value .*wt), digits=3)
σrνAverage = round(sqrt(sum(rνFit.σ .^2))/length(rνFit.σ), digits=5)
annotate!(rνPlot, maximumAH-10, 0.05, 
    latexstring("<r_\\nu> = $(rνAverage) \\pm $(σrνAverage)")
)
savefig(rνPlot, "outputData/$(inputDir[begin:end-4])/rNu_$(noFitLines)Segments.png")
CSV.write(
    "outputData/$(inputDir[begin:end-4])/rNu_$(noFitLines)Segments.dat",
    DataFrame(
        AH_begin = [fitResult.objects[i].x_begin for i in eachindex(fitResult.objects)],
        AH_end = [fitResult.objects[i].x_end for i in eachindex(fitResult.objects)],
        rNu_begin = round.([fitResult.objects[i].y_begin for i in eachindex(fitResult.objects)], digits=noDigits),
        rNu_end = round.([(fitResult.objects[i].x_end-fitResult.objects[i].x_begin)*fitResult.objects[i].slope+fitResult.objects[i].y_begin for i in eachindex(fitResult.objects)], digits=noDigits),
        rNu_err = round.([fitResult.objects[i].σSlope for i in eachindex(fitResult.objects)], digits=noDigits)
    ),
    writeheader=true, newline="\r\n", delim=' '
)
#Compute FitObject for RT from calculated rν and then plot
RTFit = ratioOfTemperatures(rνFit, pZA, dataMassExcess)
wt = 1 ./(RTFit.σ.^2)
if isempty(wt[.!isinf.(wt)])
    wt .= 1/length(wt)
else
    wt[isinf.(wt)] .= sqrt(sum(wt[.!isinf.(wt)].^2))/length(wt[.!isinf.(wt)]); wt ./= sum(wt)
end
averageRT = round(sum(RTFit.value .*wt), digits=3)
averageRTσ = round(sqrt(sum(RTFit.σ .^2))/length(RTFit.σ), digits=5)
if iseven(A₀)
    initial_RT = 1.0
else
    initial_RT = sum(RTFit.value[RTFit.A .== AH_min] .*wt[RTFit.A .== AH_min])
end
RTFit = ratioOfTemperatures(rνFitted, pZA, dataMassExcess)
RTFit.σ .= 0.0
initialFitObject = FitObject[]
for i in 0:noFitLines-2
    push!(initialFitObject, FitObject(AH_min+i, AH_min+i+1, initial_RT, 0.0, 0.0, length(RTFit.A)))
end
push!(initialFitObject, FitObject(AH_min+noFitLines-1, last(RTFit.A), initial_RT, 0.0, 0.0, length(RTFit.A)))
fitResult = FitResult(initialFitObject, length(RTFit.A)*noFitLines)
fitResult = findLinFits(fitResult, RTFit, 1, maximum(RTFit.A)-noFitLines+1)
RTFitted = DataFrame(A = Int[], value = Float64[], σ = Float64[])
for i in eachindex(fitResult.objects)
    for A in fitResult.objects[i].x_begin:fitResult.objects[i].x_end
        Value = fitResult.objects[i].slope*(A - fitResult.objects[i].x_begin) + fitResult.objects[i].y_begin
        σ = fitResult.objects[i].σSlope
        push!(RTFitted, (A, Value, σ))
    end
end
hline!(ratioTemperaturesPlot, [averageRT], linestyle=:dashdot, linecolor=:black, label="")
annotate!(ratioTemperaturesPlot, maximumAH-10,averageRT*1.1, latexstring("<R_T> = $(averageRT) \\pm $(averageRTσ)"))
plot!(ratioTemperaturesPlot, RTFitted.A, RTFitted.value, ribbon=RTFitted.σ, label="", linewidth=2)
savefig(ratioTemperaturesPlot, "outputData/$(inputDir[begin:end-4])/RT_$(noFitLines)Segments.png")
CSV.write(
    "outputData/$(inputDir[begin:end-4])/RT_$(noFitLines)Segments.dat",
    DataFrame(
        AH_begin = [fitResult.objects[i].x_begin for i in eachindex(fitResult.objects)],
        AH_end = [fitResult.objects[i].x_end for i in eachindex(fitResult.objects)],
        RT_begin = round.([fitResult.objects[i].y_begin for i in eachindex(fitResult.objects)], digits=noDigits),
        RT_end = round.([(fitResult.objects[i].x_end-fitResult.objects[i].x_begin)*fitResult.objects[i].slope+fitResult.objects[i].y_begin for i in eachindex(fitResult.objects)], digits=noDigits),
        RT_err = round.([fitResult.objects[i].σSlope for i in eachindex(fitResult.objects)], digits=noDigits)
    ),
    writeheader=true, newline="\r\n", delim=' '
)
println("*program execution successful!")