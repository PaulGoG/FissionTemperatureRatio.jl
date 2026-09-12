#Level density parameter Egidy-Bucurescu
function levelDensityParameter(A::Int, Z::Int, dm::DataFrame)
    if !isempty(dm.D[(dm.A .== A) .& (dm.Z .== Z)]) && !isempty(dm.D[(dm.A .== A+2) .& (dm.Z .== Z+1)]) && !isempty(dm.D[(dm.A .== A-2) .& (dm.Z .== Z-1)])
        D = only(dm.D[(dm.A .== A) .& (dm.Z .== Z)])
        W_exp = Z *Dᵖ + (A-Z) *Dⁿ - D
        a_sim = A *(27.72 - 25.6 *A^(-1/3))
        η = (A - 2 *Z)/A
        W_LDM = 15.65 *A - 17.63 *A^(2/3) - (0.864/1.233) *Z^2 *A^(-1/3) - a_sim *η^2
        δW₀ = W_LDM - W_exp
        D_plus = only(dm.D[(dm.A .== A+2) .& (dm.Z .== Z+1)])
        D_minus = only(dm.D[(dm.A .== A-2) .& (dm.Z .== Z-1)])
        P_d = (D_plus - 2 *D + D_minus) /4
        δW = δW₀ + P_d
        a = (1.99e-1 + 9.6e-3 *δW) *A^(8.69e-1)
        if a > 0
            return a
        else 
            return NaN
        end
    else
        return NaN
    end
end
#Constructor for the fragmentation domain
function fragmentationDomain(isobaricChgDistrib::DataFrame; A_H_range=AH_min:maximumAH)
    fragmdomain = DataFrame(A = Int[], Z = Int[], value = Float64[])
    for A_H in A_H_range
        A_L = A₀ - A_H
        if !isempty(isobaricChgDistrib.A[isobaricChgDistrib.A .== A_H])
            RMS = only(isobaricChgDistrib.rms[isobaricChgDistrib.A .== A_H])
            ΔZ = only(isobaricChgDistrib.DeltaZ[isobaricChgDistrib.A .== A_H])
        else
            RMS = 0.6
            ΔZ = -0.5
        end
        Zₚ = ΔZ + A_H *Z₀/A₀
        Z_H_min = Int(round(Zₚ) - (noZperA - 1)/2)
        Z_H_max = Z_H_min + noZperA - 1
        for Z_H in Z_H_min:Z_H_max
            push!(fragmdomain, (A_H, Z_H, exp(-(Z_H - Zₚ)^2 /(2 *RMS^2)) /(sqrt(2*π) *RMS)))
            if A_L != first(A_H_range)
                Z_L = Z₀ - Z_H
                push!(fragmdomain, (A_L, Z_L, exp(-(Z_L - Z₀ + Zₚ)^2 /(2 *RMS^2)) /(sqrt(2*π) *RMS)))
            end
        end
    end
    #Sort by mass number in ascending order
    aux_A = sort(fragmdomain.A)
    aux_Z = [sort(fragmdomain.Z[fragmdomain.A .== A]) for A in unique(aux_A)]
    aux_Z = reduce(vcat, aux_Z)
    aux_value = zeros(length(fragmdomain.value))
    for index in eachindex(aux_A)
        aux_value[index] = only(fragmdomain.value[(fragmdomain.A .== aux_A[index]) .& (fragmdomain.Z .== aux_Z[index])])
    end
    fragmdomain.A .= aux_A
    fragmdomain.Z .= aux_Z
    fragmdomain.value .= aux_value
    return fragmdomain
end
#Average q(A,Z) over p(Z,A) so it becomes q(A)
function averageOverZ(q_AZ::DataFrame, fragmdomain::DataFrame)
    q_A = DataFrame(A = Int[], value = Float64[])
    for A in sort(unique(q_AZ.A))
        denominator = 0.0
        numerator = 0.0
        for Z in q_AZ.Z[q_AZ.A .== A]
            Value = only(q_AZ.value[(q_AZ.A .== A) .& (q_AZ.Z .== Z)])
            if !isnan(Value) && !isempty(fragmdomain.value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)])
                P_ZA = only(fragmdomain.value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)])
                denominator += P_ZA
                numerator += P_ZA *Value
            end
        end
        if denominator > 0
            push!(q_A, (A, numerator/denominator))
        end
    end
    return q_A
end
#Constructor for ratioNU(AH) = νH/νPair
function ratioNU(nuA::DataFrame)
    ratio = DataFrame(A = Int[], value = Float64[], σ = Float64[])
    for A_H in AH_min:maximum(nuA.A)
        if !isempty(nuA.nu[nuA.A .== A₀-A_H]) && !isempty(nuA.nu[nuA.A .== A_H])
            νL = only(nuA.nu[nuA.A .== A₀-A_H])
            νH = only(nuA.nu[nuA.A .== A_H])
            R = νL/νH; r = 1/(1 + R)
            if r > 1
                r = 1.0
            end
            if r > 0
                σνL = only(nuA.errnu[nuA.A .== A₀-A_H])
                σνH = only(nuA.errnu[nuA.A .== A_H])
                σR = sqrt(σνL^2 + (R*σνH)^2)/νH; σr = σR/(1 + R)^2
                push!(ratio, (A_H, r, σr))
            end
        end
    end
    return ratio
end
#Constructor for RT(AH)
function ratioOfTemperatures(rν::DataFrame, pZA::DataFrame, dm::DataFrame)
    RT_AH = DataFrame(A = Int[], value = Float64[], σ = Float64[])
    for index in eachindex(rν.A)
        A_H = rν.A[index]
        Z_H = pZA.Z[pZA.A .== A_H]
        levDensP = DataFrame(
            A = [A_H for Z_H in Z_H], 
            Z = [Z_H for Z_H in Z_H], 
            value = [levelDensityParameter(A₀-A_H, Z₀-Z_H, dm)/levelDensityParameter(A_H, Z_H, dm) for Z_H in Z_H]
        )
        Ra = averageOverZ(levDensP, pZA)
        if !isempty(Ra) && only(Ra.value) > 0 && rν.value[index] < 1 && rν.value[index] > 0
            Ra = only(Ra.value)
            RT = sqrt((1-rν.value[index])/(Ra*rν.value[index]))
            σRT = rν.σ[index]/(2*sqrt(Ra)*rν.value[index]) *sqrt(1/(rν.value[index]*(1-rν.value[index])))
            push!(RT_AH, (A_H, RT, σRT))
        end
    end
    return RT_AH
end
#Least square fit for linear function with propagation of errors
function fitLinearLsq(y_init::Number, xData::SubArray, yData::SubArray, σData::SubArray)
    x_init = first(xData)
    linearModel(t, p) = p[1].*(t .- x_init) .+ y_init
    p0 = [1.0]
    weights = 1 ./(σData.^2)
    if isempty(weights[.!isinf.(weights)])
        weights .= 1.0/length(weights)
    else
        weights[isinf.(weights)] .= sqrt(sum(weights[.!isinf.(weights)].^2))/length(weights[.!isinf.(weights)])
        weights ./= sum(weights)
    end
    linearFit = curve_fit(linearModel, xData, yData, weights, p0; autodiff=:forward)
    if linearFit.converged
        slope = only(linearFit.param)
        σSlope = only(stderror(linearFit))
        rSquared = rss(linearFit)
        return slope, σSlope, rSquared
    else
        return NaN, NaN, NaN
    end
end
#Find optimal AH points for linear fits
function findLinFits(fit::FitResult, rν::DataFrame, lineIndex::Int, maxAH::Int)
    if lineIndex < noFitLines
        if lineIndex != 1
            Δx =  fit.objects[lineIndex-1].x_end - fit.objects[lineIndex-1].x_begin
            y_init = fit.objects[lineIndex-1].y_begin + fit.objects[lineIndex-1].slope*Δx
            minAH = fit.objects[lineIndex-1].x_end
        else
            y_init = fit.objects[1].y_begin
            minAH = fit.objects[1].x_begin
        end
        for pivotAH in minAH+1:maxAH
            xData = @view rν.A[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
            yData = @view rν.value[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
            σData = @view rν.σ[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
            if pivotAH - minAH >= minimumLineLength && !isempty(xData)
                slope, σSlope, rSquared = fitLinearLsq(y_init, xData, yData, σData)
                isnan(rSquared) && continue
                auxFit = FitResult([fit.objects[i] for i in eachindex(fit.objects)], 0.0)
                auxFit.objects[lineIndex] = FitObject(minAH, pivotAH, y_init, slope, σSlope, rSquared)
                auxFit.score = sum([auxFit.objects[i].R² for i in eachindex(auxFit.objects)])/length(auxFit.objects)
                auxFit = findLinFits(auxFit, rν, lineIndex+1, maxAH+1)
                if auxFit.score < fit.score
                    if auxFit.objects[lineIndex+1].x_begin == auxFit.objects[lineIndex].x_end
                        if abs(auxFit.objects[lineIndex+1].y_begin - (y_init + (pivotAH - minAH)*slope)) <= 1e-3
                            fit = auxFit
                        end
                    end
                end
            end
        end
    elseif lineIndex == noFitLines
        minAH = fit.objects[noFitLines-1].x_end
        pivotAH = maxAH
        Δx =  fit.objects[noFitLines-1].x_end - fit.objects[noFitLines-1].x_begin
        y_init = fit.objects[noFitLines-1].y_begin + fit.objects[noFitLines-1].slope*Δx
        xData = @view rν.A[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
        yData = @view rν.value[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
        σData = @view rν.σ[(rν.A .>= minAH) .& (rν.A .<= pivotAH)]
        if pivotAH - minAH >= minimumLineLength && !isempty(xData)
            slope, σSlope, rSquared = fitLinearLsq(y_init, xData, yData, σData)
            auxFit = FitResult([fit.objects[i] for i in eachindex(fit.objects)], 0.0)
            auxFit.objects[noFitLines] = FitObject(minAH, pivotAH, y_init, slope, σSlope, rSquared)
            auxFit.score = sum([auxFit.objects[i].R² for i in eachindex(auxFit.objects)])/length(auxFit.objects)
            if auxFit.score < fit.score
                fit = auxFit
            end
        end
    end
    return fit
end