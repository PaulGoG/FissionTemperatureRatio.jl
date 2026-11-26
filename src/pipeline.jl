# The pipeline: from experimental multiplicity data to a segmented temperature ratio.

"""
    run_pipeline(configuration; leave_one_out = true) -> ExtractionResult

Execute the extraction for one configuration.

The steps are: build the fragmentation domain from the mass table and the isobaric charge
distribution, and the level density model; read the experimental multiplicity data; form the
multiplicity ratio `r_ν = ν_H/(ν_L + ν_H)` of every dataset and the temperature ratio from it; and
describe the multiplicity ratio by joined straight segments, once per dataset and once more
following its systematic behaviour.

A run therefore returns several segmented curves. Where the datasets of a fissioning nucleus
disagree, as they do markedly for some, a curve fitted through all of them describes none of them,
and a prompt emission code has nothing to validate it against. They are
offered as alternatives: the code takes one as input, and which one describes reality is settled
by comparing the multiplicity distributions and yields it produces against experiment. The
systematic-trend curve is the one to use where a dataset is too sparse or too scattered to
determine the shape on its own.

The temperature ratio is obtained by inverting the fitted multiplicity ratio, not by fitting
segments to the temperature ratio a second time. The inversion is exact, whereas a second fit
would discard the uncertainty of the first and impose a piecewise-linear shape on a quantity that
is not piecewise-linear.

`leave_one_out` asks for the refits of the systematic trend with pooled datasets left out, the
datasets of one experiment together: one with each left out in turn, [`LeaveOneOut`](@ref), and
those with half of them left out behind its jackknife uncertainty, [`JackknifeDesign`](@ref).
With at least two of the latter that gave a curve the trend carries their jackknife: its
uncertainties are the larger of the propagated one and the jackknife. The refits cost one trend
fit per pooled dataset or experiment and at most `segments.jackknife_subsets` more, and can be
skipped where only the curves are wanted; the result then carries none, records that they were
not asked for, and the trend carries the propagated uncertainty alone.

Nothing is written. [`write_results`](@ref) writes the tables and the manifest of a result into a
directory of the caller's choosing; `scripts/run.jl` names that directory by the run identifier
and adds the provenance record and the figures.
"""
function run_pipeline(configuration::Configuration; leave_one_out::Bool = true)
    @info "reading input" configuration = configuration.source
    retrievals = _read_retrievals(_retrieval_directories(configuration))
    inversion = _inversion(configuration, retrievals)
    relation = (inversion.averaging, inversion.model, inversion.domain)
    input = _multiplicity_input(configuration, inversion, retrievals)
    per_dataset = _dataset_curves(configuration, relation, input)
    pool = _pooled_trend(configuration, relation, input)
    trend = pool.fitted.trend
    isempty(per_dataset.curves) &&
        !(trend isa ExtractedCurve) &&
        throw(
            ArgumentError(
                "no dataset supports a segmented curve with min_points_per_segment = \
             $(configuration.segments.min_points_per_segment)"
            ),
        )
    yields = _mass_yield_input(configuration, inversion, retrievals)
    refits = if leave_one_out && trend isa ExtractedCurve
        _leave_one_out(pool.units, input.datasets, relation, configuration, yields.averaged)
    else
        (; entries = LeaveOneOut[], curves = ExtractedCurve[], windows = Bool[])
    end
    # The trend takes the jackknife over its refits before anything reads its uncertainty.
    jackknife = nothing
    if leave_one_out && trend isa ExtractedCurve
        resampled = _jackknife_refits(
            pool.units,
            relation,
            configuration,
            trend,
            pool.fitted.windows,
            refits.curves,
            refits.windows,
        )
        if resampled !== nothing && resampled.design.refits ≥ 2
            design = resampled.design
            propagated = trend.R_T.σ
            trend = ExtractedCurve(
                trend, resampled.curves; units = design.units, deleted = design.deleted
            )
            jackknife = design
            @info "the uncertainty of the systematic trend is the larger of the propagated one \
                   and the jackknife over the refits with pooled datasets left out" units =
                design.units deleted = design.deleted refits = design.refits masses = length(
                trend.R_T
            ) raised = count(trend.R_T.σ .> propagated)
        end
    end
    segmented_curves = copy(per_dataset.curves)
    trend isa ExtractedCurve && push!(segmented_curves, trend)
    for curve in segmented_curves
        @info "segmented curve" dataset = curve.label segments = segment_count(curve.fit) breakpoints =
            curve.fit.breakpoints reduced_chi_squared = curve.fit.wrss / curve.fit.dof
    end
    symmetry = _symmetry_diagnostics(configuration, relation, segmented_curves)
    for message in symmetry.warnings
        @warn message
    end

    heavy_masses = inversion.domain.heavy_masses
    total_average_R_T = _total_averages(segmented_curves, yields.averaged, heavy_masses)
    for curve in segmented_curves, distribution in yields.averaged
        entry = get(get(total_average_R_T, curve.label, Dict()), distribution.label, nothing)
        entry === nothing && continue
        @info "total average" dataset = curve.label yield = distribution.label R_T = entry.value uncertainty =
            entry.uncertainty
    end
    sensitivity = if trend isa ExtractedCurve
        _segment_count_sensitivity(trend, pool.fitted, relation, configuration, yields.averaged)
    else
        Dict{String,Vector{Tuple{Int,Float64}}}()
    end

    datasets = input.datasets
    return ExtractionResult(
        configuration,
        datasets,
        input.r_ν,
        input.R_T,
        inversion.masses,
        inversion.model,
        inversion.charge,
        inversion.domain,
        inversion.averaging,
        inversion.mean_kinetic_energy,
        segmented_curves,
        Dict(c.label => range_mean(c) for c in segmented_curves),
        yields.mass_yields,
        yields.records,
        yields.coverage,
        total_average_R_T,
        pool.fitted.pooled.curve,
        input.diagnostics,
        per_dataset.outcomes,
        symmetry,
        sensitivity,
        pool.fitted.correlogram,
        pool.fitted.autocorrelation,
        Dict(datasets[i].label => datasets[j].label for (i, j) in pool.superseded),
        [
            [datasets[i].label for i in members] for
            members in pool.units.members if length(members) > 1
        ],
        leave_one_out,
        refits.entries,
        jackknife,
        per_dataset.flags,
        retrievals,
    )
end

# The fragmentation domain and the relation between R_T and E*_H/TXE on it: the mass table, the
# level density model, the charge distribution, the ⟨TKE⟩(A) dataset where one is configured, and
# the averaging that reduces the relation over the charges.
function _inversion(configuration::Configuration, retrievals::RetrievalDocuments)
    settings = configuration.level_density
    masses = build_mass_table(settings)
    model = build_level_density_model(settings, masses)
    charge = build_charge_model(configuration, masses)
    range = A_H_range(configuration)
    domain = fragmentation_domain(
        configuration.system,
        charge,
        range;
        charges_per_mass = configuration.fragmentation.charges_per_mass,
        zero_polarization_at_symmetry = configuration.fragmentation.zero_polarization_at_symmetry,
    )
    @info "fragmentation domain" domain charge_model = charge_model_label(charge) level_density =
        settings.model
    mean_kinetic_energy = if settings.mean_kinetic_energy_file === nothing
        nothing
    else
        file = settings.mean_kinetic_energy_file
        records = retrieval_records(_documents(retrievals, dirname(file)))
        read_mean_kinetic_energy(
            file,
            configuration.system,
            range;
            record = get(records, basename(file), nothing),
        )
    end
    # As `load_configuration` requires, for a configuration built by hand: the dataset's
    # provenance and qualifiers are on record, or it weights nothing.
    mean_kinetic_energy === nothing ||
        mean_kinetic_energy.record !== nothing ||
        throw(
            ArgumentError("level_density.mean_kinetic_energy_file: no retrieval record beside \
                 $(settings.mean_kinetic_energy_file) lists it")
        )
    averaging = _averaging(settings.ratio_averaging, masses, domain, mean_kinetic_energy)
    return (; masses, model, charge, domain, mean_kinetic_energy, averaging)
end

# The multiplicity datasets and what is formed from each before anything is fitted: r_ν and R_T
# point by point, and the diagnostics. The exclusions are checked against the datasets read, and
# the qualifiers and scales the retrieval records are reported, first.
function _multiplicity_input(
    configuration::Configuration, inversion, retrievals::RetrievalDocuments
)
    directory = configuration.multiplicity_directory
    datasets = read_multiplicity_directory(
        directory; records = retrieval_records(_documents(retrievals, directory))
    )
    # An exclusion that names no dataset is a configuration the run cannot honour.
    # `load_configuration` has checked it; a configuration built by hand has not.
    _excluded_labels(configuration.excluded_datasets, datasets)
    if configuration.yield_directory !== nothing
        held = _directory_keys(
            configuration.yield_directory,
            retrieval_records(_documents(retrievals, configuration.yield_directory)),
        )
        for key in sort!(collect(keys(configuration.excluded_mass_yields)))
            haskey(held, key) || throw(
                ArgumentError(
                    "yield.exclude excludes $(repr(key)), which names no distribution of \
                     $(configuration.yield_directory); an exclusion is keyed by the EXFOR \
                     accession of its dataset"
                ),
            )
        end
    end
    for data in datasets
        flagged = _flagged(_qualifiers(data.record))
        isempty(flagged) ||
            @warn "the retrieval records a qualifier on this dataset; it is used, not \
                   corrected" dataset = data.label qualifiers = flagged
        scale = _pair_sum(data.record)
        scale !== nothing &&
            scale.consistent === false &&
            @warn "the retrieval records the pair sum of this dataset off the scale of ν̄; a \
                   uniform scale cancels in r_ν, so it is used, not corrected" dataset =
                data.label deviation = scale.deviation uncertainty = scale.uncertainty
    end

    A₀ = configuration.system.compound.A
    range = A_H_range(configuration)
    r_ν = [multiplicity_ratio(data, A₀, range) for data in datasets]
    R_T = [
        temperature_ratio(inversion.averaging, inversion.model, inversion.domain, curve) for
        curve in r_ν
    ]
    diagnostics = [diagnose(datasets[i], r_ν[i], A₀, range) for i in eachindex(datasets)]
    any(!isempty, r_ν) || throw(
        ArgumentError(
            "no dataset provides both fragments of any pair within $(first(range)):$(last(range))",
        ),
    )
    return (; datasets, r_ν, R_T, diagnostics)
end

# One curve per dataset: the alternatives a prompt emission code chooses between. A dataset below
# the coverage floor is read, diagnosed and pooled, but offers no curve of its own: the breakpoint
# search cannot place the minimum where the data do not reach. Returned with what became of every
# dataset and, for a curve that resolves no minimum where its range covers a required window, the
# reason; such a curve is offered all the same.
function _dataset_curves(configuration::Configuration, relation, input)
    settings = configuration.segments
    A₀ = configuration.system.compound.A
    range = A_H_range(configuration)
    windows = settings.windows_apply_to_datasets ? settings.required_windows : UnitRange{Int}[]
    curves = ExtractedCurve[]
    outcomes = Dict{String,String}()
    for (index, data) in enumerate(input.datasets)
        curve = input.r_ν[index]
        if isempty(curve)
            outcomes[data.label] = "no complete fragment pair within $(first(range)):$(last(range))"
            continue
        end
        coverage = input.diagnostics[index].coverage
        if coverage < settings.min_pair_coverage
            outcomes[data.label] = "coverage $(round(coverage; digits = 3)) below the floor \
                                    $(settings.min_pair_coverage); pooled only"
            @info "no segmented curve for this dataset" dataset = data.label coverage floor =
                settings.min_pair_coverage
            continue
        end
        # Counted in measurements: an interpolated dataset's points count for raw over written.
        measured = fill(_measured_fraction(data.record), length(curve))
        segmented = _segment(curve, relation, settings, curve.label, windows, A₀; measured)
        if segmented isa ExtractedCurve
            push!(curves, segmented)
            outcomes[data.label] = "segmented curve"
        else
            outcomes[data.label] = "no fit: $(segmented.reason)"
        end
    end
    unpinned = [c.label for c in curves if c.fit.pinned_value === nothing]
    settings.pin_symmetric_split &&
        !isempty(unpinned) &&
        @info "fitted without the pin: no complete pair at the symmetric split" datasets =
            unpinned
    flags = Dict{String,String}()
    for c in curves
        reason = unresolved_minimum(c.fit, settings.required_windows)
        reason === nothing && continue
        flags[c.label] = reason
        @info "dataset curve flagged" dataset = c.label reason
    end
    return (; curves, outcomes, flags)
end

# The systematic trend: the minimum at the heavy magic fragment placed by the required windows,
# the rise above the most probable fragmentation fitted through the pooled data. Datasets the
# configuration excludes stay out of the pool but are still fitted, written and diagnosed; of one
# result published twice the superseding dataset alone enters, and the datasets of one experiment
# are combined into one before they are pooled.
function _pooled_trend(configuration::Configuration, relation, input)
    datasets = input.datasets
    candidates = _candidates(configuration, datasets, input.r_ν)
    isempty(candidates) && throw(
        ArgumentError("every usable dataset is excluded from the pooling by configuration")
    )
    superseded = _superseded(datasets, candidates)
    admitted = [i for i in candidates if !haskey(superseded, i)]
    for (i, j) in superseded
        if j in candidates
            @info "dataset superseded by a republication: pooled through its successor alone" dataset =
                datasets[i].label superseded_by = datasets[j].label
        else
            @warn "dataset superseded by a republication that is not pooled itself; neither \
                   enters the pool" dataset = datasets[i].label superseded_by =
                datasets[j].label
        end
    end
    units = _pool_units(datasets, input.r_ν, admitted)
    for members in units.members
        length(members) > 1 && @info "datasets of one experiment pooled as one" datasets = [
            datasets[i].label for i in members
        ]
    end
    # Interpolated datasets count by their measured points, not their written rows: on its own
    # at its measured fraction, combined with others of its experiment through the measured
    # fraction of each point of the experiment's curve.
    for members in units.members, i in members
        fraction = _measured_fraction(datasets[i].record)
        fraction < 1 || continue
        if length(members) == 1
            @info "interpolated dataset pooled at reduced weight" dataset = datasets[i].label weight =
                fraction
        else
            @info "interpolated dataset combined into the curve of its experiment, which carries \
                   the measured fraction of each point" dataset = datasets[i].label fraction
        end
    end
    fitted = _fit_trend(
        units.curves,
        units.weights,
        relation,
        configuration.segments,
        configuration.system.compound.A,
    )
    return (; superseded, units, fitted)
end

# The mass yield distributions: the primary one of the system, averaged over alone or, beside a
# directory of distributions, the reference their coverage is measured against; the retrieval
# record of each; the coverage of each in yield, against the primary, since a set's own yields
# cannot measure it; and those averaged over, the others being below the coverage floor or
# excluded by the configuration, reported but not averaged over.
function _mass_yield_input(
    configuration::Configuration, inversion, retrievals::RetrievalDocuments
)
    A₀ = configuration.system.compound.A
    reference = if configuration.yield_file === nothing
        nothing
    else
        read_mass_yield(
            configuration.yield_file;
            label = _dataset_label(basename(configuration.yield_file)),
        )
    end
    mass_yields = if configuration.yield_directory !== nothing
        read_mass_yield_directory(
            configuration.yield_directory;
            records = retrieval_records(_documents(retrievals, configuration.yield_directory)),
        )
    elseif reference !== nothing
        [reference]
    else
        MassYield[]
    end
    # Read before the symmetrization, which keeps the source and the label of each.
    records = _records_by_label(mass_yields, retrievals)
    if configuration.symmetrize_yields && !isempty(mass_yields)
        mass_yields = [symmetrized_yield(y, A₀) for y in mass_yields]
        reference = symmetrized_yield(reference, A₀)
        @info "mass yields symmetrized: Y(A) and Y(A₀ - A) averaged where both are measured, \
               one wing standing for the other where it alone is"
    end
    for distribution in mass_yields
        flagged = _flagged(_qualifiers(get(records, distribution.label, nothing)))
        isempty(flagged) ||
            @warn "the retrieval records a qualifier on this yield distribution; it is used, \
                   not corrected" distribution = distribution.label qualifiers = flagged
    end
    energies = inversion.mean_kinetic_energy
    if energies !== nothing
        for distribution in mass_yields
            offset = mean_kinetic_energy_offset(energies, distribution, configuration.system)
            @info "⟨TKE⟩ of the excitation weights against the energy standard" dataset =
                energies.label yield = distribution.label mean_MeV = offset.mean standard_MeV =
                offset.standard offset_MeV = offset.offset
        end
    end
    heavy_masses = inversion.domain.heavy_masses
    coverage = Dict(
        y.label => mass_yield_coverage(y, reference, heavy_masses) for y in mass_yields
    )
    excluded = if configuration.yield_directory === nothing
        Dict{String,String}()
    else
        _excluded_labels(configuration.excluded_mass_yields, mass_yields, records)
    end
    averaged = filter(mass_yields) do distribution
        @info "yield coverage against the primary distribution" yield = distribution.label coverage = coverage[distribution.label] reference =
            reference.label
        if haskey(excluded, distribution.label)
            @info "no total average over this yield distribution: excluded by configuration" yield =
                distribution.label reason = excluded[distribution.label]
            return false
        end
        coverage[distribution.label] >= configuration.min_yield_coverage && return true
        @warn "no total average over this yield distribution" yield = distribution.label coverage = coverage[distribution.label] floor =
            configuration.min_yield_coverage
        return false
    end
    return (; mass_yields, records, coverage, averaged)
end

# The trend refitted with one and two segments more than selected: how far ⟨R_T⟩ rests on the
# selected order, for each yield distribution averaged over.
function _segment_count_sensitivity(
    trend::ExtractedCurve,
    fitted,
    relation,
    configuration::Configuration,
    averaged::Vector{MassYield},
)
    sensitivity = Dict{String,Vector{Tuple{Int,Float64}}}()
    isempty(averaged) && return sensitivity
    selected = segment_count(trend.fit)
    settings = configuration.segments
    alternatives = ExtractedCurve[trend]
    for extra in 1:2
        order = selected + extra
        # The limit on the search holds for these refits as well. One above it is left out, as
        # an order the data cannot support is, rather than stopping a run whose own searches
        # stayed within it. The count of the single order bounds what the fit counts.
        sets = breakpoint_set_count(
            fitted.combined.A_H,
            order - 1;
            min_points_per_segment = settings.min_points_per_segment,
            required_windows = fitted.windows,
        )
        if sets > settings.max_breakpoint_sets
            @info "the trend is not refitted with $(order) segments: the search would examine \
                   $(sets) breakpoint sets, above max_breakpoint_sets" max_breakpoint_sets =
                settings.max_breakpoint_sets
            continue
        end
        refit = _segment(
            fitted.combined,
            relation,
            settings,
            SYSTEMATIC_TREND_LABEL,
            fitted.windows,
            configuration.system.compound.A;
            measured = fitted.pooled.measured,
            order = order,
        )
        refit isa ExtractedCurve && push!(alternatives, refit)
    end
    for distribution in averaged
        _overlaps(trend.R_T, distribution) || continue
        sensitivity[distribution.label] = [
            (segment_count(c.fit), first(total_average(c, distribution))) for c in alternatives
        ]
        @info "trend ⟨R_T⟩ against the number of segments" yield = distribution.label values = sensitivity[distribution.label]
    end
    return sensitivity
end

# The averaging a configuration names. The charge-resolved inversion is weighted by the mean total
# excitation of every fragmentation where a ⟨TKE⟩(A) dataset is given, and by p(Z, A_H) alone
# otherwise; the run says which.
function _averaging(
    label::AbstractString,
    masses::MassExcessTable,
    domain::FragmentationDomain,
    mean_kinetic_energy::Union{MeanKineticEnergy,Nothing},
)
    label == "charge_resolved" || return ratio_averaging(label)
    if mean_kinetic_energy === nothing
        @info "no ⟨TKE⟩(A) dataset: the charge-resolved inversion weights each fragmentation by \
               p(Z, A_H) alone, not by its excitation"
        return ChargeResolved()
    end
    flagged = _flagged(_qualifiers(mean_kinetic_energy.record))
    isempty(flagged) ||
        @warn "the retrieval records a qualifier on the ⟨TKE⟩(A) dataset; it is used, not \
               corrected" dataset = mean_kinetic_energy.label qualifiers = flagged
    @info "the charge-resolved inversion is weighted by ⟨TXE⟩ = Q + E*_CN - ⟨TKE⟩(A_H)" dataset =
        mean_kinetic_energy.label interpolated = mean_kinetic_energy.interpolated extrapolated =
        mean_kinetic_energy.extrapolated
    return ChargeResolved(mean_total_excitation(masses, domain, mean_kinetic_energy.values))
end

# Whether a yield distribution carries positive weight at any mass number of the curve.
function _overlaps(curve::RatioCurve, distribution::MassYield)
    total = 0.0
    for mass in curve.A_H
        entry = mass_yield(distribution, mass)
        entry === nothing || (total += entry[1])
    end
    return total > 0
end

# The total average of every segmented curve over every yield distribution. A pair that shares no
# mass number is omitted rather than reported as zero: a distribution covering only the light wing
# says nothing about a ratio defined on the heavy one.
function _total_averages(
    curves::Vector{ExtractedCurve},
    mass_yields::Vector{MassYield},
    heavy_masses::AbstractUnitRange{<:Integer},
)
    averages = Dict{String,Dict{String,TotalAverage}}()
    isempty(mass_yields) && return averages
    for curve in curves
        per_distribution = Dict{String,TotalAverage}()
        for distribution in mass_yields
            if _overlaps(curve.R_T, distribution)
                per_distribution[distribution.label] = TotalAverage(
                    curve, distribution, heavy_masses
                )
            else
                @warn "no total average" dataset = curve.label yield = distribution.label reason = "no mass number in common carries a positive yield"
            end
        end
        isempty(per_distribution) || (averages[curve.label] = per_distribution)
    end
    return averages
end

# Fit one ratio curve and carry it through to the temperature ratio, or the reason it cannot be
# fitted at all; `quiet` keeps that reason from the log where the caller reports it otherwise.
#
# `fit_segments` pins at the first abscissa it is given. r_ν = 1/2 is an identity at A₀/2 only, so
# the pin is applied to a curve whose first complete pair is the symmetric split and to no other;
# a dataset starting above it is fitted unpinned.
function _segment(
    curve::RatioCurve,
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    settings::SegmentSettings,
    label::AbstractString,
    windows::Vector{UnitRange{Int}},
    A₀::Integer;
    measured::AbstractVector{<:Real} = ones(length(curve)),
    order::Union{Integer,Nothing} = nothing,
    correlation::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    quiet::Bool = false,
)
    pinned = settings.pin_symmetric_split && !isempty(curve) && 2 * first(curve.A_H) == A₀
    fit = try
        fit_segments(
            curve.A_H,
            curve.ratio,
            curve.σ;
            max_segments = order === nothing ? settings.max_segments : order,
            min_segments = order === nothing ? 1 : order,
            measured = measured,
            min_points_per_segment = settings.min_points_per_segment,
            pinned_value = pinned ? 0.5 : nothing,
            required_windows = windows,
            bounds = (0.0, 1.0),
            correlation = correlation,
            max_breakpoint_sets = settings.max_breakpoint_sets,
        )
    catch exception
        exception isa InsufficientDataError || rethrow()
        order === nothing &&
            !quiet &&
            @warn "no segmented curve for this dataset" dataset = label reason =
                exception.msg
        return NoFit(exception.msg)
    end
    return ExtractedCurve(label, fit, curve, relation...)
end

# Checks of the identities that hold by construction at the symmetric split.
function _symmetry_diagnostics(
    configuration::Configuration,
    relation::Tuple{RatioAveraging,LevelDensityModel,FragmentationDomain},
    curves::Vector{ExtractedCurve},
)
    averaging, model, domain = relation
    warnings = String[]
    A₀ = configuration.system.compound.A
    invariant = something(symmetric_charge_set_is_invariant(domain), missing)
    if invariant === false
        push!(
            warnings,
            "the charge numbers retained at the symmetric split are not invariant under \
             Z -> Z₀ - Z, so r_ν = 1/2 does not return R_T = 1 there exactly; this occurs when \
             the charge distribution is polarized at A₀/2, or when Z₀ is odd",
        )
    end

    R_T_symmetric = missing
    if has_symmetric_split(configuration) && (A₀ ÷ 2) in domain.heavy_masses
        value = temperature_ratio(averaging, model, domain, A₀ ÷ 2, 0.5)
        if value !== nothing
            R_T_symmetric = value
            deviation = abs(value - 1)
            deviation < 1e-8 || push!(
                warnings,
                "r_ν = 1/2 inverts to R_T = $(round(value; sigdigits = 8)) at the symmetric \
                 split under ratio_averaging = \"$(ratio_averaging_label(averaging))\", \
                 $(round(deviation; sigdigits = 3)) from unity where the two fragments are the \
                 same nuclide",
            )
        end
    end

    pinned = String[]
    for c in curves
        c.fit.pinned_value === nothing && continue
        push!(pinned, c.label)
        # A pin anywhere but the symmetric split asserts r_ν = 1/2 where no identity holds.
        2 * c.fit.x₀ == A₀ || push!(
            warnings,
            "the curve \"$(c.label)\" is pinned at A_H = $(c.fit.x₀), which is not the \
             symmetric split of A₀ = $(A₀)",
        )
    end
    return SymmetryDiagnostics(invariant, R_T_symmetric, pinned, warnings)
end
