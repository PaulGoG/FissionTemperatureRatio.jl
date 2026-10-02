# Pipeline configuration: parsing and validation of the TOML input.

"""
    FragmentationSettings

The fragmentation domain: how many charge numbers are taken per mass number, how far the
heavy-fragment mass range extends, and where the isobaric charge distribution comes from.

`charge_distribution` is `"wahl"`, the `Zₚ` model FissionFragmentsDomain resolves for the system
(Wahl's 1988 parameters for the four reactions evaluated there, the 2002 systematics otherwise,
the conventional means beyond them); `"mean"`, `ΔZ = -0.5` and `σ_Z = 0.6` at every mass; or the
path of a tabulated `A dZ sigma_Z`. `zero_polarization_at_symmetry` takes `ΔZ = 0` at `A₀/2` for
the last two, as the published extraction did; Wahl's model sets it there itself.
"""
struct FragmentationSettings
    charges_per_mass::Int
    heavy_mass_min::Int
    heavy_mass_max::Int
    charge_distribution::String
    zero_polarization_at_symmetry::Bool
end

"""
    LevelDensitySettings

The level density model, the tables it is evaluated from, and how the parameter ratio of
complementary fragments enters the inversion.

- `model`: `"BSFG"` or `"GC"`.
- `mass_excess_file`: `"ame2020"`, the evaluation FissionFragmentsDomain ships, or the path of
  another table in its layout.
- `shell_correction_file`: `"gilbert_cameron_1965"`, Table III of Gilbert and Cameron as
  FissionFragmentsDomain ships it, or the path of another table in its layout; read by `"GC"`
  only.
- `deformed_branch`: whether Gilbert-Cameron applies its deformed-nucleus formula where it
  prescribes it; inert for `"BSFG"`.
- `ratio_averaging`: `"charge_resolved"`, `"ratio_of_means"` or `"mean_of_ratios"`.
- `mean_kinetic_energy_file`: the `⟨TKE⟩(A)` dataset that weights the charge-resolved inversion
  by excitation, or `nothing` for weights `p(Z, A_H)` alone.

The model is named rather than constructed here, because constructing it means reading its data;
[`build_level_density_model`](@ref) does that when a run starts.
"""
struct LevelDensitySettings
    model::String
    mass_excess_file::String
    shell_correction_file::String
    deformed_branch::Bool
    ratio_averaging::String
    mean_kinetic_energy_file::Union{String,Nothing}
end

"""
    SHIPPED_MASS_TABLE

The configuration value naming the atomic mass evaluation FissionFragmentsDomain ships,
`"ame2020"`.
"""
const SHIPPED_MASS_TABLE = "ame2020"

"""
    SHIPPED_SHELL_CORRECTIONS

The configuration value naming the Gilbert-Cameron shell corrections FissionFragmentsDomain
ships, `"gilbert_cameron_1965"`: Table III of Can. J. Phys. 43, 1446 (1965).
"""
const SHIPPED_SHELL_CORRECTIONS = "gilbert_cameron_1965"

"""
    build_mass_table(settings) -> MassExcessTable

The mass excess table `settings` names: the shipped 2020 evaluation for `"ame2020"`, otherwise
the table at the given path.
"""
function build_mass_table(settings::LevelDensitySettings)
    path = if settings.mass_excess_file == SHIPPED_MASS_TABLE
        String(AME2020_MASS_EXCESS_FILE)
    else
        settings.mass_excess_file
    end
    return read_mass_excess_table(path)
end

"""
    build_level_density_model(settings, masses) -> LevelDensityModel

Construct the level density model `settings` names: the back-shifted Fermi gas on the mass table
`masses`, or Gilbert-Cameron on the shell corrections it names — the shipped Table III for
`"gilbert_cameron_1965"` — with or without its deformed branch.

Throws an `ArgumentError` if the model is unknown. [`load_configuration`](@ref) already rejects
it, so this guards settings assembled by hand.
"""
function build_level_density_model(settings::LevelDensitySettings, masses::MassExcessTable)
    if settings.model == "BSFG"
        return BackShiftedFermiGas(masses)
    elseif settings.model == "GC"
        path = if settings.shell_correction_file == SHIPPED_SHELL_CORRECTIONS
            String(GILBERT_CAMERON_SHELL_CORRECTION_FILE)
        else
            settings.shell_correction_file
        end
        return GilbertCameron(
            read_shell_correction_table(path); deformed_branch = settings.deformed_branch
        )
    end
    return throw(ArgumentError("unknown level density model: $(settings.model)"))
end

"""
    SegmentSettings

Controls of the piecewise-linear parameterization: the largest number of segments examined, the
smallest number of data points and the smallest extent in mass units a segment may have, whether
the ratio is pinned to one half at the symmetric split, mass-number windows that must each
contain a breakpoint, and the pair coverage a dataset needs to be offered as a curve at all.

`required_windows` places a breakpoint where physics says there is one — the minimum at the heavy
magic fragment, `A_H` near 130, where the `Z = 50`, `N = 82` shell closure fixes the sharing. It
constrains the systematic-trend curve by default and the per-dataset parameterizations only when
`windows_apply_to_datasets` is set, since a dataset that resolves the feature on its own should be
left to do so.

`min_pair_coverage` is the smallest fraction of the mass numbers of the fragmentation range at
which a dataset must provide a complete pair. Below it the dataset is read, diagnosed and pooled,
but no segmented curve is fitted to it alone: with pairs at few mass numbers the breakpoint search
cannot place the minimum where the data do not reach, and the curve it returns asserts structure
between the measurements that a consuming code could not tell from a measured feature.
"""
struct SegmentSettings
    max_segments::Int
    min_points_per_segment::Int
    min_segment_span::Int
    pin_symmetric_split::Bool
    required_windows::Vector{UnitRange{Int}}
    windows_apply_to_datasets::Bool
    min_pair_coverage::Float64
end

"""
    OutputSettings

Rounding of tabulated output.

`significant_digits` is significant figures, not decimal places: an uncertainty of `1.2e-5` and a
ratio of `1.1` are both written to the same number of meaningful digits, which a fixed number of
decimals cannot do for two quantities of such different magnitude.
"""
struct OutputSettings
    significant_digits::Int
end

"""
    Configuration

A complete, validated pipeline configuration. Construct with [`load_configuration`](@ref) rather
than directly, so that the constraints documented in the TOML file are enforced.

`data_directory` is the root every input path was resolved against; the run identifier and the
metadata spell input paths relative to it, so that a run made on another machine names the same
inputs.

`yield_file` is the primary mass yield distribution of the system. Alone, it is the one the total
average is taken over; with `yield_directory`, every distribution of the directory is averaged
over and `yield_file` is the reference their coverage is measured against
([`mass_yield_coverage`](@ref)). `excluded_mass_yields` names distributions of the directory kept
out of every average, each with its reason; they are still read and reported. Both
`excluded_datasets` and `excluded_mass_yields` are keyed by the EXFOR accession of each dataset,
or by its label where it carries none.

`symmetrize_yields` imposes the pre-neutron identity `Y(A) = Y(A₀ - A)` on every mass yield
distribution before the total average: where both complements are measured each takes their mean.
The published extraction averaged over the distributions as measured.

`min_yield_coverage` is the share of the primary distribution's heavy-fragment yield,
[`mass_yield_coverage`](@ref), a distribution must cover to be averaged over; below it the total
average would describe those masses, not the fission yield.

`min_retrieval_version` is the lowest ExforFissionData version whose retrieval records the run
accepts, [`check_retrieval_versions`](@ref).
"""
struct Configuration
    system::FissioningSystem
    fragmentation::FragmentationSettings
    level_density::LevelDensitySettings
    multiplicity_directory::String
    excluded_datasets::Dict{String,String}
    yield_directory::Union{String,Nothing}
    yield_file::Union{String,Nothing}
    excluded_mass_yields::Dict{String,String}
    symmetrize_yields::Bool
    min_yield_coverage::Float64
    segments::SegmentSettings
    output::OutputSettings
    min_retrieval_version::VersionNumber
    source::String
    data_directory::String
end

"""
    A_H_range(configuration) -> UnitRange{Int}

Heavy-fragment mass numbers spanned by the fragmentation range, `heavy_mass_min:heavy_mass_max`.
"""
function A_H_range(configuration::Configuration)
    return configuration.fragmentation.heavy_mass_min:configuration.fragmentation.heavy_mass_max
end

"""
    has_symmetric_split(configuration) -> Bool

Whether the fissioning nucleus has a symmetric split, that is, whether `A₀` is even. The
multiplicity ratio equals one half there by the identity of the two fragments, which is what
pinning the parameterization relies on.
"""
has_symmetric_split(configuration::Configuration) = iseven(configuration.system.compound.A)

"""
    build_charge_model(configuration, masses) -> ChargeModel

The isobaric charge distribution the configuration names: Wahl's `Zₚ` model as
FissionFragmentsDomain resolves it for the system, the conventional means, or a tabulated
distribution read with its reader.

A tabulated distribution is labelled by its path relative to the data directory, which is how the
run record names it; the label then does not depend on where the data directory sits.
"""
function build_charge_model(configuration::Configuration, masses::MassExcessTable)
    source = configuration.fragmentation.charge_distribution
    source == "wahl" && return charge_model(masses, configuration.system)
    source == "mean" && return mean_charge_distribution()
    table = read_charge_distribution(source)
    return ChargeDistribution(
        table.ΔZ,
        table.σ_Z,
        table.default_ΔZ,
        table.default_σ_Z,
        relpath(source, configuration.data_directory),
    )
end

const CHANNELS = ("sf", "nth", "nres", "nfast")
const MODELS = ("BSFG", "GC")
const CHARGE_MODELS = ("wahl", "mean")

function _section(document::AbstractDict, name::String, source::String)
    haskey(document, name) ||
        throw(ArgumentError("configuration $(source) is missing the [$(name)] section"))
    section = document[name]
    section isa AbstractDict ||
        throw(ArgumentError("configuration $(source): [$(name)] must be a table"))
    return section
end

# A key the loader does not read is refused rather than ignored: a misspelt or retired key would
# otherwise leave a run silently on the default, and the file would claim a setting the run did
# not honour.
function _refuse_unknown(table::AbstractDict, allowed, path::String)
    unknown = sort!(collect(setdiff(keys(table), allowed)))
    isempty(unknown) ||
        throw(ArgumentError("$(path) has no key $(join(map(repr, unknown), ", "))"))
    return nothing
end

const SECTIONS = (
    "system",
    "fragmentation",
    "level_density",
    "multiplicity",
    "yield",
    "segments",
    "retrieval",
    "output",
)
const SYSTEM_KEYS = ("target_A", "target_Z", "channel", "incident_energy")
const FRAGMENTATION_KEYS = (
    "charges_per_mass",
    "heavy_mass_min",
    "heavy_mass_max",
    "charge_distribution_file",
    "zero_polarization_at_symmetry",
)
const LEVEL_DENSITY_KEYS = (
    "model",
    "mass_excess_file",
    "shell_correction_file",
    "deformed_branch",
    "ratio_averaging",
    "mean_kinetic_energy_file",
)
const MULTIPLICITY_KEYS = ("subdirectory", "exclude")
const YIELD_KEYS = (
    "subdirectory", "mass_yield_file", "exclude", "symmetrize", "min_yield_coverage"
)
const SEGMENT_KEYS = (
    "max_segments",
    "min_points_per_segment",
    "min_segment_span",
    "pin_symmetric_split",
    "required_windows",
    "windows_apply_to_datasets",
    "min_pair_coverage",
)
const RETRIEVAL_KEYS = ("min_package_version",)
const OUTPUT_KEYS = ("significant_digits",)

function _value(section::AbstractDict, key::String, ::Type{T}, path::String) where {T}
    haskey(section, key) ||
        throw(ArgumentError("configuration is missing the required key $(path)"))
    value = section[key]
    value isa T || throw(
        ArgumentError("$(path) must be of type $(T), got $(typeof(value)) ($(repr(value)))")
    )
    return value
end

function _value(
    section::AbstractDict, key::String, ::Type{T}, path::String, default::T
) where {T}
    haskey(section, key) || return default
    return _value(section, key, T, path)
end

function _in_bounds(value::Real, path::String; min = -Inf, max = Inf, exclusive_min = false)
    lower_ok = exclusive_min ? value > min : value ≥ min
    lower_ok && value ≤ max || throw(
        ArgumentError("$(path) must lie in $(exclusive_min ? "(" : "[")$(min), $(max)], \
                       got $(value)")
    )
    return value
end

function _one_of(value::String, options, path::String)
    value in options || throw(
        ArgumentError("$(path) must be one of $(join(options, " | ")), got $(repr(value))")
    )
    return value
end

# Datasets kept out of the pooling, each with the reason written down. A reason is required:
# excluding a measurement is a judgement, and an unexplained one is indistinguishable from a
# mistake to anyone reading the configuration later. An exclusion names its dataset by accession,
# which does not change when a second dataset of the same author and year is retrieved, as a label
# does.
function _exclusions(section::AbstractDict, path::String)
    entries = @NamedTuple{by::Symbol, name::String, reason::String}[]
    haskey(section, "exclude") || return entries
    raw = section["exclude"]
    raw isa AbstractVector || throw(
        ArgumentError("$(path) must be an array of tables, each with `accession` and `reason`"),
    )
    for (index, entry) in enumerate(raw)
        well_formed =
            entry isa AbstractDict &&
            haskey(entry, "reason") &&
            entry["reason"] isa String &&
            count(key -> haskey(entry, key), ("accession", "dataset")) == 1 &&
            all(key -> !haskey(entry, key) || entry[key] isa String, ("accession", "dataset"))
        well_formed || throw(
            ArgumentError(
                "$(path)[$(index)] must be a table with a string `reason` and one of the \
                 string keys `accession` and `dataset`"
            ),
        )
        _refuse_unknown(entry, ("accession", "dataset", "reason"), "$(path)[$(index)]")
        isempty(strip(entry["reason"])) &&
            throw(ArgumentError("$(path)[$(index)] must give a non-empty reason"))
        by = haskey(entry, "accession") ? :accession : :dataset
        name = entry[String(by)]
        by === :accession &&
            !occursin(ACCESSION_PATTERN, name) &&
            throw(
                ArgumentError(
                    "$(path)[$(index)].accession must be an EXFOR dataset identifier, 8 or 9 \
                 characters from 0-9 and A-Z, got $(repr(name))"
                ),
            )
        push!(entries, (by = by, name = name, reason = entry["reason"]))
    end
    return entries
end

# The exclusions of one section against the data files of its directory. Each must name exactly
# one dataset held there, and is keyed as `_exclusion_key` keys that dataset: by its accession,
# or by its label where it carries none.
function _resolve_exclusions(entries, directory::AbstractString, path::String)
    files = _data_files(directory)
    _refuse_shared_accession(directory, files)
    labels = _unique_labels(files)
    accessions = [_dataset_accession(joinpath(directory, file)) for file in files]
    exclusions = Dict{String,String}()
    for (index, entry) in enumerate(entries)
        key = if entry.by === :accession
            entry.name in accessions || throw(
                ArgumentError(
                    "$(path)[$(index)] excludes the accession $(repr(entry.name)), which no \
                     dataset of $(directory) carries; it holds \
                     $(join(filter(!isempty, accessions), ", "))",
                ),
            )
            entry.name
        else
            position = findfirst(==(entry.name), labels)
            position === nothing && throw(
                ArgumentError(
                    "$(path)[$(index)] excludes the dataset labelled $(repr(entry.name)), \
                     which is not among those of $(directory): $(join(labels, ", "))",
                ),
            )
            accession = accessions[position]
            if isempty(accession)
                entry.name
            else
                @warn "$(path)[$(index)] names a dataset by its label, which changes when a \
                       second dataset of the same author and year is retrieved; this form is \
                       deprecated for a dataset that has an accession: write \
                       accession = \"$(accession)\"" dataset = entry.name
                accession
            end
        end
        haskey(exclusions, key) &&
            throw(ArgumentError("$(path) names $(repr(key)) more than once"))
        exclusions[key] = entry.reason
    end
    return exclusions
end

function _windows(section::AbstractDict, path::String)
    haskey(section, "required_windows") || return UnitRange{Int}[]
    raw = section["required_windows"]
    raw isa AbstractVector ||
        throw(ArgumentError("$(path) must be an array of two-element arrays [first, last]"))
    windows = UnitRange{Int}[]
    for (index, entry) in enumerate(raw)
        entry isa AbstractVector && length(entry) == 2 && all(e -> e isa Integer, entry) ||
            throw(ArgumentError("$(path)[$(index)] must be a two-element array of integers \
                           [first, last], got $(repr(entry))"))
        first(entry) ≤ last(entry) ||
            throw(ArgumentError("$(path)[$(index)] must be ordered, got $(repr(entry))"))
        push!(windows, Int(first(entry)):Int(last(entry)))
    end
    return windows
end

# A path under the data root. `subdirectory` names a folder under a root; `directory` names a root.
function _subdirectory(
    section::AbstractDict, path::String, root::AbstractString; required::Bool = true
)
    if !required && !haskey(section, "subdirectory")
        return nothing
    end
    resolved = joinpath(root, _value(section, "subdirectory", String, path))
    isdir(resolved) || throw(ArgumentError("$(path) does not exist: $(resolved)"))
    return resolved
end

function _input_file(
    section::AbstractDict,
    key::String,
    path::String,
    root::AbstractString;
    required::Bool = true,
)
    relative =
        required ? _value(section, key, String, path) : _value(section, key, String, path, "")
    if isempty(relative)
        required && throw(ArgumentError("$(path) must not be empty"))
        return nothing
    end
    resolved = joinpath(root, relative)
    isfile(resolved) || throw(ArgumentError("$(path) does not exist: $(resolved)"))
    return resolved
end

"""
    load_configuration(path; data_directory = datadir()) -> Configuration

Read and validate a pipeline configuration.

Every key is checked for presence, type, enumerated choice and numerical range, every input
file named by the configuration is checked for existence, and a section or key the loader does
not know is refused, before the pipeline is allowed to start. Failures throw an `ArgumentError`
naming the offending key, so that a configuration the pipeline cannot honour never begins a run.
Every retrieval record beside the inputs must state a parser version at or above
`retrieval.min_package_version`.

`data_directory` is the root the `subdirectory` and file keys are resolved against; it exists so
that tests can point at a fixture directory.

# Example

```julia
configuration = load_configuration(joinpath(projectdir(), "config", "U233_nth.toml"))
```
"""
function load_configuration(path::AbstractString; data_directory::AbstractString = datadir())
    isfile(path) || throw(ArgumentError("configuration file not found: $(path)"))
    document = try
        TOML.parsefile(path)
    catch exception
        throw(ArgumentError("configuration $(path) is not valid TOML: $(exception)"))
    end
    source = String(path)
    _refuse_unknown(document, SECTIONS, "configuration $(source)")

    system_section = _section(document, "system", source)
    _refuse_unknown(system_section, SYSTEM_KEYS, "[system]")
    target_A = Int(
        _in_bounds(
            _value(system_section, "target_A", Integer, "system.target_A"),
            "system.target_A";
            min = 2,
            max = 400,
        ),
    )
    target_Z = Int(
        _in_bounds(
            _value(system_section, "target_Z", Integer, "system.target_Z"),
            "system.target_Z";
            min = 1,
            max = 118,
        ),
    )
    target_Z < target_A ||
        throw(ArgumentError("system.target_Z must be smaller than system.target_A, \
             got target_Z = $(target_Z), target_A = $(target_A)"))
    channel = _one_of(
        _value(system_section, "channel", String, "system.channel"), CHANNELS, "system.channel"
    )
    incident_energy = Float64(
        _in_bounds(
            _value(system_section, "incident_energy", Real, "system.incident_energy", 0.0),
            "system.incident_energy";
            min = 0,
        ),
    )
    # Spontaneous fission has no incident particle, so an energy for it is a contradiction rather
    # than a harmless extra.
    channel == "sf" &&
        incident_energy != 0 &&
        throw(ArgumentError("system.incident_energy must be zero for spontaneous fission, \
             got $(incident_energy) MeV with channel \"sf\""))
    target = Nuclide(target_Z, target_A)
    system = if channel == "sf"
        spontaneous_fission(target)
    else
        neutron_induced_fission(target, incident_energy, channel)
    end
    A₀ = system.compound.A
    symmetric_split = cld(A₀, 2)

    fragmentation_section = _section(document, "fragmentation", source)
    _refuse_unknown(fragmentation_section, FRAGMENTATION_KEYS, "[fragmentation]")
    charges_per_mass = Int(
        _in_bounds(
            _value(
                fragmentation_section,
                "charges_per_mass",
                Integer,
                "fragmentation.charges_per_mass",
            ),
            "fragmentation.charges_per_mass";
            min = 1,
            max = 15,
        ),
    )
    isodd(charges_per_mass) || throw(
        ArgumentError("fragmentation.charges_per_mass must be odd, got $(charges_per_mass)")
    )
    # The symmetric split is the smallest heavy-fragment mass there is, and the default: below it
    # the heavy fragment would be the lighter of the pair.
    heavy_mass_min = Int(
        _in_bounds(
            _value(
                fragmentation_section,
                "heavy_mass_min",
                Integer,
                "fragmentation.heavy_mass_min",
                symmetric_split,
            ),
            "fragmentation.heavy_mass_min";
            min = symmetric_split,
            max = A₀ - 1,
        ),
    )
    heavy_mass_max = Int(
        _in_bounds(
            _value(
                fragmentation_section, "heavy_mass_max", Integer, "fragmentation.heavy_mass_max"
            ),
            "fragmentation.heavy_mass_max";
            min = heavy_mass_min + 1,
            max = A₀ - 1,
        ),
    )
    charge_distribution = _value(
        fragmentation_section,
        "charge_distribution_file",
        String,
        "fragmentation.charge_distribution_file",
        "wahl",
    )
    if !(charge_distribution in CHARGE_MODELS)
        resolved = joinpath(data_directory, charge_distribution)
        isfile(resolved) || throw(
            ArgumentError(
                "fragmentation.charge_distribution_file must be one of \
                 $(join(CHARGE_MODELS, " | ")) or a table under the data directory; \
                 $(resolved) does not exist"
            ),
        )
        charge_distribution = resolved
    end
    zero_polarization = _value(
        fragmentation_section,
        "zero_polarization_at_symmetry",
        Bool,
        "fragmentation.zero_polarization_at_symmetry",
        false,
    )
    # Wahl's model fixes ΔZ at the symmetric split itself; the key would claim a setting the run
    # does not apply.
    zero_polarization &&
        charge_distribution == "wahl" &&
        throw(
            ArgumentError(
                "fragmentation.zero_polarization_at_symmetry applies to a tabulated or a mean \
                 charge distribution; Wahl's model sets the polarization at A0/2 itself"
            ),
        )

    level_density_section = _section(document, "level_density", source)
    _refuse_unknown(level_density_section, LEVEL_DENSITY_KEYS, "[level_density]")
    model = _one_of(
        _value(level_density_section, "model", String, "level_density.model"),
        MODELS,
        "level_density.model",
    )
    mass_excess_file = _value(
        level_density_section,
        "mass_excess_file",
        String,
        "level_density.mass_excess_file",
        SHIPPED_MASS_TABLE,
    )
    if mass_excess_file != SHIPPED_MASS_TABLE
        mass_excess_file = _input_file(
            level_density_section,
            "mass_excess_file",
            "level_density.mass_excess_file",
            data_directory,
        )
    end
    shell_correction_file = _value(
        level_density_section,
        "shell_correction_file",
        String,
        "level_density.shell_correction_file",
        SHIPPED_SHELL_CORRECTIONS,
    )
    if shell_correction_file != SHIPPED_SHELL_CORRECTIONS
        shell_correction_file = _input_file(
            level_density_section,
            "shell_correction_file",
            "level_density.shell_correction_file",
            data_directory,
        )
    end
    deformed_branch = _value(
        level_density_section, "deformed_branch", Bool, "level_density.deformed_branch", true
    )
    averaging = _one_of(
        _value(
            level_density_section,
            "ratio_averaging",
            String,
            "level_density.ratio_averaging",
            "charge_resolved",
        ),
        RATIO_AVERAGINGS,
        "level_density.ratio_averaging",
    )
    mean_kinetic_energy_file = _input_file(
        level_density_section,
        "mean_kinetic_energy_file",
        "level_density.mean_kinetic_energy_file",
        data_directory;
        required = false,
    )
    if mean_kinetic_energy_file !== nothing
        averaging == "charge_resolved" || throw(
            ArgumentError(
                "level_density.mean_kinetic_energy_file weights the charge-resolved inversion \
                 only; level_density.ratio_averaging is $(repr(averaging))",
            ),
        )
        # The dataset must be one a retrieval produced, so that its provenance and its
        # qualifiers are on record, and it must reach the heavy-mass range.
        retrieval_record(mean_kinetic_energy_file) === nothing && throw(
            ArgumentError(
                "level_density.mean_kinetic_energy_file: no retrieval record beside \
                 $(mean_kinetic_energy_file) lists it"
            ),
        )
        read_mean_kinetic_energy(
            mean_kinetic_energy_file, system, heavy_mass_min:heavy_mass_max
        )
    end

    multiplicity_section = _section(document, "multiplicity", source)
    _refuse_unknown(multiplicity_section, MULTIPLICITY_KEYS, "[multiplicity]")
    multiplicity_directory = _subdirectory(
        multiplicity_section, "multiplicity.subdirectory", data_directory
    )
    excluded_datasets = _resolve_exclusions(
        _exclusions(multiplicity_section, "multiplicity.exclude"),
        multiplicity_directory,
        "multiplicity.exclude",
    )

    # Optional. Without it the run reports the mean over the fragment mass range only; with it,
    # the total average over each yield distribution, which is the quantity the literature quotes.
    yield_directory = nothing
    yield_file = nothing
    excluded_mass_yields = Dict{String,String}()
    if haskey(document, "yield")
        yield_section = _section(document, "yield", source)
        _refuse_unknown(yield_section, YIELD_KEYS, "[yield]")
        # The primary distribution is always named: averaged over alone, or the reference the
        # coverage of every distribution of a directory is measured against.
        haskey(yield_section, "mass_yield_file") || throw(
            ArgumentError(
                "[yield] needs yield.mass_yield_file, the primary Y(A) of the system: it is \
                 averaged over alone, or with yield.subdirectory is the reference every \
                 distribution's coverage is measured against",
            ),
        )
        yield_file = _input_file(
            yield_section, "mass_yield_file", "yield.mass_yield_file", data_directory
        )
        if haskey(yield_section, "subdirectory")
            yield_directory = _subdirectory(yield_section, "yield.subdirectory", data_directory)
            excluded_mass_yields = _resolve_exclusions(
                _exclusions(yield_section, "yield.exclude"), yield_directory, "yield.exclude"
            )
        else
            # Without a directory the primary distribution alone is averaged over, and an
            # exclusion has nothing to act on and nothing to be checked against. It is read for
            # its form and then dropped, so that neither the run identifier nor the metadata
            # records an exclusion that was not applied.
            isempty(_exclusions(yield_section, "yield.exclude")) ||
                @warn "yield.exclude has no effect without yield.subdirectory: the primary \
                       distribution alone is averaged over, and the exclusions are neither \
                       applied nor recorded"
        end
    end
    symmetrize_yields = if haskey(document, "yield")
        _value(document["yield"], "symmetrize", Bool, "yield.symmetrize", true)
    else
        true
    end
    min_yield_coverage = if haskey(document, "yield")
        Float64(
            _in_bounds(
                _value(
                    document["yield"],
                    "min_yield_coverage",
                    Real,
                    "yield.min_yield_coverage",
                    0.3,
                ),
                "yield.min_yield_coverage";
                min = 0,
                max = 1,
            ),
        )
    else
        0.3
    end

    segments_section = _section(document, "segments", source)
    haskey(segments_section, "min_dataset_coverage") && throw(
        ArgumentError(
            "segments.min_dataset_coverage is retired, having gated two measures: set \
             segments.min_pair_coverage, the fraction of the range's heavy masses a ν(A) \
             dataset must pair to offer a curve, and yield.min_yield_coverage, the share of \
             the primary distribution's heavy-fragment yield a Y(A) must cover to be \
             averaged over",
        ),
    )
    _refuse_unknown(segments_section, SEGMENT_KEYS, "[segments]")
    max_segments = Int(
        _in_bounds(
            _value(segments_section, "max_segments", Integer, "segments.max_segments", 6),
            "segments.max_segments";
            min = 1,
            max = 12,
        ),
    )
    min_points = Int(
        _in_bounds(
            _value(
                segments_section,
                "min_points_per_segment",
                Integer,
                "segments.min_points_per_segment",
                4,
            ),
            "segments.min_points_per_segment";
            min = 2,
            max = 50,
        ),
    )
    min_span = Int(
        _in_bounds(
            _value(
                segments_section, "min_segment_span", Integer, "segments.min_segment_span", 3
            ),
            "segments.min_segment_span";
            min = 1,
            max = 50,
        ),
    )
    pin = _value(
        segments_section, "pin_symmetric_split", Bool, "segments.pin_symmetric_split", true
    )
    if pin
        isodd(A₀) && throw(
            ArgumentError(
                "segments.pin_symmetric_split requires an even fissioning mass number, since the \
                 ratio equals one half only where the two fragments are identical; got \
                 A0 = $(A₀)",
            ),
        )
        # The fit is pinned at the first abscissa of the range, so that abscissa has to be the
        # symmetric split; pinning a range that starts above it would fix the ratio to one half
        # where nothing says it is.
        heavy_mass_min == symmetric_split || throw(
            ArgumentError(
                "segments.pin_symmetric_split requires fragmentation.heavy_mass_min to be the \
                 symmetric split $(symmetric_split), got $(heavy_mass_min)",
            ),
        )
    end
    windows = _windows(segments_section, "segments.required_windows")
    windows_apply = _value(
        segments_section,
        "windows_apply_to_datasets",
        Bool,
        "segments.windows_apply_to_datasets",
        false,
    )
    min_pair_coverage = Float64(
        _in_bounds(
            _value(
                segments_section, "min_pair_coverage", Real, "segments.min_pair_coverage", 0.3
            ),
            "segments.min_pair_coverage";
            min = 0,
            max = 1,
        ),
    )
    for (index, window) in enumerate(windows)
        issubset(window, heavy_mass_min:heavy_mass_max) || throw(
            ArgumentError("segments.required_windows[$(index)] = $(window) lies outside the \
                           fragmentation range $(heavy_mass_min):$(heavy_mass_max)"),
        )
    end

    output_section = _section(document, "output", source)
    _refuse_unknown(output_section, OUTPUT_KEYS, "[output]")
    significant_digits = Int(
        _in_bounds(
            _value(
                output_section, "significant_digits", Integer, "output.significant_digits", 6
            ),
            "output.significant_digits";
            min = 1,
            max = 15,
        ),
    )

    min_package_version = if haskey(document, "retrieval")
        retrieval_section = _section(document, "retrieval", source)
        _refuse_unknown(retrieval_section, RETRIEVAL_KEYS, "[retrieval]")
        _value(
            retrieval_section,
            "min_package_version",
            String,
            "retrieval.min_package_version",
            "0.2.3",
        )
    else
        "0.2.3"
    end
    min_retrieval_version = tryparse(VersionNumber, min_package_version)
    min_retrieval_version === nothing && throw(
        ArgumentError(
            "retrieval.min_package_version must be a version number such as \"0.2.3\", got \
             $(repr(min_package_version))"
        ),
    )
    check_retrieval_versions(
        String[
            directory for directory in (
                multiplicity_directory,
                yield_directory,
                yield_file === nothing ? nothing : dirname(yield_file),
                if mean_kinetic_energy_file === nothing
                    nothing
                else
                    dirname(mean_kinetic_energy_file)
                end,
            ) if directory !== nothing
        ],
        min_retrieval_version,
    )

    return Configuration(
        system,
        FragmentationSettings(
            charges_per_mass,
            heavy_mass_min,
            heavy_mass_max,
            charge_distribution,
            zero_polarization,
        ),
        LevelDensitySettings(
            model,
            mass_excess_file,
            shell_correction_file,
            deformed_branch,
            averaging,
            mean_kinetic_energy_file,
        ),
        multiplicity_directory,
        excluded_datasets,
        yield_directory,
        yield_file,
        excluded_mass_yields,
        symmetrize_yields,
        min_yield_coverage,
        SegmentSettings(
            max_segments, min_points, min_span, pin, windows, windows_apply, min_pair_coverage
        ),
        OutputSettings(significant_digits),
        min_retrieval_version,
        source,
        String(data_directory),
    )
end
