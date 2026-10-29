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

Controls of the segmented description: the largest number of segments examined, the smallest
number of data points a segment may hold, whether the ratio is pinned to one half at the
symmetric split, mass-number windows that must each contain a breakpoint, and the pair coverage a
dataset needs to be offered as a curve at all.

`required_windows` places a breakpoint where physics says there is one — the minimum at the heavy
magic fragment, `A_H` near 130, where the `Z = 50`, `N = 82` shell closure fixes the sharing. It
constrains the systematic-trend curve by default and the dataset curves only when
`windows_apply_to_datasets` is set, since a dataset that resolves the feature on its own should be
left to do so.

`min_pair_coverage` is the smallest fraction of the mass numbers of the fragmentation range at
which a dataset must provide a complete pair. Below it the dataset is read, diagnosed and pooled,
but no segmented curve is fitted to it alone: with pairs at few mass numbers the breakpoint search
cannot place the minimum where the data do not reach, and the curve it returns asserts structure
between the measurements that a consuming code could not tell from a measured feature.

`autocorrelation_lags` is the number of lags of the deviation correlogram the autocorrelation of
the systematic trend is fitted to, [`autocorrelation_decay`](@ref).
"""
struct SegmentSettings
    max_segments::Int
    min_points_per_segment::Int
    pin_symmetric_split::Bool
    required_windows::Vector{UnitRange{Int}}
    windows_apply_to_datasets::Bool
    min_pair_coverage::Float64
    autocorrelation_lags::Int
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

A complete, validated pipeline configuration, as [`load_configuration`](@ref) reads it from a TOML
file. The constructor enforces the rules that tie one setting to another — the heavy-mass range
to the fissioning nucleus, the pin to the symmetric split, the required windows to the range, the
zero polarization to a charge model that takes it, the excitation weights to the charge-resolved
inversion, a yield directory to its primary distribution and an exclusion of yields to a
directory — so that a configuration assembled by hand cannot hold a combination the pipeline
cannot honour either. What is checked against a single key, its type, bounds and choices, is the
loader's.

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

    function Configuration(
        system::FissioningSystem,
        fragmentation::FragmentationSettings,
        level_density::LevelDensitySettings,
        multiplicity_directory::AbstractString,
        excluded_datasets::AbstractDict{String,String},
        yield_directory::Union{AbstractString,Nothing},
        yield_file::Union{AbstractString,Nothing},
        excluded_mass_yields::AbstractDict{String,String},
        symmetrize_yields::Bool,
        min_yield_coverage::Real,
        segments::SegmentSettings,
        output::OutputSettings,
        min_retrieval_version::VersionNumber,
        source::AbstractString,
        data_directory::AbstractString,
    )
        _check_rules(
            system,
            fragmentation,
            level_density,
            yield_directory,
            yield_file,
            excluded_mass_yields,
            segments,
        )
        return new(
            system,
            fragmentation,
            level_density,
            multiplicity_directory,
            excluded_datasets,
            yield_directory,
            yield_file,
            excluded_mass_yields,
            symmetrize_yields,
            min_yield_coverage,
            segments,
            output,
            min_retrieval_version,
            source,
            data_directory,
        )
    end
end

# The rules that tie one setting to another, each named by the configuration keys it binds.
function _check_rules(
    system::FissioningSystem,
    fragmentation::FragmentationSettings,
    level_density::LevelDensitySettings,
    yield_directory::Union{AbstractString,Nothing},
    yield_file::Union{AbstractString,Nothing},
    excluded_mass_yields::AbstractDict,
    segments::SegmentSettings,
)
    A₀ = system.compound.A
    # The symmetric split is the smallest heavy-fragment mass there is: below it the heavy
    # fragment would be the lighter of the pair.
    symmetric_split = cld(A₀, 2)
    lower, upper = fragmentation.heavy_mass_min, fragmentation.heavy_mass_max
    _in_interval(lower, (symmetric_split, A₀ - 1), "fragmentation.heavy_mass_min")
    _in_interval(upper, (lower + 1, A₀ - 1), "fragmentation.heavy_mass_max")
    # Wahl's model fixes ΔZ at the symmetric split itself; the key would claim a setting the run
    # does not apply.
    fragmentation.zero_polarization_at_symmetry &&
        fragmentation.charge_distribution == "wahl" &&
        throw(
            ArgumentError(
                "fragmentation.zero_polarization_at_symmetry applies to a tabulated or a mean \
                 charge distribution; Wahl's model sets the polarization at A0/2 itself"
            ),
        )
    level_density.mean_kinetic_energy_file === nothing ||
        level_density.ratio_averaging == "charge_resolved" ||
        throw(
            ArgumentError(
                "level_density.mean_kinetic_energy_file weights the charge-resolved inversion \
                 only; level_density.ratio_averaging is $(repr(level_density.ratio_averaging))"
            ),
        )
    if segments.pin_symmetric_split
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
        lower == symmetric_split || throw(
            ArgumentError(
                "segments.pin_symmetric_split requires fragmentation.heavy_mass_min to be the \
                 symmetric split $(symmetric_split), got $(lower)",
            ),
        )
    end
    for (index, window) in enumerate(segments.required_windows)
        issubset(window, lower:upper) || throw(
            ArgumentError("segments.required_windows[$(index)] = $(window) lies outside the \
                           fragmentation range $(lower):$(upper)"),
        )
    end
    # A directory of distributions is measured against the primary one, and an exclusion of
    # distributions has a directory to act on or nothing.
    yield_directory === nothing ||
        yield_file !== nothing ||
        throw(
            ArgumentError(
                "yield.subdirectory needs yield.mass_yield_file, the primary Y(A) every \
                 distribution's coverage is measured against"
            ),
        )
    isempty(excluded_mass_yields) ||
        yield_directory !== nothing ||
        throw(
            ArgumentError(
                "yield.exclude needs yield.subdirectory: without a directory the primary \
                 distribution alone is averaged over, and an exclusion has nothing to act on"
            ),
        )
    return nothing
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
multiplicity ratio equals one half there by the identity of the two fragments, which is what the
pin relies on.
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

# The lowest ExforFissionData version whose retrieval records a run accepts where the
# configuration does not say: the first whose records name the datasets of one experiment
# (`correlated_with`, `correlation_relation`). The pooling of an experiment as one rests on them;
# inputs retrieved before would pool its datasets side by side, under the same run identifier.
const DEFAULT_MIN_RETRIEVAL_VERSION = "0.2.7"

struct Required end
const REQUIRED = Required()

# One configuration key: the type its value must have; its default, `REQUIRED` where it has none
# and `nothing` where the key is optional; and the closed interval or the choices it is held to,
# or `nothing`.
struct Key
    type::Type
    default::Any
    allowed::Union{Nothing,Tuple{Real,Real},Tuple{Vararg{String}}}
end

# Every key of every section. The loader checks each value against its entry here and nothing
# else; what depends on another key is checked where the settings are built, and by the
# `Configuration` constructor.
const CONFIGURATION_KEYS = (
    system = (
        target_A = Key(Integer, REQUIRED, (2, 400)),
        target_Z = Key(Integer, REQUIRED, (1, 118)),
        channel = Key(String, REQUIRED, CHANNELS),
        incident_energy = Key(Real, 0.0, (0, Inf)),
    ),
    fragmentation = (
        charges_per_mass = Key(Integer, REQUIRED, (1, 15)),
        # The symmetric split where absent; its bounds follow from the fissioning nucleus.
        heavy_mass_min = Key(Integer, nothing, nothing),
        heavy_mass_max = Key(Integer, REQUIRED, nothing),
        charge_distribution_file = Key(String, "wahl", nothing),
        zero_polarization_at_symmetry = Key(Bool, false, nothing),
    ),
    level_density = (
        model = Key(String, REQUIRED, MODELS),
        mass_excess_file = Key(String, SHIPPED_MASS_TABLE, nothing),
        shell_correction_file = Key(String, SHIPPED_SHELL_CORRECTIONS, nothing),
        deformed_branch = Key(Bool, true, nothing),
        ratio_averaging = Key(String, "charge_resolved", RATIO_AVERAGINGS),
        mean_kinetic_energy_file = Key(String, nothing, nothing),
    ),
    multiplicity = (
        subdirectory = Key(String, REQUIRED, nothing),
        exclude = Key(AbstractVector, nothing, nothing),
    ),
    yield = (
        subdirectory = Key(String, nothing, nothing),
        mass_yield_file = Key(String, REQUIRED, nothing),
        exclude = Key(AbstractVector, nothing, nothing),
        symmetrize = Key(Bool, true, nothing),
        min_yield_coverage = Key(Real, 0.3, (0, 1)),
    ),
    segments = (
        max_segments = Key(Integer, 6, (1, 12)),
        min_points_per_segment = Key(Integer, 4, (2, 50)),
        pin_symmetric_split = Key(Bool, true, nothing),
        required_windows = Key(AbstractVector, nothing, nothing),
        windows_apply_to_datasets = Key(Bool, false, nothing),
        min_pair_coverage = Key(Real, 0.3, (0, 1)),
        autocorrelation_lags = Key(
            Integer, DEFAULT_AUTOCORRELATION_LAGS, (1, CORRELOGRAM_LAGS)
        ),
    ),
    retrieval = (min_package_version = Key(String, DEFAULT_MIN_RETRIEVAL_VERSION, nothing),),
    output = (significant_digits = Key(Integer, 6, (1, 15)),),
)

# Sections a configuration may leave out: without [yield] no total average is taken, and without
# [retrieval] the default version floor holds.
const OPTIONAL_SECTIONS = (:yield, :retrieval)

# A key the loader does not read is refused rather than ignored: a misspelt or retired key would
# otherwise leave a run silently on the default, and the file would claim a setting the run did
# not honour.
function _refuse_unknown(table::AbstractDict, allowed, path::String)
    unknown = sort!(collect(setdiff(keys(table), allowed)))
    isempty(unknown) ||
        throw(ArgumentError("$(path) has no key $(join(map(repr, unknown), ", "))"))
    return nothing
end

# TOML reads `true` as a Bool, which Julia counts among the integers and the reals; a numeric key
# takes a number, not a truth value.
_conforms(value, ::Type{T}) where {T} = value isa T && (T === Bool || !(value isa Bool))

function _in_interval(value::Real, (lower, upper)::Tuple{Real,Real}, path::AbstractString)
    isfinite(value) && lower ≤ value ≤ upper || throw(
        ArgumentError("$(path) must be a finite number in [$(lower), $(upper)], got $(value)"),
    )
    return value
end

_check_allowed(value, ::Nothing, ::AbstractString) = value
function _check_allowed(value::Real, interval::Tuple{Real,Real}, path::AbstractString)
    return _in_interval(value, interval, path)
end
function _check_allowed(
    value::AbstractString, choices::Tuple{Vararg{String}}, path::AbstractString
)
    value in choices || throw(
        ArgumentError("$(path) must be one of $(join(choices, " | ")), got $(repr(value))")
    )
    return value
end

function _key_value(section::AbstractDict, key::Symbol, specification::Key, path::String)
    name = String(key)
    if !haskey(section, name)
        specification.default === REQUIRED &&
            throw(ArgumentError("configuration is missing the required key $(path)"))
        return specification.default
    end
    value = section[name]
    _conforms(value, specification.type) || throw(
        ArgumentError("$(path) must be of type $(specification.type), got $(typeof(value)) \
             ($(repr(value)))"),
    )
    return _check_allowed(value, specification.allowed, path)
end

# The values of one section, every key at its value or its default, as a named tuple; `nothing`
# for an optional section the document does not hold.
function _section(document::AbstractDict, name::Symbol, source::String)
    if !haskey(document, String(name))
        name in OPTIONAL_SECTIONS && return nothing
        throw(ArgumentError("configuration $(source) is missing the [$(name)] section"))
    end
    section = document[String(name)]
    section isa AbstractDict ||
        throw(ArgumentError("configuration $(source): [$(name)] must be a table"))
    specification = CONFIGURATION_KEYS[name]
    _refuse_unknown(section, String.(keys(specification)), "[$(name)]")
    return NamedTuple{keys(specification)}(
        Tuple(
            _key_value(section, key, specification[key], "$(name).$(key)") for
            key in keys(specification)
        ),
    )
end

# Datasets kept out of the pooling, each with its reason; a reason is required. An exclusion names
# its dataset by accession, which does not change when a second dataset of the same author and
# year is retrieved; the label form is for a dataset that carries no accession.
function _exclusions(raw::Union{AbstractVector,Nothing}, path::String)
    entries = @NamedTuple{by::Symbol, name::String, reason::String}[]
    raw === nothing && return entries
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
# or by its label where it carries none. A label that names a dataset with an accession is
# refused, with the accession to write in its place.
function _resolve_exclusions(entries, directory::AbstractString, path::String)
    files = _data_files(directory)
    records = retrieval_records(directory)
    _refuse_shared_accession(directory, files, records)
    labels = _unique_labels(files)
    accessions = [_accession(get(records, file, nothing), file) for file in files]
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
            isempty(accession) || throw(
                ArgumentError(
                    "$(path)[$(index)] names the dataset $(repr(entry.name)) by its label; \
                     a dataset that has an accession is excluded by it: write \
                     accession = \"$(accession)\"",
                ),
            )
            entry.name
        end
        haskey(exclusions, key) &&
            throw(ArgumentError("$(path) names $(repr(key)) more than once"))
        exclusions[key] = entry.reason
    end
    return exclusions
end

function _windows(raw::Union{AbstractVector,Nothing}, path::String)
    windows = UnitRange{Int}[]
    raw === nothing && return windows
    for (index, entry) in enumerate(raw)
        entry isa AbstractVector &&
        length(entry) == 2 &&
        all(e -> _conforms(e, Integer), entry) ||
            throw(ArgumentError("$(path)[$(index)] must be a two-element array of integers \
                                 [first, last], got $(repr(entry))"))
        first(entry) ≤ last(entry) ||
            throw(ArgumentError("$(path)[$(index)] must be ordered, got $(repr(entry))"))
        push!(windows, Int(first(entry)):Int(last(entry)))
    end
    return windows
end

# A file or a directory under the data directory, by the spelling the configuration gives.
function _input_file(root::AbstractString, relative::AbstractString, path::String)
    isempty(relative) && throw(ArgumentError("$(path) must name a file; omit it for none"))
    resolved = joinpath(root, relative)
    isfile(resolved) || throw(ArgumentError("$(path) does not exist: $(resolved)"))
    return resolved
end

function _input_directory(root::AbstractString, relative::AbstractString, path::String)
    resolved = joinpath(root, relative)
    isdir(resolved) || throw(ArgumentError("$(path) does not exist: $(resolved)"))
    return resolved
end

# A named source, or a file under the data directory.
function _source(root::AbstractString, value::AbstractString, named, path::String)
    value in named && return String(value)
    resolved = joinpath(root, value)
    isfile(resolved) || throw(
        ArgumentError(
            "$(path) must be one of $(join(named, " | ")) or a table under the data \
             directory; $(resolved) does not exist"
        ),
    )
    return resolved
end

"""
    load_configuration(path; data_directory) -> Configuration

Read and validate a pipeline configuration.

Every key is checked for presence, type, enumerated choice and numerical range, every input
file named by the configuration is checked for existence, and a section or key the loader does
not know is refused, before the pipeline is allowed to start; the rules that tie one key to
another are the [`Configuration`](@ref) constructor's. Failures throw an `ArgumentError` naming
the offending key, so that a configuration the pipeline cannot honour never begins a run. A
`true` or `false` where a number stands is refused, not read as one or zero. Every retrieval
record beside the inputs must state a parser version at or above
`retrieval.min_package_version`.

`data_directory` is the root the `subdirectory` and file keys are resolved against: `data/` of
the repository for the shipped configurations, a fixture directory in the tests. It is required,
so that what a configuration reads never depends on which project happens to be active.

# Example

```julia
configuration = load_configuration("config/U233_nth.toml"; data_directory = "data")
```
"""
function load_configuration(path::AbstractString; data_directory::AbstractString)
    isfile(path) || throw(ArgumentError("configuration file not found: $(path)"))
    document = try
        TOML.parsefile(path)
    catch exception
        throw(ArgumentError("configuration $(path) is not valid TOML: $(exception)"))
    end
    source = String(path)
    _refuse_unknown(document, String.(keys(CONFIGURATION_KEYS)), "configuration $(source)")
    sections = map(name -> _section(document, name, source), keys(CONFIGURATION_KEYS))
    sections = NamedTuple{keys(CONFIGURATION_KEYS)}(sections)
    root = data_directory

    entry = sections.system
    entry.target_Z < entry.target_A ||
        throw(ArgumentError("system.target_Z must be smaller than system.target_A, \
             got target_Z = $(entry.target_Z), target_A = $(entry.target_A)"))
    # Spontaneous fission has no incident particle, so an energy for it is a contradiction rather
    # than a harmless extra.
    entry.channel == "sf" &&
        entry.incident_energy != 0 &&
        throw(ArgumentError("system.incident_energy must be zero for spontaneous fission, \
             got $(entry.incident_energy) MeV with channel \"sf\""))
    target = Nuclide(Int(entry.target_Z), Int(entry.target_A))
    system = if entry.channel == "sf"
        spontaneous_fission(target)
    else
        neutron_induced_fission(target, Float64(entry.incident_energy), entry.channel)
    end

    entry = sections.fragmentation
    isodd(entry.charges_per_mass) || throw(
        ArgumentError(
            "fragmentation.charges_per_mass must be odd, got $(entry.charges_per_mass)"
        ),
    )
    fragmentation = FragmentationSettings(
        entry.charges_per_mass,
        something(entry.heavy_mass_min, cld(system.compound.A, 2)),
        entry.heavy_mass_max,
        _source(
            root,
            entry.charge_distribution_file,
            CHARGE_MODELS,
            "fragmentation.charge_distribution_file",
        ),
        entry.zero_polarization_at_symmetry,
    )

    entry = sections.level_density
    level_density = LevelDensitySettings(
        entry.model,
        if entry.mass_excess_file == SHIPPED_MASS_TABLE
            SHIPPED_MASS_TABLE
        else
            _input_file(root, entry.mass_excess_file, "level_density.mass_excess_file")
        end,
        if entry.shell_correction_file == SHIPPED_SHELL_CORRECTIONS
            SHIPPED_SHELL_CORRECTIONS
        else
            _input_file(
                root, entry.shell_correction_file, "level_density.shell_correction_file"
            )
        end,
        entry.deformed_branch,
        entry.ratio_averaging,
        if entry.mean_kinetic_energy_file === nothing
            nothing
        else
            _input_file(
                root, entry.mean_kinetic_energy_file, "level_density.mean_kinetic_energy_file"
            )
        end,
    )

    multiplicity_directory = _input_directory(
        root, sections.multiplicity.subdirectory, "multiplicity.subdirectory"
    )
    # Optional. Without it the run reports the mean over the fragment mass range only; with it,
    # the total average over each yield distribution, which is the quantity the literature quotes.
    # The primary distribution is always named: averaged over alone, or the reference the
    # coverage of every distribution of a directory is measured against.
    yields = something(sections.yield, _defaults(:yield))
    yield_file = if sections.yield === nothing
        nothing
    else
        _input_file(root, yields.mass_yield_file, "yield.mass_yield_file")
    end
    yield_directory = if yields.subdirectory === nothing
        nothing
    else
        _input_directory(root, yields.subdirectory, "yield.subdirectory")
    end

    floor = something(sections.retrieval, _defaults(:retrieval)).min_package_version
    min_retrieval_version = tryparse(VersionNumber, floor)
    min_retrieval_version === nothing && throw(
        ArgumentError(
            "retrieval.min_package_version must be a version number such as \"0.2.7\", got \
             $(repr(floor))"
        ),
    )
    # Before any entry of a record is read: a directory written by a release below the floor is
    # refused for that, and not for an entry an older release wrote in another form.
    energies = level_density.mean_kinetic_energy_file
    check_retrieval_versions(
        String[
            directory for directory in (
                multiplicity_directory,
                yield_directory,
                yield_file === nothing ? nothing : dirname(yield_file),
                energies === nothing ? nothing : dirname(energies),
            ) if directory !== nothing
        ],
        min_retrieval_version,
    )

    excluded_datasets = _resolve_exclusions(
        _exclusions(sections.multiplicity.exclude, "multiplicity.exclude"),
        multiplicity_directory,
        "multiplicity.exclude",
    )
    yield_exclusions = _exclusions(yields.exclude, "yield.exclude")
    excluded_mass_yields = Dict{String,String}()
    if yield_directory !== nothing
        excluded_mass_yields = _resolve_exclusions(
            yield_exclusions, yield_directory, "yield.exclude"
        )
    elseif !isempty(yield_exclusions)
        # Read for its form first, so that a malformed entry reports its own error.
        throw(
            ArgumentError(
                "yield.exclude needs yield.subdirectory: without a directory the primary \
                 distribution alone is averaged over, and an exclusion has nothing to act on"
            ),
        )
    end

    entry = sections.segments
    segments = SegmentSettings(
        entry.max_segments,
        entry.min_points_per_segment,
        entry.pin_symmetric_split,
        _windows(entry.required_windows, "segments.required_windows"),
        entry.windows_apply_to_datasets,
        Float64(entry.min_pair_coverage),
        entry.autocorrelation_lags,
    )

    configuration = Configuration(
        system,
        fragmentation,
        level_density,
        multiplicity_directory,
        excluded_datasets,
        yield_directory,
        yield_file,
        excluded_mass_yields,
        yields.symmetrize,
        Float64(yields.min_yield_coverage),
        segments,
        OutputSettings(sections.output.significant_digits),
        min_retrieval_version,
        source,
        String(root),
    )

    # What rests on the file itself. The ⟨TKE⟩(A) dataset must be one a retrieval produced, so that
    # its provenance and its qualifiers are on record, and it must reach the heavy-mass range.
    if energies !== nothing
        record = retrieval_record(energies)
        record === nothing && throw(
            ArgumentError(
                "level_density.mean_kinetic_energy_file: no retrieval record beside \
                 $(energies) lists it"
            ),
        )
        read_mean_kinetic_energy(energies, system, A_H_range(configuration); record = record)
    end
    return configuration
end

# The values an optional section takes where the document leaves it out.
_defaults(name::Symbol) = map(key -> key.default, CONFIGURATION_KEYS[name])
