# Pipeline configuration: parsing and validation of the TOML input.

"""
    SystemSpecification

The fissioning system: what was irradiated, through which entrance channel, and at what energy.

The target and the entrance channel are what a configuration declares; everything else follows and
is **derived**, never declared. The reaction code comes from the channel through
[`CHANNEL_REACTION`](@ref), `A₀` and `Z₀` are therefore the compound nucleus for neutron-induced
fission and the parent for spontaneous fission and cannot contradict the target they came from,
and `label` is derived by [`system_label`](@ref).

Incident energy is carried even though the extraction does not yet use it, because it is part of
what identifies a system: the same target and channel at two energies are two systems, and a
downstream code matching on the label alone would conflate them.

# Fields

- `target_A`, `target_Z`: the nuclide irradiated, or the one that fissions spontaneously.
- `channel`: the entrance channel, one of the keys of [`CHANNEL_REACTION`](@ref).
- `reaction`: the reaction code the channel implies, one of the keys of [`REACTIONS`](@ref);
  derived.
- `incident_energy`: in MeV; zero for spontaneous fission.
- `A₀`, `Z₀`: the fissioning nucleus, derived.
- `label`: the canonical system identifier, derived.
"""
struct SystemSpecification
    target_A::Int
    target_Z::Int
    channel::String
    reaction::String
    incident_energy::Float64
    A₀::Int
    Z₀::Int
    label::String
end

"""
    REACTIONS

The reactions a fissioning system can be formed by, and the number of neutrons each absorbs:
spontaneous fission, `"0,f"`, and neutron-induced fission, `"n,f"`. The compound nucleus follows —
neutron-induced fission adds one mass unit to the target, spontaneous fission none — so `A₀` and
`Z₀` are derived rather than declared, and cannot disagree with the target they came from.
"""
const REACTIONS = Dict("0,f" => 0, "n,f" => 1)

"""
    CHANNEL_REACTION

The reaction code each entrance channel is formed through: spontaneous fission `sf`, and
neutron-induced fission from the thermal, resonance and fast regions, `nth`, `nres` and `nfast`.

The channel is what a configuration declares, because the reaction code alone cannot tell a
thermal run of a system from a resonance run of the same system — both are `"n,f"` — while the
system identifier must. The incident-energy window still does the selecting, and the channel has
to agree with it; the boundaries between the regions are conventions rather than constants, so
that agreement is a matter of stating the channel correctly and is not enforced here.
"""
const CHANNEL_REACTION = Dict("sf" => "0,f", "nth" => "n,f", "nres" => "n,f", "nfast" => "n,f")

# Element symbols by proton number, for deriving the system identifier from the target.
const ELEMENT_SYMBOLS = split(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As \
     Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd \
     Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am \
     Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og",
)

"""
    element_symbol(Z) -> String

Chemical symbol of the element with proton number `Z`.
"""
function element_symbol(Z::Integer)
    1 ≤ Z ≤ length(ELEMENT_SYMBOLS) || throw(ArgumentError("no element symbol for Z = $(Z)"))
    return String(ELEMENT_SYMBOLS[Z])
end

"""
    system_label(target_A, target_Z, channel) -> String

The canonical identifier of a fissioning system, `<symbol><target A>_<channel>` — `Cf252_sf`,
`U235_nth`. It names the data directory, the configuration file and the output subdirectory, and
it is the key a downstream code matches on.

Derived from the target and the channel rather than declared, so that two runs of one system
cannot be labelled differently and a label cannot contradict the nuclide it names.

This is the ASCII token; [`system_notation`](@ref) is the typeset form for a figure.

# Examples

```jldoctest
julia> system_label(252, 98, "sf")
"Cf252_sf"

julia> system_label(235, 92, "nth")
"U235_nth"
```
"""
function system_label(target_A::Integer, target_Z::Integer, channel::AbstractString)
    return string(element_symbol(target_Z), target_A, "_", channel)
end

const SUPERSCRIPT_DIGITS = Dict(
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
)

# The channel as the fission literature typesets it, inside the reaction parentheses.
const CHANNEL_NOTATION = Dict(
    "sf" => "sf", "nth" => "nth,f", "nres" => "nres,f", "nfast" => "nfast,f"
)

"""
    system_notation(target_A, target_Z, channel) -> String

The typeset form of a fissioning system, for a figure label or a caption — `²⁵²Cf(sf)`,
`²³³U(nth,f)`, with the mass number superscripted as the literature writes it.

Distinct from [`system_label`](@ref) on purpose: one name may not mean both the token that goes
into a path and the notation that goes into a figure.

# Examples

```jldoctest
julia> system_notation(252, 98, "sf")
"²⁵²Cf(sf)"

julia> system_notation(233, 92, "nth")
"²³³U(nth,f)"
```
"""
function system_notation(target_A::Integer, target_Z::Integer, channel::AbstractString)
    haskey(CHANNEL_NOTATION, channel) ||
        throw(ArgumentError("no notation for the entrance channel $(repr(channel))"))
    mass = map(digit -> SUPERSCRIPT_DIGITS[digit], string(target_A))
    return string(mass, element_symbol(target_Z), "(", CHANNEL_NOTATION[channel], ")")
end

"""
    system_notation(system) -> String

The typeset form of the system a [`SystemSpecification`](@ref) describes.
"""
function system_notation(system::SystemSpecification)
    return system_notation(system.target_A, system.target_Z, system.channel)
end

"""
    FragmentationSettings

Construction of the fragmentation range: how many charge numbers are taken per mass number, how
far the heavy-fragment mass range extends, and where the charge polarization and dispersion of
the isobaric charge distribution come from. Paths are relative to the data directory.
"""
struct FragmentationSettings
    charges_per_mass::Int
    heavy_mass_min::Int
    heavy_mass_max::Int
    charge_distribution_file::Union{String,Nothing}
    fallback_charge_polarization::Float64
    fallback_charge_dispersion::Float64
end

"""
    LevelDensitySettings

Which level density model to use, the tabulated data it is evaluated from, and the order in which
the parameter ratio of complementary fragments is averaged over the isobaric charge distribution.

The model is named rather than constructed here, because constructing it means reading its data;
[`build_level_density_model`](@ref) does that when a run starts.
"""
struct LevelDensitySettings
    model::Symbol
    mass_excess_file::String
    shell_correction_file::Union{String,Nothing}
    ratio_averaging::RatioAveraging
end

"""
    build_level_density_model(settings) -> LevelDensityModel

Read the tabulated data named by `settings` and construct the level density model from it.

Throws an `ArgumentError` if the model is unknown, or if the data it needs was not named.
[`load_configuration`](@ref) already rejects both, so this guards settings assembled by hand.
"""
function build_level_density_model(settings::LevelDensitySettings)
    if settings.model === :BSFG
        return BackShiftedFermiGas(read_mass_excess_table(settings.mass_excess_file))
    elseif settings.model === :GC
        path = settings.shell_correction_file
        path === nothing &&
            throw(ArgumentError("the Gilbert-Cameron model needs a shell correction file"))
        return GilbertCameron(read_shell_correction_table(path))
    end
    return throw(ArgumentError("unknown level density model: $(settings.model)"))
end

"""
    SegmentSettings

Controls of the piecewise-linear parameterization: the largest number of segments examined, the
smallest number of data points a segment may contain, whether the ratio is pinned to one half at
the symmetric split, and mass-number windows that must each contain a breakpoint.

`required_windows` places a breakpoint where physics says there is one — the minimum at the heavy
magic fragment, `A_H` near 130, where the `Z = 50`, `N = 82` shell closure fixes the sharing. It
constrains the systematic-trend curve by default and the per-dataset parameterizations only when
`windows_apply_to_datasets` is set, since a dataset that resolves the feature on its own should be
left to do so.

`parsimony` biases the choice of order towards fewer segments by multiplying the penalty the
selection criterion charges per parameter. One is the criterion as published.
"""
struct SegmentSettings
    max_segments::Int
    min_points_per_segment::Int
    pin_symmetric_split::Bool
    required_windows::Vector{UnitRange{Int}}
    windows_apply_to_datasets::Bool
    parsimony::Float64
end

"""
    OutputSettings

Rounding of tabulated output and the subdirectory, under the results and plots directories, that
a run writes into.

`significant_digits` is significant figures, not decimal places: an uncertainty of `1.2e-5` and a
ratio of `1.1` are both written to the same number of meaningful digits, which a fixed number of
decimals cannot do for two quantities of such different magnitude.
"""
struct OutputSettings
    significant_digits::Int
    subdirectory::String
end

"""
    Configuration

A complete, validated pipeline configuration. Construct with [`load_configuration`](@ref) rather
than directly, so that the constraints documented in the TOML file are enforced.
"""
struct Configuration
    system::SystemSpecification
    fragmentation::FragmentationSettings
    level_density::LevelDensitySettings
    multiplicity_directory::String
    excluded_datasets::Dict{String,String}
    yield_directory::Union{String,Nothing}
    segments::SegmentSettings
    output::OutputSettings
    source::String
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
has_symmetric_split(configuration::Configuration) = iseven(configuration.system.A₀)

const LEVEL_DENSITY_MODELS = Dict{String,Symbol}("BSFG" => :BSFG, "GC" => :GC)
const AVERAGINGS = Dict{String,RatioAveraging}(
    "ratio_of_means" => RatioOfMeans(), "mean_of_ratios" => MeanOfRatios()
)

function _section(document::AbstractDict, name::String, source::String)
    haskey(document, name) ||
        throw(ArgumentError("configuration $(source) is missing the [$(name)] section"))
    section = document[name]
    section isa AbstractDict ||
        throw(ArgumentError("configuration $(source): [$(name)] must be a table"))
    return section
end

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

function _one_of(value::String, options::AbstractDict, path::String)
    haskey(options, value) || throw(
        ArgumentError("$(path) must be one of $(join(sort(collect(keys(options))), " | ")), \
                       got $(repr(value))"),
    )
    return options[value]
end

# Datasets kept out of the pooling, each with the reason written down. A reason is required:
# excluding a measurement is a judgement, and an unexplained one is indistinguishable from a
# mistake to anyone reading the configuration later.
function _exclusions(section::AbstractDict, path::String)
    haskey(section, "exclude") || return Dict{String,String}()
    raw = section["exclude"]
    raw isa AbstractVector || throw(
        ArgumentError("$(path) must be an array of tables, each with `dataset` and `reason`"),
    )
    exclusions = Dict{String,String}()
    for (index, entry) in enumerate(raw)
        entry isa AbstractDict &&
        haskey(entry, "dataset") &&
        haskey(entry, "reason") &&
        entry["dataset"] isa String &&
        entry["reason"] isa String || throw(
            ArgumentError(
                "$(path)[$(index)] must be a table with string keys `dataset` and `reason`"
            ),
        )
        isempty(strip(entry["reason"])) &&
            throw(ArgumentError("$(path)[$(index)] must give a non-empty reason"))
        haskey(exclusions, entry["dataset"]) &&
            throw(ArgumentError("$(path) names $(repr(entry["dataset"])) more than once"))
        exclusions[entry["dataset"]] = entry["reason"]
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

Every key is checked for presence, type, enumerated choice and numerical range, and every input
file named by the configuration is checked for existence, before the pipeline is allowed to
start. Failures throw an `ArgumentError` naming the offending key, so that a configuration the
pipeline cannot honour never begins a run.

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

    system_section = _section(document, "system", source)
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
    channel = _value(system_section, "channel", String, "system.channel")
    reaction = _one_of(channel, CHANNEL_REACTION, "system.channel")
    neutrons_absorbed = REACTIONS[reaction]
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
    A₀ = target_A + neutrons_absorbed
    Z₀ = target_Z
    label = system_label(target_A, target_Z, channel)
    symmetric_split = cld(A₀, 2)

    fragmentation_section = _section(document, "fragmentation", source)
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
            min = heavy_mass_min,
            max = A₀ - 1,
        ),
    )
    fallback_ΔZ = Float64(
        _in_bounds(
            _value(
                fragmentation_section,
                "fallback_charge_polarization",
                Real,
                "fragmentation.fallback_charge_polarization",
                -0.5,
            ),
            "fragmentation.fallback_charge_polarization";
            min = -5,
            max = 5,
        ),
    )
    fallback_σ_Z = Float64(
        _in_bounds(
            _value(
                fragmentation_section,
                "fallback_charge_dispersion",
                Real,
                "fragmentation.fallback_charge_dispersion",
                0.6,
            ),
            "fragmentation.fallback_charge_dispersion";
            min = 0,
            max = 5,
            exclusive_min = true,
        ),
    )
    charge_file = _input_file(
        fragmentation_section,
        "charge_distribution_file",
        "fragmentation.charge_distribution_file",
        data_directory;
        required = false,
    )

    level_density_section = _section(document, "level_density", source)
    model = _one_of(
        _value(level_density_section, "model", String, "level_density.model"),
        LEVEL_DENSITY_MODELS,
        "level_density.model",
    )
    mass_excess_file = _input_file(
        level_density_section,
        "mass_excess_file",
        "level_density.mass_excess_file",
        data_directory,
    )
    shell_correction_file = _input_file(
        level_density_section,
        "shell_correction_file",
        "level_density.shell_correction_file",
        data_directory;
        required = false,
    )
    model === :GC &&
        shell_correction_file === nothing &&
        throw(ArgumentError("level_density.shell_correction_file is required when \
             level_density.model is \"GC\""))
    averaging = _one_of(
        _value(
            level_density_section,
            "ratio_averaging",
            String,
            "level_density.ratio_averaging",
            "ratio_of_means",
        ),
        AVERAGINGS,
        "level_density.ratio_averaging",
    )

    multiplicity_section = _section(document, "multiplicity", source)
    multiplicity_directory = _subdirectory(
        multiplicity_section, "multiplicity.subdirectory", data_directory
    )
    excluded_datasets = _exclusions(multiplicity_section, "multiplicity.exclude")

    # Optional. Without it the run reports the mean over the fragment mass range only; with it,
    # the total average over each yield distribution, which is the quantity the literature quotes.
    yield_directory = if haskey(document, "yield")
        _subdirectory(document["yield"], "yield.subdirectory", data_directory)
    else
        nothing
    end

    segments_section = _section(document, "segments", source)
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
    parsimony = Float64(
        _in_bounds(
            _value(segments_section, "parsimony", Real, "segments.parsimony", 1.0),
            "segments.parsimony";
            min = 0,
            max = 10,
            exclusive_min = true,
        ),
    )
    for (index, window) in enumerate(windows)
        issubset(window, heavy_mass_min:heavy_mass_max) || throw(
            ArgumentError("segments.required_windows[$(index)] = $(window) lies outside the \
                           fragmentation range $(heavy_mass_min):$(heavy_mass_max)"),
        )
    end

    output_section = _section(document, "output", source)
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
    subdirectory = _value(output_section, "subdirectory", String, "output.subdirectory", label)
    isempty(subdirectory) && throw(ArgumentError("output.subdirectory must not be empty"))

    return Configuration(
        SystemSpecification(
            target_A, target_Z, channel, reaction, incident_energy, A₀, Z₀, label
        ),
        FragmentationSettings(
            charges_per_mass,
            heavy_mass_min,
            heavy_mass_max,
            charge_file,
            fallback_ΔZ,
            fallback_σ_Z,
        ),
        LevelDensitySettings(model, mass_excess_file, shell_correction_file, averaging),
        multiplicity_directory,
        excluded_datasets,
        yield_directory,
        SegmentSettings(max_segments, min_points, pin, windows, windows_apply, parsimony),
        OutputSettings(significant_digits, subdirectory),
        source,
    )
end
