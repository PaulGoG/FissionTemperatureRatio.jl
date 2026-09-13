# Pipeline configuration: parsing and validation of the TOML input.

"""
    SystemSpecification

The fissioning system: what was irradiated, with what, and at what energy.

The target and the reaction are what a configuration declares; the nucleus that actually undergoes
fission follows from them and is **derived**, never declared. `A₀` and `Z₀` are therefore the
compound nucleus for neutron-induced fission and the parent for spontaneous fission, and cannot
contradict the target they came from. `label` is derived too, by [`case_label`](@ref).

Incident energy is carried even though the extraction does not yet use it, because it is part of
what identifies a system: the same target and reaction at two energies are two systems, and a
downstream code matching on the label alone would conflate them.

# Fields

- `target_A`, `target_Z`: the nuclide irradiated, or the one that fissions spontaneously.
- `reaction`: one of [`REACTIONS`](@ref).
- `incident_energy`: in MeV; zero for spontaneous fission.
- `A₀`, `Z₀`: the fissioning nucleus, derived.
- `label`: the canonical case identifier, derived.
"""
struct SystemSpecification
    target_A::Int
    target_Z::Int
    reaction::String
    incident_energy::Float64
    A₀::Int
    Z₀::Int
    label::String
end

"""
    REACTIONS

The reactions a fissioning system can be formed by: spontaneous fission, `"0,f"`, and
neutron-induced fission, `"n,f"`. The compound nucleus follows — neutron-induced fission adds one
mass unit to the target, spontaneous fission none — so `A₀` and `Z₀` are derived rather than
declared, and cannot disagree with the target they came from.
"""
const REACTIONS = Dict("0,f" => 0, "n,f" => 1)

# Element symbols by proton number, for deriving the case label from the target.
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
    case_label(system) -> String

The canonical identifier of a fissioning system, `<symbol><target A>_<reaction>` with the comma
dropped — `Cf252_0f`, `U235_nf`. Derived from the target and the reaction rather than declared, so
that two runs of the same system cannot be labelled differently, and a label cannot contradict the
nuclide it names. It is the key a downstream code matches on.
"""
function case_label(target_A::Integer, target_Z::Integer, reaction::AbstractString)
    return string(element_symbol(target_Z), target_A, "_", replace(reaction, "," => ""))
end

"""
    FragmentationSettings

Construction of the fragmentation range: how many charge numbers are taken per mass number, how
far the heavy-fragment mass range extends, and where the charge polarization and dispersion of
the isobaric charge distribution come from. Paths are relative to the data directory.
"""
struct FragmentationSettings
    charges_per_mass::Int
    A_H_max::Int
    charge_distribution_file::Union{String,Nothing}
    fallback_charge_polarization::Float64
    fallback_rms::Float64
end

"""
    LevelDensitySettings

Which level density prescription to use, the tabulated data it is evaluated from, and the order in
which the parameter ratio of complementary fragments is averaged over the isobaric charge
distribution.

The prescription is named rather than constructed here, because constructing it means reading its
data; [`build_prescription`](@ref) does that when a run starts.
"""
struct LevelDensitySettings
    prescription::Symbol
    mass_excess_file::String
    shell_correction_file::Union{String,Nothing}
    ratio_averaging::RatioAveraging
end

"""
    build_prescription(settings) -> LevelDensityPrescription

Read the tabulated data named by `settings` and construct the level density prescription from it.

Throws an `ArgumentError` if the prescription is unknown, or if the data it needs was not named.
[`load_configuration`](@ref) already rejects both, so this guards settings assembled by hand.
"""
function build_prescription(settings::LevelDensitySettings)
    if settings.prescription === :BSFG
        return BackShiftedFermiGas(read_mass_excess(settings.mass_excess_file))
    elseif settings.prescription === :GC
        path = settings.shell_correction_file
        path === nothing && throw(
            ArgumentError("the Gilbert-Cameron prescription needs a shell correction file")
        )
        return GilbertCameron(read_shell_corrections(path))
    end
    return throw(ArgumentError("unknown level density prescription: $(settings.prescription)"))
end

"""
    SegmentSettings

Controls of the piecewise-linear parameterization: the largest number of segments examined, the
smallest number of data points a segment may contain, whether the ratio is pinned to one half at
the symmetric split, and mass-number windows that must each contain a breakpoint.

`required_windows` places a breakpoint where physics says there is one — the minimum at the heavy
magic fragment, `A_H` near 130, where the `Z = 50`, `N = 82` shell closure fixes the sharing. It
constrains the systematic-trend curve by default and the per-data-set parameterizations only when
`windows_apply_to_data_sets` is set, since a set that resolves the feature on its own should be
left to do so.

`parsimony` biases the choice of order towards fewer segments by multiplying the penalty the
selection criterion charges per parameter. One is the criterion as published.
"""
struct SegmentSettings
    max_segments::Int
    min_points_per_segment::Int
    pin_symmetric_split::Bool
    required_windows::Vector{UnitRange{Int}}
    windows_apply_to_data_sets::Bool
    parsimony::Float64
end

"""
    OutputSettings

Rounding of tabulated output and the subdirectory, under the results and plots directories, that
a run writes into.
"""
struct OutputSettings
    digits::Int
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
    excluded_sets::Dict{String,String}
    yield_directory::Union{String,Nothing}
    segments::SegmentSettings
    output::OutputSettings
    source::String
end

"""
    A_H_min(configuration) -> Int

Smallest heavy-fragment mass number, the symmetric split rounded up. For an odd `A₀` there is no
symmetric split and this is the lighter member of the most symmetric pair.
"""
A_H_min(configuration::Configuration) = cld(configuration.system.A₀, 2)

"""
    A_H_range(configuration) -> UnitRange{Int}

Heavy-fragment mass numbers spanned by the fragmentation range.
"""
function A_H_range(configuration::Configuration)
    return A_H_min(configuration):configuration.fragmentation.A_H_max
end

"""
    has_symmetric_split(configuration) -> Bool

Whether the fissioning nucleus has a symmetric split, that is, whether `A₀` is even. The
multiplicity ratio equals one half there by the identity of the two fragments, which is what
pinning the parameterization relies on.
"""
has_symmetric_split(configuration::Configuration) = iseven(configuration.system.A₀)

const PRESCRIPTIONS = Dict{String,Symbol}("BSFG" => :BSFG, "GC" => :GC)
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

# Data sets kept out of the pooling, each with the reason written down. A reason is required:
# excluding a measurement is a judgement, and an unexplained one is indistinguishable from a
# mistake to anyone reading the configuration later.
function _exclusions(section::AbstractDict, path::String)
    haskey(section, "exclude") || return Dict{String,String}()
    raw = section["exclude"]
    raw isa AbstractVector ||
        throw(ArgumentError("$(path) must be an array of tables, each with `set` and `reason`"))
    exclusions = Dict{String,String}()
    for (index, entry) in enumerate(raw)
        entry isa AbstractDict &&
        haskey(entry, "set") &&
        haskey(entry, "reason") &&
        entry["set"] isa String &&
        entry["reason"] isa String || throw(
            ArgumentError(
                "$(path)[$(index)] must be a table with string keys `set` and `reason`"
            ),
        )
        isempty(strip(entry["reason"])) &&
            throw(ArgumentError("$(path)[$(index)] must give a non-empty reason"))
        haskey(exclusions, entry["set"]) &&
            throw(ArgumentError("$(path) names $(repr(entry["set"])) more than once"))
        exclusions[entry["set"]] = entry["reason"]
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

"""
    load_configuration(path; data_directory = datadir()) -> Configuration

Read and validate a pipeline configuration.

Every key is checked for presence, type, enumerated choice and numerical range, and every input
file named by the configuration is checked for existence, before the pipeline is allowed to
start. Failures throw an `ArgumentError` naming the offending key, so that a configuration the
pipeline cannot honour never begins a run.

`data_directory` is where the relative paths in the configuration are resolved against; it exists
so that tests can point at a fixture directory.

# Example

```julia
configuration = load_configuration(joinpath(projectdir(), "config", "U233_nf.toml"))
```
"""
function load_configuration(path::AbstractString; data_directory::AbstractString = datadir())
    isfile(path) || throw(ArgumentError("configuration file not found: $(path)"))
    document = try
        TOML.parsefile(path)
    catch err
        throw(ArgumentError("configuration $(path) is not valid TOML: $(err)"))
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
    reaction = _value(system_section, "reaction", String, "system.reaction")
    neutrons_absorbed = _one_of(reaction, REACTIONS, "system.reaction")
    incident_energy = Float64(
        _in_bounds(
            _value(system_section, "incident_energy", Real, "system.incident_energy", 0.0),
            "system.incident_energy";
            min = 0,
        ),
    )
    # Spontaneous fission has no incident particle, so an energy for it is a contradiction rather
    # than a harmless extra.
    reaction == "0,f" &&
        incident_energy != 0 &&
        throw(ArgumentError("system.incident_energy must be zero for spontaneous fission, \
             got $(incident_energy) MeV with reaction \"0,f\""))
    A₀ = target_A + neutrons_absorbed
    Z₀ = target_Z
    label = case_label(target_A, target_Z, reaction)

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
    A_H_max = Int(
        _in_bounds(
            _value(fragmentation_section, "A_H_max", Integer, "fragmentation.A_H_max"),
            "fragmentation.A_H_max";
            min = 2,
            max = A₀ - 1,
        ),
    )
    A_H_max ≥ cld(A₀, 2) ||
        throw(ArgumentError("fragmentation.A_H_max must be at least the symmetric split \
                       $(cld(A₀, 2)), got $(A_H_max)"))
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
    fallback_rms = Float64(
        _in_bounds(
            _value(
                fragmentation_section, "fallback_rms", Real, "fragmentation.fallback_rms", 0.6
            ),
            "fragmentation.fallback_rms";
            min = 0,
            max = 5,
            exclusive_min = true,
        ),
    )
    charge_file = let
        relative = _value(
            fragmentation_section,
            "charge_distribution_file",
            String,
            "fragmentation.charge_distribution_file",
            "",
        )
        if isempty(relative)
            nothing
        else
            resolved = joinpath(data_directory, relative)
            isfile(resolved) || throw(
                ArgumentError(
                    "fragmentation.charge_distribution_file does not exist: $(resolved)"
                ),
            )
            resolved
        end
    end

    level_density_section = _section(document, "level_density", source)
    prescription = _one_of(
        _value(level_density_section, "prescription", String, "level_density.prescription"),
        PRESCRIPTIONS,
        "level_density.prescription",
    )
    mass_excess_file = joinpath(
        data_directory,
        _value(
            level_density_section, "mass_excess_file", String, "level_density.mass_excess_file"
        ),
    )
    isfile(mass_excess_file) || throw(
        ArgumentError("level_density.mass_excess_file does not exist: $(mass_excess_file)")
    )
    shell_correction_file = let
        relative = _value(
            level_density_section,
            "shell_correction_file",
            String,
            "level_density.shell_correction_file",
            "",
        )
        if isempty(relative)
            prescription === :GC && throw(
                ArgumentError("level_density.shell_correction_file is required when \
                     level_density.prescription is \"GC\"")
            )
            nothing
        else
            resolved = joinpath(data_directory, relative)
            isfile(resolved) || throw(
                ArgumentError(
                    "level_density.shell_correction_file does not exist: $(resolved)"
                ),
            )
            resolved
        end
    end
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
    multiplicity_directory = joinpath(
        data_directory,
        _value(multiplicity_section, "directory", String, "multiplicity.directory"),
    )
    isdir(multiplicity_directory) ||
        throw(ArgumentError("multiplicity.directory does not exist: $(multiplicity_directory)"))
    excluded_sets = _exclusions(multiplicity_section, "multiplicity.exclude")

    # Optional. Without it the run reports the mean over the fragment mass range only; with it,
    # the total average over each yield distribution, which is the quantity the literature quotes.
    yield_directory = if haskey(document, "yield")
        directory = joinpath(
            data_directory,
            _value(document["yield"], "directory", String, "yield.directory"),
        )
        isdir(directory) ||
            throw(ArgumentError("yield.directory does not exist: $(directory)"))
        directory
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
    pin &&
        isodd(A₀) &&
        throw(
            ArgumentError(
                "segments.pin_symmetric_split requires an even fissioning mass number, since the ratio \
                       equals one half only where the two fragments are identical; got \
                       A0 = $(A₀)",
            ),
        )
    windows = _windows(segments_section, "segments.required_windows")
    windows_apply = _value(
        segments_section,
        "windows_apply_to_data_sets",
        Bool,
        "segments.windows_apply_to_data_sets",
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
        issubset(window, cld(A₀, 2):A_H_max) || throw(
            ArgumentError("segments.required_windows[$(index)] = $(window) lies outside the \
                           fragmentation range $(cld(A₀, 2)):$(A_H_max)"),
        )
    end

    output_section = _section(document, "output", source)
    digits = Int(
        _in_bounds(
            _value(output_section, "digits", Integer, "output.digits", 6),
            "output.digits";
            min = 1,
            max = 15,
        ),
    )
    subdirectory = _value(output_section, "subdirectory", String, "output.subdirectory", label)
    isempty(subdirectory) && throw(ArgumentError("output.subdirectory must not be empty"))

    return Configuration(
        SystemSpecification(target_A, target_Z, reaction, incident_energy, A₀, Z₀, label),
        FragmentationSettings(
            charges_per_mass, A_H_max, charge_file, fallback_ΔZ, fallback_rms
        ),
        LevelDensitySettings(prescription, mass_excess_file, shell_correction_file, averaging),
        multiplicity_directory,
        excluded_sets,
        yield_directory,
        SegmentSettings(max_segments, min_points, pin, windows, windows_apply, parsimony),
        OutputSettings(digits, subdirectory),
        source,
    )
end
