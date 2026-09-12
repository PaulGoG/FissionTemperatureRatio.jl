# Pipeline configuration: parsing and validation of the TOML input.

"""
    SystemSpecification

The fissioning nucleus. `label` identifies the case in output paths and figure legends; `A₀` and
`Z₀` are the mass and charge numbers of the nucleus undergoing fission, that is, of the compound
nucleus for neutron-induced fission and of the parent for spontaneous fission.
"""
struct SystemSpecification
    label::String
    A₀::Int
    Z₀::Int
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

`required_windows` constrains the systematic-trend curve only, not the per-data-set
parameterizations, which are left to follow their own data. Its purpose is to place the minimum
at the heavy magic fragment where a data set is too sparse or too scattered to resolve it.
"""
struct SegmentSettings
    max_segments::Int
    min_points_per_segment::Int
    pin_symmetric_split::Bool
    required_windows::Vector{UnitRange{Int}}
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
    label = _value(system_section, "label", String, "system.label")
    isempty(label) && throw(ArgumentError("system.label must not be empty"))
    A₀ = Int(
        _in_bounds(
            _value(system_section, "A0", Integer, "system.A0"), "system.A0"; min = 2, max = 400
        ),
    )
    Z₀ = Int(
        _in_bounds(
            _value(system_section, "Z0", Integer, "system.Z0"), "system.Z0"; min = 1, max = 120
        ),
    )
    Z₀ < A₀ || throw(ArgumentError("system.Z0 must be smaller than system.A0, \
                                    got Z0 = $(Z₀), A0 = $(A₀)"))

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
                "segments.pin_symmetric_split requires an even system.A0, since the ratio \
                       equals one half only where the two fragments are identical; got \
                       A0 = $(A₀)"
            ),
        )
    windows = _windows(segments_section, "segments.required_windows")
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
        SystemSpecification(label, A₀, Z₀),
        FragmentationSettings(
            charges_per_mass, A_H_max, charge_file, fallback_ΔZ, fallback_rms
        ),
        LevelDensitySettings(prescription, mass_excess_file, shell_correction_file, averaging),
        multiplicity_directory,
        SegmentSettings(max_segments, min_points, pin, windows),
        OutputSettings(digits, subdirectory),
        source,
    )
end
