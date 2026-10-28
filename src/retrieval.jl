# The records the EXFOR retrieval writes beside the files it produced: where every input dataset
# came from, read once per directory and checked for what this package takes from them.

"""
    FLAGGED_QUALIFIERS

Reaction-code qualifiers that mark a dataset as other than a direct measurement of the quantity
at the entrance channel's energy: `DERIV`, derived from other data, and `SPA`, averaged over an
unspecified neutron spectrum; and `preliminary`, which marks a dataset whose publication calls
its results preliminary. A dataset carrying any of them is used, but flagged in the log, in the
dataset diagnostics and in the run metadata.
"""
const FLAGGED_QUALIFIERS = ("DERIV", "SPA", "preliminary")

# An EXFOR dataset identifier as the manifest admits it: entry, subentry and, where the subentry
# holds several datasets, a pointer character.
const ACCESSION_PATTERN = r"^[0-9A-Z]{8,9}$"

"""
    PairSumScale

The scale of a multiplicity dataset as the retrieval states it for one read by the complement
test: `deviation`, the relative departure of its pair sum `ν(A) + ν(A₀ - A)`, weighted with the
light-fragment yield, from `ν̄`; its `uncertainty`; and `consistent`, whether the two agree within
three standard deviations. `uncertainty` and `consistent` are `missing` where the dataset states
no uncertainty.

A uniform scale cancels in `r_ν = ν_H/(ν_L + ν_H)`, so a dataset off the scale of `ν̄` is used
and flagged, not corrected. A scale error that varies with mass does not cancel, and the pair sum
cannot detect it; the structural diagnostics of the dataset are what bear on it.
"""
struct PairSumScale
    deviation::Float64
    uncertainty::Union{Float64,Missing}
    consistent::Union{Bool,Missing}
end

"""
    RetrievalRecord

The entry of one dataset in a `retrieval*.toml` the EXFOR retrieval writes beside the files it
produced, with the run of that retrieval, read by [`retrieval_records`](@ref).

# Fields

- `path`: the record file.
- `accession`: the EXFOR dataset identifier the entry states as `identifier`; empty where it
  states none.
- `qualifiers`: the reaction-code qualifiers recorded for the dataset, as tags (`"DERIV"`,
  `"MXW"`, …).
- `measured_fraction`: the share of a measurement each row of the file amounts to: raw points
  over rows written, `mass_values_non_integer / rows_written`, for a dataset whose retrieval
  interpolated non-integer masses onto the integers (`mass_treatment = "interpolated"`), at most
  one, and one otherwise. Neighbouring interpolated rows share their bracketing points, so a
  dataset written at more masses than it measured carries no more information than its
  measurements; one written at fewer, several raw masses falling between the same two integers,
  has a measurement behind every row.
- `correlated_with`: the EXFOR identifiers of the other datasets of the experiment the dataset
  belongs to, runs or analyses on one apparatus, whose values are not independent measurements.
- `correlation_relation`: how it is related to them — `"republication"`, `"alternative_analysis"`,
  `"repeated_run"`, `"complementary_range"`, … — or empty where the entry states none.
- `pair_sum`: the [`PairSumScale`](@ref) the entry states, or `nothing` where it forms none.
- `entry`: the dataset's `[[accepted]]` table as written, for the provenance record.
- `run`: the record's `[run]` table: configuration, parser version and revision, timestamp;
  empty for a record that carries none.
"""
struct RetrievalRecord
    path::String
    accession::String
    qualifiers::Vector{String}
    measured_fraction::Float64
    correlated_with::Vector{String}
    correlation_relation::String
    pair_sum::Union{PairSumScale,Nothing}
    entry::Dict{String,Any}
    run::Dict{String,Any}
end

# Every `retrieval*.toml` in a directory: one per query that contributed files to it.
function _retrieval_documents(directory::AbstractString)
    isdir(directory) || return Pair{String,Dict{String,Any}}[]
    names = sort!(
        filter(f -> startswith(f, "retrieval") && endswith(f, ".toml"), readdir(directory))
    )
    return [
        joinpath(directory, name) => TOML.parsefile(joinpath(directory, name)) for name in names
    ]
end

# "MXW: Maxwellian-averaged" → "MXW"; the retrieval writes the tag and its meaning.
_qualifier_tag(qualifier::AbstractString) = String(strip(first(split(qualifier, ':'))))

# A number that is not a Bool: TOML reads `true` as a Bool, which Julia counts among the integers.
_is_count(value) = value isa Integer && !(value isa Bool)
_is_number(value) = value isa Real && !(value isa Bool)

# One `[[accepted]]` entry, checked for every key this package reads. A key of the wrong type or
# out of range is a defect of the record and is refused, naming the record and the file, rather
# than read as absent: a count that is not an integer would otherwise pool an interpolated
# dataset at full weight, and a malformed `correlated_with` would pool one experiment twice.
function _record_entry(path::AbstractString, entry::AbstractDict, run::AbstractDict)
    file = entry["file"]
    defect(what) =
        throw(ArgumentError("the retrieval record $(path), entry $(repr(file)): $(what)"))

    accession = get(entry, "identifier", "")
    accession isa AbstractString ||
        defect("identifier must be a string, got $(repr(accession))")
    isempty(accession) ||
        occursin(ACCESSION_PATTERN, accession) ||
        defect("identifier $(repr(accession)) is not an EXFOR dataset identifier")

    listed = get(entry, "qualifiers", String[])
    listed isa AbstractVector && all(q -> q isa AbstractString, listed) ||
        defect("qualifiers must be an array of strings, got $(repr(listed))")

    treatment = get(entry, "mass_treatment", "")
    treatment isa AbstractString ||
        defect("mass_treatment must be a string, got $(repr(treatment))")
    measured_fraction = 1.0
    if treatment == "interpolated"
        points = get(entry, "mass_values_non_integer", nothing)
        rows = get(entry, "rows_written", nothing)
        _is_count(points) && _is_count(rows) && points > 0 && rows > 0 || defect(
            "an interpolated dataset states its raw points and rows written as positive \
             integers, mass_values_non_integer and rows_written, got $(repr(points)) and \
             $(repr(rows))"
        )
        measured_fraction = min(points / rows, 1.0)
    end

    partners = get(entry, "correlated_with", String[])
    partners isa AbstractVector &&
    all(p -> p isa AbstractString && occursin(ACCESSION_PATTERN, p), partners) ||
        defect("correlated_with must be an array of EXFOR dataset identifiers, got \
                 $(repr(partners))")
    relation = get(entry, "correlation_relation", "")
    relation isa AbstractString ||
        defect("correlation_relation must be a string, got $(repr(relation))")

    pair_sum = if haskey(entry, "pair_sum_deviation")
        deviation = entry["pair_sum_deviation"]
        uncertainty = get(entry, "pair_sum_deviation_uncertainty", missing)
        consistent = get(entry, "scale_consistent", missing)
        _is_number(deviation) || defect("pair_sum_deviation must be a number")
        ismissing(uncertainty) ||
            _is_number(uncertainty) ||
            defect("pair_sum_deviation_uncertainty must be a number")
        ismissing(consistent) ||
            consistent isa Bool ||
            defect("scale_consistent must be true or false")
        PairSumScale(
            Float64(deviation),
            ismissing(uncertainty) ? missing : Float64(uncertainty),
            consistent,
        )
    else
        nothing
    end

    return RetrievalRecord(
        String(path),
        String(accession),
        String[_qualifier_tag(q) for q in listed],
        measured_fraction,
        String[String(p) for p in partners],
        String(relation),
        pair_sum,
        Dict{String,Any}(entry),
        Dict{String,Any}(run),
    )
end

"""
    retrieval_records(directory) -> Dict{String,RetrievalRecord}

The [`RetrievalRecord`](@ref) of every data file the `retrieval*.toml` records in `directory`
list, keyed by file name; empty where the directory holds no record. Each record is read once.

Throws an `ArgumentError` naming the record where an entry states a key this package reads with
the wrong type or out of range, where a `[run]` table is not a table, and where two entries list
one file.
"""
function retrieval_records(directory::AbstractString)
    records = Dict{String,RetrievalRecord}()
    for (path, document) in _retrieval_documents(directory)
        run = get(document, "run", Dict{String,Any}())
        run isa AbstractDict ||
            throw(ArgumentError("the retrieval record $(path) has a [run] that is not a table"))
        accepted = get(document, "accepted", Any[])
        accepted isa AbstractVector || throw(
            ArgumentError(
                "the retrieval record $(path) has an `accepted` that is not an array"
            ),
        )
        for entry in accepted
            entry isa AbstractDict && get(entry, "file", nothing) isa AbstractString || throw(
                ArgumentError("the retrieval record $(path) has an accepted entry without a \
                               file name"),
            )
            file = entry["file"]
            haskey(records, file) && throw(
                ArgumentError("$(file) is listed by both $(records[file].path) and $(path)")
            )
            records[file] = _record_entry(path, entry, run)
        end
    end
    return records
end

"""
    retrieval_record(path) -> Union{RetrievalRecord,Nothing}

The retrieval record of the data file at `path`, from the records of its directory, or `nothing`
where none there lists the file. Reads every record of the directory; a caller holding several
files of one directory takes [`retrieval_records`](@ref) once instead.
"""
function retrieval_record(path::AbstractString)
    return get(retrieval_records(dirname(path)), basename(path), nothing)
end

"""
    check_retrieval_versions(directories, floor) -> Dict{String,VersionNumber}

The parser version, `[run] package_version`, of every `retrieval*.toml` in `directories`, keyed by
record path. A release of the retrieval can change what a dataset holds, its scale or whether it
is admitted at all, so inputs written before the version a run was validated on must not re-enter
silently.

Throws an `ArgumentError` naming the record for one that states no version or one below `floor`;
warns when the records were written by different versions. A directory holding no record, a
tabulation not produced by the retrieval, is not judged.
"""
function check_retrieval_versions(directories, floor::VersionNumber)
    versions = Dict{String,VersionNumber}()
    for directory in unique(directories)
        for (record, document) in _retrieval_documents(directory)
            run = get(document, "run", nothing)
            stated = run isa AbstractDict ? get(run, "package_version", nothing) : nothing
            version = stated isa AbstractString ? tryparse(VersionNumber, stated) : nothing
            version === nothing && throw(
                ArgumentError(
                    "the retrieval record $(record) states no [run] package_version, so it \
                     cannot be held to retrieval.min_package_version = \"$(floor)\"; retrieve \
                     the inputs again with ExforFissionData $(floor) or later",
                ),
            )
            version < floor && throw(
                ArgumentError(
                    "the retrieval record $(record) was written by ExforFissionData \
                     $(version), below retrieval.min_package_version = \"$(floor)\"; retrieve \
                     the inputs again with $(floor) or later",
                ),
            )
            versions[record] = version
        end
    end
    length(unique(values(versions))) > 1 &&
        @warn "the input directories of this run were written by different versions of the \
               retrieval" versions = Dict(record => string(v) for (record, v) in versions)
    return versions
end

# The EXFOR accession of a dataset: the identifier its retrieval record states, else the one its
# file name leads with, else none, for a tabulation that came from no archive.
function _accession(record::Union{RetrievalRecord,Nothing}, file::AbstractString)
    record === nothing || isempty(record.accession) || return record.accession
    candidate = first(split(basename(file), '_'))
    return occursin(ACCESSION_PATTERN, candidate) ? String(candidate) : ""
end

_qualifiers(::Nothing) = String[]
_qualifiers(record::RetrievalRecord) = record.qualifiers

_measured_fraction(::Nothing) = 1.0
_measured_fraction(record::RetrievalRecord) = record.measured_fraction

_pair_sum(::Nothing) = nothing
_pair_sum(record::RetrievalRecord) = record.pair_sum

_correlated_with(::Nothing) = String[]
_correlated_with(record::RetrievalRecord) = record.correlated_with

_correlation_relation(::Nothing) = ""
_correlation_relation(record::RetrievalRecord) = record.correlation_relation

function _flagged(qualifiers::AbstractVector{<:AbstractString})
    return [q for q in qualifiers if q in FLAGGED_QUALIFIERS]
end
