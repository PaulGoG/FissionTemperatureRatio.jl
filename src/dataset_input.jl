# The multiplicity datasets of a directory as read: their labels, their accessions and the
# exclusions a configuration names among them.

"""
    read_multiplicity_directory(directory) -> Vector{Multiplicity}

Read every `.dat` file in `directory`, sorted by name so that colours and markers are assigned
reproducibly. Other files are ignored, so a run record can sit beside the data it describes.
`records` are the retrieval records of the directory, read from it by default; a caller that
holds them passes them.

The label of each dataset is its file name stripped of the leading archive identifier and the
extension, with underscores replaced by spaces, which is how the measurements are named in the
literature. Two files of one author and year — two subentries of one measurement — would share
that label, and a label selects a curve downstream, so each of them carries its archive identifier
in parentheses instead.

Throws an `ArgumentError` where the directory does not exist, holds no data file, or holds two
files of one EXFOR accession: one measurement would be read, fitted and pooled twice.
"""
function read_multiplicity_directory(
    directory::AbstractString;
    records::AbstractDict{String,RetrievalRecord} = retrieval_records(directory),
)
    isdir(directory) || throw(ArgumentError("multiplicity directory not found: $(directory)"))
    files = _data_files(directory)
    isempty(files) &&
        throw(ArgumentError("multiplicity directory holds no data files: $(directory)"))
    _refuse_shared_accession(directory, files, records)
    return [
        read_multiplicity(
            joinpath(directory, file); label = label, record = get(records, file, nothing)
        ) for (file, label) in zip(files, _unique_labels(files))
    ]
end

# Data files carry the `.dat` extension; anything else in the directory — a run record, a note —
# is not data and is skipped rather than parsed and rejected.
function _data_files(directory::AbstractString)
    return sort!(filter!(f -> endswith(f, ".dat") && !startswith(f, "."), readdir(directory)))
end

# `<accession>_<Author>_<year>.dat`. The accession has eight characters, or nine where the dataset
# is a pointer within a subentry; entries outside the neutron areas of EXFOR lead with a letter.
function _dataset_label(file::AbstractString)
    stem = replace(splitext(file)[1], ACCESSION_FILE_PREFIX => "")
    label = replace(stem, '_' => ' ')
    # Initials are written without a space in the file names; restore it for the legend.
    return replace(label, r"(?<=\b[A-Z])\.(?=[A-Z][a-z])" => ". ")
end

function _accession(file::AbstractString)
    found = match(ACCESSION_FILE_PREFIX, file)
    return found === nothing ? splitext(file)[1] : found.captures[1]
end

# The labels of the data files of one directory, the archive identifier in parentheses after each
# label two files share.
function _unique_labels(files::AbstractVector{<:AbstractString})
    labels = _dataset_label.(files)
    return [
        count(==(label), labels) == 1 ? label : "$(label) ($(_accession(file)))" for
        (file, label) in zip(files, labels)
    ]
end

_accession(data::Multiplicity) = _accession(data.record, data.source)
_accession(energies::MeanKineticEnergy) = _accession(energies.record, energies.source)

# Two files of one directory carrying one accession are one measurement held twice — a leftover
# of an earlier retrieval under another spelling of the author, a copy — and would be read, fitted
# and pooled twice, while an exclusion by that accession could name only one of them.
function _refuse_shared_accession(
    directory::AbstractString, files, records::AbstractDict{String,RetrievalRecord}
)
    held = Dict{String,String}()
    for file in files
        accession = _accession(get(records, file, nothing), file)
        isempty(accession) && continue
        haskey(held, accession) && throw(
            ArgumentError(
                "$(directory) holds two files of the accession $(accession), \
                 $(held[accession]) and $(file); one measurement would be read twice"
            ),
        )
        held[accession] = file
    end
    return nothing
end

# The exclusion key of every data file of a directory, with its label, from the file names and the
# retrieval record alone.
function _directory_keys(
    directory::AbstractString, records::AbstractDict{String,RetrievalRecord}
)
    isdir(directory) || throw(ArgumentError("directory not found: $(directory)"))
    files = _data_files(directory)
    _refuse_shared_accession(directory, files, records)
    return Dict(
        _exclusion_key(get(records, file, nothing), file, label) => label for
        (file, label) in zip(files, _unique_labels(files))
    )
end

# The key a dataset is excluded under: its accession, or its label where it carries none.
function _exclusion_key(
    record::Union{RetrievalRecord,Nothing}, source::AbstractString, label::AbstractString
)
    accession = _accession(record, source)
    return isempty(accession) ? String(label) : accession
end
_exclusion_key(data::Multiplicity) = _exclusion_key(data.record, data.source, data.label)
function _exclusion_key(data::MassYield, records::AbstractDict{String,RetrievalRecord})
    return _exclusion_key(get(records, data.label, nothing), data.source, data.label)
end

# The exclusions that apply to the datasets `keyed` lists as exclusion key => label, by label,
# each with its reason. One that names no dataset read is an error: the configuration would claim
# an exclusion the run did not make.
function _excluded_labels(
    exclusions::Dict{String,String}, keyed::AbstractVector, path::AbstractString
)
    keys_read = Dict{String,String}()
    for (key, label) in keyed
        haskey(keys_read, key) && throw(
            ArgumentError("$(repr(keys_read[key])) and $(repr(label)) are both held under \
                 $(repr(key)); an exclusion could not tell them apart")
        )
        keys_read[key] = label
    end
    for key in sort!(collect(keys(exclusions)))
        haskey(keys_read, key) || throw(
            ArgumentError(
                "$(path) excludes $(repr(key)), which names no dataset read; an exclusion \
                 is keyed by the EXFOR accession of its dataset"
            ),
        )
    end
    return Dict(keys_read[key] => reason for (key, reason) in exclusions)
end
function _excluded_labels(exclusions::Dict{String,String}, datasets::Vector{Multiplicity})
    return _excluded_labels(
        exclusions, [_exclusion_key(d) => d.label for d in datasets], "multiplicity.exclude"
    )
end
function _excluded_labels(
    exclusions::Dict{String,String},
    distributions::Vector{MassYield},
    records::AbstractDict{String,RetrievalRecord},
)
    return _excluded_labels(
        exclusions,
        [_exclusion_key(y, records) => y.label for y in distributions],
        "yield.exclude",
    )
end

# One field of the pair-sum scale the retrieval records for a dataset, or `missing`.
function _scale_field(data::Multiplicity, field::Symbol, digits::Integer)
    scale = _pair_sum(data.record)
    scale === nothing && return missing
    value = getfield(scale, field)
    return value isa AbstractFloat ? round(value; sigdigits = digits) : value
end

# The retrieval record of each yield distribution, keyed by label, where one lists it; the
# records of each directory are read once.
function _records_by_label(mass_yields::Vector{MassYield}, retrievals::RetrievalDocuments)
    by_directory = Dict(
        directory => retrieval_records(_documents(retrievals, directory)) for
        directory in unique(dirname(y.source) for y in mass_yields)
    )
    records = Dict{String,RetrievalRecord}()
    for y in mass_yields
        record = get(by_directory[dirname(y.source)], basename(y.source), nothing)
        record === nothing || (records[y.label] = record)
    end
    return records
end
