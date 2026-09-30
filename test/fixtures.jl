# Shared test inputs.
#
# The mass table and the Gilbert-Cameron shell corrections ship with FissionFragmentsDomain, so the
# fragmentation domain, both level density models and every inversion on them are tested on a bare
# clone. The experimental input data is
# held locally and is not shipped, so tests of a run on real measurements are skipped when it is
# absent; `data/README.md` records what the files are and where they come from.

using Statistics: median

const DATA_DIRECTORY = joinpath(pkgdir(FissionTemperatureRatio), "data")
const CONFIG_DIRECTORY = joinpath(pkgdir(FissionTemperatureRatio), "config")

# The measured multiplicities, yields and mean kinetic energies: a local copy of public data,
# retrieved from EXFOR.
const DATA_AVAILABLE =
    all(isdir(joinpath(DATA_DIRECTORY, s, "nu_vs_A")) for s in ("Cf252_sf", "U233_nth")) &&
    isdir(joinpath(DATA_DIRECTORY, "U233_nth", "TKE_vs_A"))

DATA_AVAILABLE || @warn "input data not present; the tests that need it are skipped" directory =
    DATA_DIRECTORY

const TEST_MASSES = read_mass_excess_table(String(AME2020_MASS_EXCESS_FILE))
const BSFG_MODEL = BackShiftedFermiGas(TEST_MASSES)
const GC_MODEL = GilbertCameron(
    read_shell_correction_table(String(GILBERT_CAMERON_SHELL_CORRECTION_FILE))
)

"The systems the package ships a configuration for, by label."
const SHIPPED_SYSTEMS = ("Cf252_sf", "U233_nth", "U235_nth", "Pu239_nth")

"The fissioning system of a configuration document, read without its data."
function system_of(document::AbstractDict)
    target = Nuclide(document["system"]["target_Z"], document["system"]["target_A"])
    document["system"]["channel"] == "sf" && return spontaneous_fission(target)
    return neutron_induced_fission(
        target, document["system"]["incident_energy"], document["system"]["channel"]
    )
end

"The fragmentation domain a shipped configuration builds, on Wahl's charge model."
function shipped_domain(label::AbstractString; charges_per_mass::Integer = 5)
    document = TOML.parsefile(joinpath(CONFIG_DIRECTORY, "$(label).toml"))
    system = system_of(document)
    fragmentation = document["fragmentation"]
    return fragmentation_domain(
        system,
        charge_model(TEST_MASSES, system),
        fragmentation["heavy_mass_min"]:fragmentation["heavy_mass_max"];
        charges_per_mass = charges_per_mass,
    )
end

"""
A smooth pre-neutron ⟨TKE⟩(A_H) of the shape the measurements share — a maximum near the doubly
magic heavy fragment, falling to both sides — for tests of the inversion that must not depend on
local data.
"""
function synthetic_kinetic_energy(domain)
    return Dict(
        A => 178.0 - 0.02 * (A - 132)^2 - 8.0 * exp(-((A - 126) / 3)^2) for
        A in domain.heavy_masses
    )
end

"A piecewise-linear reference shaped like the multiplicity ratio: a fall to a minimum, then a rise."
reference_ratio(A_H) = A_H ≤ 130 ? 0.5 - 0.03 * (A_H - 126) : 0.38 + 0.0125 * (A_H - 130)

"""
    variant_configuration(label, directory; set, remove) -> Configuration

A shipped configuration with the keys in `set` (section => key => value) replaced and those in
`remove` (section => keys) dropped, written into `directory` and loaded against the local data.
"""
function variant_configuration(
    label::AbstractString, directory::AbstractString; set = Dict(), remove = Dict()
)
    document = TOML.parsefile(joinpath(CONFIG_DIRECTORY, "$(label).toml"))
    for (section, entries) in set, (key, value) in entries
        document[section][key] = value
    end
    for (section, keys) in remove, key in keys
        delete!(document[section], key)
    end
    path = joinpath(directory, "$(label).toml")
    open(io -> TOML.print(io, document), path, "w")
    return load_configuration(path; data_directory = DATA_DIRECTORY)
end

"""
The settings of the published extraction: the ratio of means over the evaluated Gaussian charge
tables with `ΔZ(A₀/2) = 0`, or over the mean values where no table is held.
"""
function published_settings(label::AbstractString)
    table = joinpath(label, "charge_distribution_vs_A.dat")
    charge = isfile(joinpath(DATA_DIRECTORY, table)) ? table : "mean"
    return (
        set = Dict(
            "fragmentation" => Dict(
                "charge_distribution_file" => charge,
                "zero_polarization_at_symmetry" => true,
            ),
            "level_density" => Dict("ratio_averaging" => "ratio_of_means"),
            "yield" => Dict("symmetrize" => false),
        ),
        remove = Dict("level_density" => ["mean_kinetic_energy_file"]),
    )
end

"The charge tables of the published extraction are held locally for these systems."
const PUBLISHED_TABLES_AVAILABLE = all(
    isfile(joinpath(DATA_DIRECTORY, s, "charge_distribution_vs_A.dat")) for
    s in ("Cf252_sf", "U235_nth", "Pu239_nth")
)

# The total averages ⟨R_T⟩ of Table 1 of Eur. Phys. J. A 60, 190 (2024),
# doi:10.1140/epja/s10050-024-01375-7, back-shifted Fermi gas, five charges per mass number, with
# the relative tolerance each row is held to. The inputs were retrieved independently of whatever
# the authors used, and the breakpoints are selected rather than placed by hand, so agreement is at
# the per-cent level, not to the digits quoted. Fraser (233-U) quotes no uncertainties and has
# fifteen usable pairs; the published value itself carries ±0.25.
const PUBLISHED_TOTAL_AVERAGES = [
    # system, ν(A) dataset, Y(A) distribution, published ⟨R_T⟩, relative tolerance
    ("U233_nth", "K. Nishio 1998", "V.M. Surin 1972", 1.1861, 0.006),
    ("U233_nth", "V.F. Apalin 1965", "V.M. Surin 1972", 1.0346, 0.008),
    ("U233_nth", "J.S. Fraser 1966", "V.M. Surin 1972", 1.3809, 0.015),
    ("Cf252_sf", "C. Budtz-jorgensen 1988", "A. Goeoek 2014", 1.0975, 0.006),
    ("Cf252_sf", "Yu.S. Zamyatnin 1979", "A. Goeoek 2014", 1.1128, 0.006),
    ("Cf252_sf", "A. Goeoek 2014", "A. Goeoek 2014", 1.1163, 0.008),
    ("Cf252_sf", "A. Al-adili 2020", "A. Goeoek 2014", 1.0860, 0.008),
    ("U235_nth", "K. Nishio 1998", "A. Al-adili 2020", 1.1586, 0.008),
    ("U235_nth", "K. Nishio 1998", "Ch.Straede 1987", 1.1644, 0.012),
    ("U235_nth", "A.S. Vorobyev 2010", "A. Al-adili 2020", 1.1186, 0.010),
    ("U235_nth", "A.S. Vorobyev 2010", "Ch.Straede 1987", 1.1221, 0.006),
]

# Widened for the shipped configurations only, which average over Y(A) symmetrized to the
# pre-neutron identity where the published table averaged over it as measured. Symmetrizing
# Göök's and Straede's distributions, the two held on both wings, moves these rows to −0.61 % and
# −1.24 %; the published settings keep the tolerances above.
const SHIPPED_TOLERANCE = Dict(
    ("Cf252_sf", "C. Budtz-jorgensen 1988", "A. Goeoek 2014") => 0.007,
    ("U235_nth", "K. Nishio 1998", "Ch.Straede 1987") => 0.013,
)
