# Pre-commit: format against .JuliaFormatter.toml with the version-bounded formatter, then run the tests.
#
#     julia check.jl            format in place, then test
#     julia check.jl --check    fail on formatting differences, do not rewrite
#
# Mirrors the gate that CI applies, so a failure is reproducible locally.

# Directories and root scripts the formatter covers; format.yml applies the same list.
const TARGETS = (
    "src", "test", "scripts", "bench", "docs", "ext", "formatter", "activate.jl", "check.jl"
)

let overwrite = !("--check" in ARGS)
    include(joinpath(@__DIR__, "formatter", "activate.jl"))
    using JuliaFormatter: format

    # One pass of the formatter is not always its fixed point: a construct it rewrites can take a
    # second pass to settle, and CI checks against the fixed point. Format until nothing changes.
    formatted = true
    settled = false
    for pass in 1:(overwrite ? 4 : 1)
        settled = true
        for target in TARGETS
            path = joinpath(@__DIR__, target)
            ispath(path) || continue
            settled &= format(path; overwrite = overwrite)
        end
        pass == 1 && (formatted = settled)
        settled && break
    end

    if !formatted
        if overwrite
            settled || error("formatting did not settle in four passes")
            @info "formatting applied"
        else
            @error "formatting differs from .JuliaFormatter.toml; run `julia check.jl`"
            exit(1)
        end
    end

    include(joinpath(@__DIR__, "activate.jl"))
    Pkg.test()
end
