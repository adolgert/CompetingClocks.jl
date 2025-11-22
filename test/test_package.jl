using Test

@safetestset track_package_compat = "Package compat" begin
using TOML

"""
    validate_compat_consistency(root_dir::String)

Recursively searches for `Project.toml` files starting from `root_dir`.
It collects all `[compat]` entries. If a package appears in multiple
Project.toml files, it ensures the version strings are identical.

Returns `(true, [])` if valid, or `(false, messages)` if mismatches are found.
"""
function check_compat_consistency()
    # Dictionary to store: PackageName => (VersionString, RelativePathToFile)
    seen_compat = Dict{String, Tuple{String, String}}()
    errors = String[]
    
    check_root = pwd()
    parent = dirname(check_root)
    while !isfile(joinpath(check_root, "LICENSE"))
        parent = dirname(check_root)
        if parent == check_root
            error("Cannot find repository root")
        end
        check_root = parent
    end
    # Walk through every directory in the repo
    projects = String[]
    for (root, dirs, files) in walkdir(check_root)
        append!(projects, joinpath.(root, filter!(f -> endswith(f, "Project.toml"), files)))
    end

    for file in projects
        data = try
            TOML.parsefile(file)
        catch e
            push!(errors, "Could not parse TOML file at $file: $e")
            continue
        end

        # Check if [compat] exists
        if haskey(data, "compat")
            for (pkg, ver) in data["compat"]
                # If we have seen this package before...
                if haskey(seen_compat, pkg)
                    prev_ver, prev_file = seen_compat[pkg]
                    
                    # COMPARE: Strict string equality
                    if prev_ver != ver
                        msg = """
                        Compat mismatch for package '$pkg':
                            - $prev_ver (in $prev_file)
                            - $ver (in $file)
                        """
                        push!(errors, msg)
                    end
                else
                    # Record it for the first time
                    seen_compat[pkg] = (ver, file)
                end
            end
        end
    end

    return isempty(errors), errors
end
pass, errors = check_compat_consistency()
if !pass
    for error in errors
        print(error)
    end
end
@test pass
end
