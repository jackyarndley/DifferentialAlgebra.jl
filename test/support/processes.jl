# Capture both streams so failures in loading and threading checks are reviewable.
function fresh_julia(script; threads = 1, timeout = 120)
    project = dirname(Base.active_project())
    command = `$(Base.julia_cmd()) --startup-file=no --threads=$threads`
    command = `$command --project=$project -e $script`
    return mktemp() do _, io
        process = run(pipeline(command; stdout = io, stderr = io); wait = false)
        expired = timedwait(() -> process_exited(process), timeout) === :timed_out
        if expired
            kill(process)
        end
        wait(process)
        flush(io)
        seekstart(io)
        output = read(io, String)
        ok = !expired && success(process)
        if !ok
            @error "Fresh Julia process failed" command expired
            print(stderr, output)
        end
        return ok
    end
end
