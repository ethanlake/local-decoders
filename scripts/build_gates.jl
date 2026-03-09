using LocalDecoders
using ArgParse, JLD2

const DATADIR = joinpath(@__DIR__, "..", "data")

"""
standalone script for building and saving gate circuits.
"""

function main(; model="tc", gate="R", l=0, adj="")
    n = model ∈ ["rep_5bit" "k2"] ? 5 : 3

    fout = joinpath(DATADIR, "$(model)$(l)_gates_$(gate)$(adj).jld2")

    gates, depth, Lx, Ly = master_gate_builder(model,l,gate)

    println("writing gates to file: $fout")
    mkpath(dirname(fout))
    f = jldopen(fout,"w")
    write(f,"model",model); write(f,"gates",gates); write(f,"gate",gate); write(f,"n",n); write(f,"l",l)
    write(f,"depth",depth); write(f,"Lx",Lx); write(f,"Ly",Ly); write(f,"L",Lx)
    close(f)
end

function parse_args()
    s = ArgParseSettings(description="Build and save gate circuits")
    @add_arg_table! s begin
        "--model"
            help = "code model: rep, rep_5bit, k2, tc, twod_rep"
            default = "tc"
        "--gate"
            help = "gate type: R, Z, Y, etc."
            default = "R"
        "--level", "-l"
            help = "level"
            arg_type = Int
            default = 0
        "--adj"
            help = "distinguishing info for output file"
            default = ""
    end
    return ArgParse.parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_args()
    main(; model=args["model"], gate=args["gate"], l=args["level"], adj=args["adj"])
end
